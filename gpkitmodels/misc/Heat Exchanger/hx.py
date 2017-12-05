from gpkit import Model, Vectorize, SignomialsEnabled, SignomialEquality
from gpkit import Variable, VarKey, units, parse_variables
from gpkit.constraints.bounded import Bounded
import numpy as np
import matplotlib.pyplot as plt


class Channel(Model):
    """single channel HX model

    Variables
    ---------
    A_r                 [m^2]     frontal area
    h                   [m]       channel height
    l                   [m]       channel length
    eta_h_ref  0.917    [-]       reference effectiveness
    Pf_ref  21.66       [-]       reference pressure drop parameter
    Re_ref  90550       [-]       reference Reynolds number
    Upper Unbounded
    ---------------
    A_r, l

    Lower Unbounded
    ---------------
    A_r, l

    """
    def setup(self):
        exec parse_variables(Channel.__doc__)
        return [A_r == h**2]

    def dynamic(self,state):
        return ChannelP(self,state)

class ChannelP(Model):
    """ single channel air HX performance model

    Variables
    ---------
    alpha             [-]         T_out/T_in
    cp         1004   [J/(kg*K)]  heat capacity of air
    dT                [-]         wall/free-stream temp ratio - 1
    eta_h             [-]         effectiveness
    fV                [-]         air velocity ratio across channel
    fr                [N/m^2]     force per frontal area
    Hdot       9.9    [J/s]       heat flow rate
    mdot              [kg/s]      air mass flow rate
    Pf                [-]         pressure drop parameter
    Re                [-]         Reynolds number
    Tr         85+273 [K]         wall temperature 
    A_e               [m^2]       flow exit area


    Upper Unbounded
    ---------------
    fV, A_e, P0_in, mu_out, mu_in

    Lower Unbounded
    ---------------
    fr, fV, V_out, P_out, V_in, rho_in, rho_out, mu_out, mu_in, Pf_ref
 
    """
    def setup(self, channel, state):
        self.channel = channel
        exec parse_variables(ChannelP.__doc__)

        mu_in = self.mu_in = state.mu_in
        mu_out = self.mu_out = state.mu_out
        rho_in = self.rho_in = state.rho_in
        rho_out = self.rho_out = state.rho_out
        V_in = self.V_in = state.V_in
        V_out = self.V_out = state.V_out
        T_in = self.T_in = state.T_in
        T_out = self.T_out = state.T_out
        P0_in = self.P0_in = state.P0_in
        P_out = self.P_out = state.P_out
        A_r = self.channel.A_r
        eta_h_ref = self.eta_h_ref = self.channel.eta_h_ref
        Pf_ref = self.Pf_ref = self.channel.Pf_ref
        Re_ref = self.channel.Re_ref
        l = self.l = self.channel.l

        Pf_rat = Pf/Pf_ref
        Re_rat = Re/Re_ref
        eta_h_rat = eta_h/eta_h_ref

        constraints = []
        constraints += [mdot == rho_in*V_in*A_r,
                        mdot == rho_out*V_out*A_e,
                        Hdot == mdot*cp*dT*eta_h*T_in,
                        fV == V_out/V_in,
                        #channel Reynolds number (geometric average)
                        Re == (rho_in*rho_out*V_in*V_out/mu_in/mu_out)**(0.5)*self.channel.l,
                        #effectiveness
                        eta_h/eta_h_ref == 0.799*Re_rat**-0.0296,
                        eta_h <= 0.844, #to make sure fit is valid, since the curve is signomial
                        #pressure drop
                        Pf_rat**0.155 >= 0.475 * Re_rat ** 0.00121 + 0.0338 * Re_rat ** -0.336,
                        Pf*(0.5*rho_in*V_in**2) == fr,
                        alpha >= 1 + dT*eta_h,
                        dT*T_in + T_in <= Tr,
                        alpha == T_out/T_in,
                        P0_in >= P_out + 0.5*rho_in*V_in**2*Pf + 0.5*rho_out*V_out**2,
                        T_out <= Tr,
                        ]
        return constraints

class HXState(Model):
    """ HX flow state model

    Variables
    ---------
    mu_in      1.7*10**-5             [kg/(m*s)]  incoming dynamic viscosity
    mu_out     1.7*10**-5             [kg/(m*s)]  exit dynamic viscosity
    P0_in      0.5*1.25*30**2+101000 [Pa]        incoming total pressure
    P_in                             [Pa]        incoming static pressure
    P_out      90000                 [Pa]        exit static pressure
    R          287.1                 [J/(kg*K)]  specific gas constant of air 
    rho_in                           [kg/(m^3)]  incoming air density
    rho_out                          [kg/(m^3)]  exiting air density
    V_in       30                    [m/s]       incoming air velocity
    V_out                            [m/s]       exiting air velocity
    T_in       -20+273               [K]         incoming air temperature
    T_out                            [K]         exiting air temperature


    Upper Unbounded
    ---------------
    V_out, T_out

    Lower Unbounded
    ---------------
    rho_out, rho_in, P_in

    """
    #calc_p0in = lambda self, c: c[self.P_in] + 0.5*c[self.rho_in] * c[self.V_in]**2

    def setup(self):
        exec parse_variables(HXState.__doc__)
        return [P0_in >= P_in + 1/2*rho_in*V_in**2,
                V_out == V_in*T_out/T_in,
                P_out == rho_out*R*T_out,
                P_in  == rho_in*R*T_in,
                T_out >= T_in]

class HX(Model):
    """ Heat eXchanger model 

    Upper Unbounded
    ---------------
    mdot, A_e, A_r

    Lower Unbounded
    ---------------
    fr, rho_in, P_in
    """
    def setup(self,state):
        self.channel = Channel()
        self.channelP = self.channel.dynamic(state)
        self.state = state
        self.mdot = self.channelP.mdot
        self.fr = self.channelP.fr
        self.A_e = self.channelP.A_e
        self.A_r = self.channel.A_r
        self.rho_in = self.state.rho_in
        self.P_in = self.state.P_in
        return self.channel, self.channelP, state

if __name__ == "__main__":
    state = HXState()
    m = HX(state)

    m.cost = m.channel.A_r*m.channelP.A_e*m.channel.l*m.channelP.fr
    sol = m.solve()