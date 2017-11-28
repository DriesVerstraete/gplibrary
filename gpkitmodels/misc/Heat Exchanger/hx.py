from gpkit import Model, Vectorize, SignomialsEnabled, SignomialEquality
from gpkit import Variable, VarKey, units, parse_variables
from gpkit.constraints.bounded import Bounded
import numpy as np
import matplotlib.pyplot as plt


class Channel(Model):
    """single channel HX model

    Variables
    ---------
    A_r        [m^2]     frontal area

    Upper Unbounded
    ---------------
    A_r

    Lower Unbounded
    ---------------
    A_r

    """
    def setup(self):
        exec parse_variables(Channel.__doc__)

    def dynamic(self,state):
        return ChannelP(self,state)

class ChannelP(Model):
    """ single channel air HX performance model

    Variables
    ---------
    alpha      [-]        T_out/T_in
    cp         [J/(kg*K)] heat capacity of air
    dT         [-]        wall/free-stream temp ratio - 1
    dTr        [K]        air/wall temperature difference
    fV         [-]        air velocity ratio across channel
    fr         [N/m^2]    force per frontal area
    eps        [-]        effectiveness
    Hdot       [J/s]      heat flow rate
    mdot       [kg/s]     air mass flow rate
    Pf         [-]        pressure drop parameter
    Tr         [K]        wall temperature 
    A_e         [m^2]      flow exit area


    Upper Unbounded
    ---------------
    Hdot, mdot, rho_out, fV, Tr, A_e, P0_in, dTr, V_out, cp

    Lower Unbounded
    ---------------
    Hdot, fr, cp, fV, V_out, P_out, mdot, eps, V_in, A_e, rho_in, Pf, dTr, dT, rho_out
 
    """
    def setup(self,channel,state):
        self.channel = channel
        exec parse_variables(ChannelP.__doc__)

        rho_in = self.rho_in = state.rho_in
        rho_out = self.rho_out = state.rho_out
        V_in = self.V_in = state.V_in
        V_out = self.V_out = state.V_out
        T_in = self.T_in = state.T_in
        T_out = self.T_out = state.T_out
        P0_in = self.P0_in = state.P0_in
        P_out = self.P_out = state.P_out
        A_r = self.channel.A_r

        constraints = []
        constraints += [mdot == state.rho_in*state.V_in*A_r,
                        mdot == rho_out*state.V_out*A_e,
                        Hdot == mdot*cp*dT*eps*state.T_in,
                        #pressure drop
                        Pf*(0.5*state.rho_in*state.V_in**2) == fr,
                        state.T_out/state.T_in >= 1 + dT*eps,
                        dT*state.T_in + state.T_in <= Tr,
                        alpha == state.T_out/state.T_in,
                        state.P0_in >= state.P_out + 0.5*state.rho_in*state.V_in**2*(Pf+alpha**2),
                        state.T_out <= Tr,
                        ]
        return constraints

class HXState(Model):
    """ HX flow state model

    Variables
    ---------
    P0_in      self.calc_p0in [Pa]        incoming total pressure
    P_in       101000         [Pa]        incoming static pressure
    P_out      101000         [Pa]        exit static pressure
    rho_in     1.25           [kg/(m^3)]  incoming air density
    rho_out                   [kg/(m^3)]  exiting air density
    V_in       2              [m/s]       incoming air velocity
    V_out                     [m/s]       exiting air velocity
    T_in       268            [K]         incoming air temperature
    T_out                     [K]         exiting air temperature


    Upper Unbounded
    ---------------
    V_out, T_out

    """
    calc_p0in = lambda self, c: c[self.P_in] + 0.5*c[self.rho_in] * c[self.V_in]**2

    def setup(self):
        exec parse_variables(HXState.__doc__)
        return [V_out == V_in*T_out/T_in,
                V_out == V_in*rho_in/rho_out,
                T_out >= T_in]

class HX(Model):
    """ Heat eXchanger model
    SKIP VERIFICATION
    
    """
    def setup(self,state):
        self.channel = Channel()
        self.channelP = self.channel.dynamic(state)
        self.state = state
        constraints = []
        return constraints, self.channel, self.channelP, state

if __name__ == "__main__":
    state = HXState()
    state.substitutions.update({

        })
    m = HX(state)
    m.substitutions.update({
        m.channelP.Hdot:9.9*units('W'),
        m.channelP.mdot: 0.35*units('kg/s'),
        m.channelP.Tr: (85+273)*units('K'),
        m.channelP.eps: 0.5,
        m.channelP.Pf: 5,
        m.channelP.cp: 1004*units('J/(kg*K)'),
        })

    m.cost = m.channel.A_r + m.channelP.A_e + m.channelP.fr*units('m^4/N')
    sol = m.solve()