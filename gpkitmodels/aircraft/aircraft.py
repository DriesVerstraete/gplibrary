# coding=utf-8
"Implements an aircraft model composed of multiple sub-models"
from gpkit import Model, Variable, SignomialsEnabled, LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from vtail import VerticalTail
from fuselage import Fuselage
from landing_gear import LandingGear
from htail import HorizontalTail
from wing import Wing
from wingbox import WingBox
from gpkit.tools import te_exp_minus1

class Aircraft(Model):
    """
    Combined fuselage, tail, and landing gear model
    """

    def __init__(self):

        D      = Variable('D', 'N', 'Total aircraft drag (cruise)')
        Dfuse  = Variable('D_{fuse}', 'N', 'Fuselage drag')
        Dht    = Variable('D_{ht}', 'N', 'Horizontal tail drag')
        Dvt    = Variable('D_{vt}', 'N', 'Vertical tail drag')
        Dwing  = Variable('D_{wing}', 'N', 'Wing drag')
        LD     = Variable('\\frac{L}{D}', '-', 'Lift/drag ratio')
        Lh     = Variable('L_h', 'N', 'Horizontal tail downforce')
        Lw     = Variable('L_w', 'N', 'Wing lift')
        R      = Variable('Range', 'nautical_miles', 'Range')
        ct     = Variable('c_t', '1/s',
                          'Thrust specific fuel consumption')
        V      = Variable('V_{\\infty}', 'm/s', 'Cruise velocity')
        W      = Variable('W', 'N', 'Total aircraft weight')
        We     = Variable('W_{empty}', 'N', 'Empty aircraft weight')
        Weng   = Variable('W_{eng}', 'N', 'Engine weight')
        Wfuse  = Variable('W_{fuse}', 'N', 'Fuselage weight')
        Wfuel  = Variable('W_{fuel}', 'N', 'Fuel weight')
        Wht    = Variable('W_{ht}', 'N', 'Horizontal tail weight')
        Wlg    = Variable('W_{lg}', 'N', 'Landing gear weight')
        Wvt    = Variable('W_{vt}', 'N', 'Vertical tail weight')
        Wwing  = Variable('W_{wing}', 'N', 'Wing weight')
        Zbre   = Variable('Z_{bre}', 'm**0.5*s/kg**0.5', 'Breguet parameter')
        g      = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
        xCGeng = Variable('x_{CG_{eng}}', 'm', 'x-location of engine CG')
        xCGfu  = Variable('x_{CG_{fu}}', 'm', 'x-location of fuselage CG')
        xCGht  = Variable('x_{CG_{ht}}', 'm', 'x-location of htail CG')
        xCGlg  = Variable('x_{CG_{lg}}', 'm', 'x-location of landing gear CG')
        xCGvt  = Variable('x_{CG_{vt}}', 'm', 'x-location of vtail CG') 
        xCGwing = Variable('x_{CG_{wing}}', 'm', 'x-location of wing CG')
        xw     = Variable('x_w', 'm', 'x-location of wing aerodynamic center')
        rho    = Variable('\\rho', 'kg/m^3', 'Air density')
        Sw     = Variable('S_w', 'm**2', 'Wing reference area')
        CL     = Variable('C_{L_w}', '-', 'Lift coefficient')
        CD     = Variable('C_D', '-', 'Drag coefficient')


        with SignomialsEnabled():

            objective = Wfuel
            hlc = [# High level constraints
                   D >= Dvt + Dfuse       + Dwing + Dht,
                   We >= Wvt + Wfuse + Wlg + Wwing + Wht + Weng,
                   TCS([Wfuel >= W - We]), # [SP]

                   # Range equation for a jet
                   TCS([R + Zbre*We**0.5 <= Zbre*W**0.5]),
                   Zbre == 2*(2/(rho*Sw))**0.5*(1/ct)*CL**0.5/CD, 
                   D == 0.5*rho*V**2*Sw*CD,
                   W == 0.5*rho*V**2*Sw*CL,
                   Lw >= W, # TODO: add Lh

                   # CG relationships
                   TCS([xCG*W >= Wvt*xCGvt + Wfuse*xCGfu + Wlg*xCGlg
                               + Wwing*xCGwing + Wht*xCGht + Weng*xCGeng],
                       reltol=1E-2, raiseerror=False),
                   TCS([0.99*xCG*W <= Wvt*xCGvt + Wfuse*xCGfu + Wlg*xCGlg
                               + Wwing*xCGwing + Wht*xCGht + Weng*xCGeng],
                       reltol=1E-2, raiseerror=False),
                   xw == xCGwing,
                   xCGeng == xCGwing,
                  ]

            # Subsystem models
            vt = VerticalTail.aircraft_737()
            vts = VerticalTail.standalone_737()
            fu = Fuselage.aircraft_737()
            fus = Fuselage.standalone_737()
            lg = LandingGear.aircraft_737()
            lgs = LandingGear.standalone_737()
            ht = HorizontalTail.aircraft_737()
            hts = HorizontalTail.standalone_737()
            wi = Wing.aircraft_737()
            wb = WingBox()

            # Need to initialize solve with solution of uncoupled models
            vt_sol = vts.localsolve(verbosity=0)
            fu_sol = fus.localsolve(verbosity=0)
            lg_sol = lgs.localsolve(verbosity=0)
            ht_sol = hts.localsolve(verbosity=0)

            init = vt_sol['variables'].copy()
            init.update(fu_sol['variables'])
            init.update(lg_sol['variables'])
            init.update(ht_sol['variables'])
            init.update({
                         'x_{CG}': 15,
                         'x_{CG_{fu}}': 15,
                         'x_{CG_{ht}}': 38,
                         'x_{CG_{lg}}': 16,
                         'x_{CG_{vt}}': 35,
                        })

            self.init = init

        substitutions = {
                         'Range': 3000,
                         'c_t': 9E-5, # [1/s] CF6 (wiki)
                         'V_{\\infty}': 234,
                         'W_{eng}': 10000,
                        }

        lc = LinkedConstraintSet([hlc, vt, fu, lg, ht, wi],
                                 exclude=[vk.name for vk in wb.varkeys])
        Model.__init__(self, objective,
                             lc,
                             substitutions)

    def test(self):
        sol = self.localsolve(x0=self.init)
        return sol

if __name__ == "__main__":
    a = Aircraft()
    sol = a.test()

