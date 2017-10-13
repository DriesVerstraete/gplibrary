" vertical tail "
from gpkit import Model, Variable
from .tail_aero import TailAero
from gpkitmodels.GP.aircraft.wing.wing import AeroSurf
from gpkitmodels.GP.aircraft.wing.constant_taper_chord import c_bar
from gpkitmodels.GP.aircraft.wing.wing_interior import WingInterior
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSkin

#pylint: disable=invalid-name, too-many-locals, unused-variable
#pylint: disable=attribute-defined-outside-init

class VerticalTail(Model):
    "vertical tail model"
    def setup(self, N=3, lam=0.8):

        Vv = Variable("V_v", "-", "vertical tail volume coefficient")
        W = Variable("W", "lbf", "vertical tail weight")
        lv = Variable("l_v", "ft", "vertical tail moment arm")
        mfac = Variable("m_{fac}", 1.1, "-", "vertical tail margin factor")

        cb, eta, deta, cbarmac = c_bar(lam, N)
        subdict = {"\\lambda": lam, "\\bar{c}": cb, "\\eta": eta,
                   "\\bar{c}_{ave}": (cb[1:]+cb[:-1])/2, "\\tau": 0.08,
                   "\\bar{c}_{MAC}": cbarmac, "d\\eta": deta,
                   "C_{L_{max}}": 1.5}

        self.surf = AeroSurf(N=N)
        self.surf.substitutions.update(subdict)

        self.skin = WingSkin(self.surf)
        self.skin.substitutions.update({"\\rho_{CFRP}": 0.049})
        self.foam = WingInterior(self.surf)
        self.foam.substitutions.update({"\\bar{A}_{jh01}": 0.0548})
        self.foam.substitutions.update({"\\rho_{foam}": 0.024})

        self.components = [self.skin, self.foam]

        constraints = [W/mfac >= sum([c["W"] for c in self.components])]

        self.flight_model = TailAero

        return constraints, self.surf, self.components
