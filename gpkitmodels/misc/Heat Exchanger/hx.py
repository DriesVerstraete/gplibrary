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

	"""
	def setup(self):
		exec parse_variables(Channel.__doc__)
		constraints = []
		return constraints

	def dynamic(self,state):
		return ChannelP(self,state)

class ChannelP(Model):
	""" single channel air HX performance model

	Variables
	---------
	cp         [J/(kg*K)] heat capacity of air
	dP         [N/m^2]    pressure drop across channel
	dT         [K]        air temperature increase across channel
	dTr        [K]        air/wall temperature difference
	fV         [-]        air velocity ratio across channel
	fr         [N/m^2]    force per frontal area
	eps        [-]        effectiveness
	Hdot       [J/s]      heat flow rate
	mdot       [kg/s]     air mass flow rate
	Pf     5   [-]        pressure drop parameter
 	Tr         [K]        wall temperature 

	"""
	def setup(self,channel,state):
		self.channel = channel
		exec parse_variables(ChannelP.__doc__)

		constraints = []
		constraints += [mdot == state.rho_in*state.V_in*self.channel['A_r'],
						Hdot == mdot*cp*dT,
						eps == dT/dTr,
						#pressure drop
						Pf == fr/(0.5*state.rho_in*state.V_in**2),
						self.channel.A_r >= Hdot/(state.rho_in*state.V_in*cp*state.T_in) /
												(dT*eps)]
		return constraints

class HXState(Model):
	""" HX flow state model

	Variables
	---------
	rho_in     [kg/(m^3)] incoming air density
	rho_out    [kg/(m^3)] exiting air density
 	V_in       [m/s]      incoming air velocity
 	V_out      [m/s]      exiting air velocity
 	T_in       [K]        incoming air temperature
  	T_out      [K]        exiting air temperature

	"""
	def setup(self):
		exec parse_variables(HXState.__doc__)
		constraints = []
		constraints += [V_out == V_in*T_out/T_in,
					    V_out == V_in*rho_in/rho_out]
		return constraints

class HX(Model):
	def setup(self,state):
		self.channel = Channel()
		self.channelP = self.channel.dynamic(state)
		constraints = []
		return constraints

if __name__ == "__main__":
	state = HXState()
	state.substitutions.update({

		})
	m = HX(state)
	m.substitutions.update({

		})