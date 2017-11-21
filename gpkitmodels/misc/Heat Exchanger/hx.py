from gpkit import Model, Vectorize, SignomialsEnabled, SignomialEquality,
from gpkit import MVariable, VarKey, units
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

	def dynamic(self,state):
		return ChannelP(self,state)

class ChannelP(Model):
	"""single channel HX performance model

	Variables
	---------
	cp     
	eps
	Hdot

	"""
	def setup(self):
		exec parse_variables(ChannelP.__doc__)

if __name__ == "__main__":
	m = Channel()
	m.substitutions.update({})