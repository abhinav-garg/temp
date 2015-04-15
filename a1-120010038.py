from math import pi
# from matplotlib import subplots

from pylab import meshgrid, quiver, figure, show, axis
from numpy import arange, zeros
import code
from pandas import DataFrame
import numpy as np

import time
import matplotlib.pyplot as plt

# Potential flow (incompressible, inviscid, irrotational) 
# Potential, stream function, and velocities add linearly
# Analysis in 2 dimensions

class Vortex(object):
	def __init__(self, strength=1, location=(0.0,0.0), vel=(0.0,0.0)):
		self.location = location
		self.strength = strength
		self.velocity = vel

	def velocityAt(self, location):
		# Evaluate value at a given location
		if location == self.location:
			return (0,0)
		func = -1j*(self.strength)/(2*pi*(complex(location[0], location[1])-complex(self.location[0], self.location[1])))
		return (func.real, -func.imag)

	def setVelocity(self, vel):
		# Holds the velocity at the object's location
		self.velocity = vel

	def setIntermediatePosition(self, position):
		# Required to hold an intermediate position sometimes
		self.intermediatePosition = position

	def move(self):
		# Update the current location with intermediate position
		self.location = self.intermediatePosition

class Tracer(object):
	def __init__(self, location=(0.0,0.0), vel=(0.0,0.0)):
		self.location = location
		self.velocity = vel

	def setVelocity(self, vel):
		# Holds the velocity at the object's location
		self.velocity = vel

	def setIntermediatePosition(self, position):
		# Required to hold an intermediate position sometimes
		self.intermediatePosition = position

	def move(self):
		# Update the current location with intermediate position
		self.location = self.intermediatePosition

def euler(vortices, points, dt=0.1):
	(u,v) = (0.0,0.0)
	t = 0
	
	# Current state -> Increment by dt
	for index, point in enumerate(points):
		# code.interact(local=dict(globals(), **locals()))	
		for vor in vortices:
			vel = vor.velocityAt(point.location)
			u += vel[0]
			v += vel[1]
			# code.interact(local=dict(globals(), **locals()))
		point.setIntermediatePosition((point.location[0]+u*dt, point.location[1]+v*dt))
		t += dt
		# Reset velocity for evaluation at next point
		(u,v) = (0.0,0.0)
		# code.interact(local=dict(globals(), **locals()))
	
	# Update all positions simultaneously - thereby ensuring no conflicts
	for point in points:
		point.move()
	

def rk2(vortices, points, dt=0.1):
	(u,v) = (0.0,0.0)
	t = 0
	
	# Current state -> Increment by dt
	for index, point in enumerate(points):
		# code.interact(local=dict(globals(), **locals()))	
		for vor in vortices:
			vel = vor.velocityAt(point.location)
			u += vel[0]
			v += vel[1]
			# code.interact(local=dict(globals(), **locals()))
		point.setVelocity((u,v))
		point.setIntermediatePosition((point.location[0]+u*dt, point.location[1]+v*dt))
		# Reset velocity for evaluation at next point
		(u,v) = (0.0,0.0)

	# code.interact(local=dict(globals(), **locals()))	

	# Need shifted vortices for intermediate evaluation
	tempVortices = [None]*len(vortices)
	for index, tempVor in enumerate(tempVortices):
		vor = vortices[index]
		tempVortices[index] = Vortex(location=vor.intermediatePosition, vel=vor.velocity)

	# Evaluate for next increment
	# dt, dy = v(tn,yn)dt
	# yn+1 = yn + 0.5( v(tn, yn)dt + v(tn+dt, yn+dy)dt )
	(u,v) = (0.0,0.0)
	for index, point in enumerate(points):
		# code.interact(local=dict(globals(), **locals()))	
		for tempVor in tempVortices:
			vel = tempVor.velocityAt(point.intermediatePosition)
			u += vel[0]
			v += vel[1]
		point.setIntermediatePosition((
			0.5 * (point.intermediatePosition[0] + point.location[0] + dt * u)
			,
			0.5 * (point.intermediatePosition[1] + point.location[1] + dt * v)
			))
		# code.interact(local=dict(globals(), **locals()))	
		# Reset velocity for evaluation at next point
		(u,v) = (0.0,0.0)

	t += dt

	# Update all positions - thereby ensuring no conflicts
	for point in points:
		point.move()
	
def main():
	# Sets interative mode on
	plt.ion()
	v1 = Vortex()
	v2 = Vortex(location=(1.0,0))
	
	t1 = Tracer(location=(0.5,0.5))

	v3 = Vortex(location=(-0.5,0.0))
	v4 = Vortex(location=(0.5,0.0))
	v5 = Vortex(location=(0.0,0.5))

	t = 0

	N = 10
	vortices = [Vortex( location=(0.0,-1 + n/N + 1/2N) ) for n in range(0,N)]
	for v_index, vortex in enumerate(vortices):
		plt.plot(vortex.location[0], vortex.location[1], color="#" + str(010101*v_index))

	code.interact(local=locals())

	# Euler
	# while True:
	# 	t+=0.1
	# 	# rk2([v1, v2], [v1, v2, t1])
	# 	# euler([v1, v2], [v1, v2, t1])
	# 	euler([v1,v2], [v1,v2])
	# 	print t
	# 	print v1.location 
	# 	print v2.location
	# 	plt.plot(v1.location[0], v1.location[1], 'ro')
	# 	plt.plot(v2.location[0], v2.location[1], 'bo')
	# 	plt.axis([-4,4,-4,4])
	# 	plt.show()
	# 	code.interact(local=dict(globals(), **locals()))

	# RK2
	while True:
		t+=0.1
		# rk2([v1, v2], [v1, v2, t1])
		# euler([v1, v2], [v1, v2, t1])
		rk2([v1,v2], [v1,v2])
		print t
		print v1.location 
		print v2.location
		plt.plot(v1.location[0], v1.location[1], 'ro')
		plt.plot(v2.location[0], v2.location[1], 'bo')
		plt.axis([-4,4,-4,4])
		plt.show()
		code.interact(local=dict(globals(), **locals()))

	# 2 vortices, 10 tracers
	while True:
		t+=0.1
		# rk2([v1, v2], [v1, v2, t1])
		# euler([v1, v2], [v1, v2, t1])
		tracers = [Tracer(location=(0.5,i*0.2)) for i in range(0,10)]
		rk2([v1,v2], [v1,v2]+tracers)
		print t
		print v1.location 
		print v2.location
		plt.plot(v1.location[0], v1.location[1], 'ro')
		plt.plot(v2.location[0], v2.location[1], 'bo')
		x = [tracer.location[0] for tracer in tracers]
		y = [tracer.location[1] for tracer in tracers]
		plt.plot(x, y, 'yo')
		plt.axis([-4,4,-4,4])
		plt.show()
		code.interact(local=dict(globals(), **locals()))

	# 3 vortices, 10 tracers
	while True:
		t+=0.1
		# rk2([v1, v2], [v1, v2, t1])
		# euler([v1, v2], [v1, v2, t1])
		tracers = [Tracer(location=(i*0.2,i*0.2)) for i in range(0,10)]
		rk2([v3,v4,v5], [v3,v4,v5]+tracers)
		print t
		print v3.location 
		print v4.location
		print v5.location
		plt.plot(v3.location[0], v3.location[1], 'ro')
		plt.plot(v4.location[0], v4.location[1], 'bo')
		plt.plot(v5.location[0], v5.location[1], 'go')
		x = [tracer.location[0] for tracer in tracers]
		y = [tracer.location[1] for tracer in tracers]
		plt.plot(x, y, 'yo')
		plt.axis([-4,4,-4,4])
		plt.show()
		code.interact(local=dict(globals(), **locals()))

	code.interact(local=dict(globals(), **locals()))

def init():
	line.set_data([],[])
	return line,

if __name__ == "__main__":
	main()








