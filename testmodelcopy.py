import numpy as np
import matplotlib.pyplot as plt

class sim():
	"""Simulator using Smoluchowski Eqn"""
	def __init__(self, no_particles, sim_length, k = 1, max_agg_size = 4):
		self.no_particles = no_particles
		self.sim_length = sim_length 
		self.k = float(k)
		self.max_agg_size = max_agg_size + 1
		self.aggregates = np.zeros(self.max_agg_size)
		self.aggregates[1] = self.no_particles

		self.data = []
		self.total_check = []
		for i in range(len(self.aggregates)):
			self.data.append([])
			self.total_check.append([])



	
	def update(self):

		for dt in range(0,self.sim_length):
			for j in range(1,len(self.aggregates)):
				for i in range(1,self.max_agg_size):
					if j == 1:
						self.aggregates[j] -= (self.k*self.aggregates[j]*self.aggregates[i])
				
					else:
						self.aggregates[j] -= (self.k*self.aggregates[j]*self.aggregates[i])
						if i < j:		
							self.aggregates[j] += (self.k/2 * (self.aggregates[i]*self.aggregates[j-i]*dt))
				
				

				self.data[j].append(self.aggregates[j]/self.no_particles)
				
				self.total_check[j] = self.aggregates[j]*j
				
				print self.total_check[j]
				# print self.aggregates[j]
				# print self.aggregates[j]/self.no_particles
		# print self.total_check
		for i in self.data:
			plt.plot(i)
		plt.show()	
		
					



sim(1000, 50, 0.0005, 4).update()

