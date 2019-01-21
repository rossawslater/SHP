import numpy as np
import matplotlib.pyplot as plt


class Smoluchowski():
	"""Simulator using Smoluchowski Eqn"""
	def __init__(self, init_no_particles, sim_length, k = 1, max_agg_size = 5):
		self.init_no_particles = init_no_particles
		self.sim_length = sim_length #number of time steps for simulation
		self.k = float(k)
		
		self.max_agg_size = (max_agg_size + 1)#+1 so we can index from 1 not 0 
		
		self.aggregates = np.zeros(self.max_agg_size)
		self.aggregates[1] = self.init_no_particles #set number of aggregate of size one to initial number of particles 

		self.data = []

		for i in range(len(self.aggregates)):
			self.data.append([]) #make data have empty list for each aggregate size, to pupulate over time

		self.tau = 2/(self.k * self.init_no_particles)
		self.ttau = [] #list for holding time/tau values  
		
		self.current_no_aggregates = [] #empty lists for keeping track of aggregate and particle numbers over time
		self.current_no_particles = []

		self.aggregate_count = self.init_no_particles #count of aggregates and particles in system is initallised as initial no particles 
		self.particle_count = self.init_no_particles

		self.dt = 0.0001 #1 is arbritary, actual value possibly physically related to k

	def produce(self,j,i):
		self.aggregates[j] += (self.k/2 * (self.aggregates[i]*self.aggregates[j-i]*self.dt))

	def decay(self,j,i):
		if i + j > self.max_agg_size:
			pass
		else:
			self.aggregates[j] -= (self.k*self.aggregates[j]*self.aggregates[i]*self.dt)

	def get_j_counts_at_step(self,j):
		self.data[j].append(self.aggregates[j]/self.init_no_particles)

	def add_aggregate_count(self):
		self.current_no_aggregates.append(self.aggregate_count/self.init_no_particles)
		self.aggregate_count = 0

	def add_particle_count(self):
		self.current_no_particles.append(self.particle_count/self.init_no_particles)
		self.particle_count = 0 

	def	count_total_aggregates(self,j):
		self.aggregate_count += self.aggregates[j]
			
	def count_total_particles(self,j):	
		self.particle_count += self.aggregates[j]*j
	
	def run_sim(self):
		for t in range(0,self.sim_length):		
			self.add_aggregate_count()
			self.add_particle_count() 

			for j in range(1,self.max_agg_size):			
				self.get_j_counts_at_step(j)

				for i in range(1,self.max_agg_size):
					self.decay(j,i)
					if j != 1 and i < j:
						self.produce(j,i)
							
				self.count_total_aggregates(j)	
				self.count_total_particles(j)

			self.ttau.append(t/self.tau)
				
	def plot(self):
		for i in range(1,5):#len(self.data)):
			plt.plot(self.ttau,self.data[i], label = "N%i"%(i))
	
		plt.plot(self.ttau, self.current_no_aggregates, label = "Total no. of aggregates") #no. aggregates as func time
		plt.plot(self.ttau, self.current_no_particles, label = "Total number of particles") #should be horizonta line at 1 but is decaying 
		plt.xlabel("t / tau")
		plt.ylabel("N/N_total")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,x2,y1,1.1))
		plt.legend()
		plt.show()	
	
					


def main():
	sim = Smoluchowski(35, 5000, 1, 100)
	sim.run_sim()
	sim.plot()
	# sim(100, 300, 0.00075, 4).update()
main()
