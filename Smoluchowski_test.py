import numpy as np
import matplotlib.pyplot as plt
import sys


class Smoluchowski():
	"""Simulator using Smoluchowski Eqn"""
	def __init__(self, init_no_particles, sim_length, k = 1, max_agg_size = 5):
		self.init_no_particles = init_no_particles
		self.sim_length = sim_length #number of time steps for simulation
		self.k = float(k)
		
		self.max_agg_size = (max_agg_size + 1)#+1 so we can index from 1 not 0 
		
		self.aggregates = np.zeros(self.max_agg_size)
		self.aggregates[1] = self.init_no_particles #set number of aggregate of size one to initial number of particles 

		self.data = np.zeros((self.sim_length, self.max_agg_size)) #+1 for tau

		self.tau = 2/(self.k * self.init_no_particles)
		self.ttau = [] #list for holding time/tau values  
		
		self.current_no_aggregates = [] #empty lists for keeping track of aggregate and particle numbers over time
		self.current_no_particles = []

		self.aggregate_count = self.init_no_particles #count of aggregates and particles in system is initallised as initial no particles 
		self.particle_count = self.init_no_particles

		self.dt = 0.0001 #0.0001 arbritary, actual value possibly physically related to k
		

	def produce(self,j,i):
		self.aggregates[j] += (self.k/2 * (self.aggregates[i]*self.aggregates[j-i]*self.dt))

	def decay(self,j,i):
		if i + j > self.max_agg_size:
			pass
		else:
			self.aggregates[j] -= (self.k*self.aggregates[j]*self.aggregates[i]*self.dt)

	def get_j_counts_at_step(self,j,t):
		self.data[t,j]=self.aggregates[j]/self.init_no_particles

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
				self.get_j_counts_at_step(j,t)

				for i in range(1,self.max_agg_size):
					self.decay(j,i)
					if j != 1 and i < j:
						self.produce(j,i)
							
				self.count_total_aggregates(j)	
				self.count_total_particles(j)

			self.ttau.append(t/self.tau)
		
			sys.stdout.write(" Simulation progress: %.1f%%   \r" %(t*100/float(self.sim_length)))
			sys.stdout.flush()
			# print t/float(self.sim_length)
	
	def run_analytical(self):
		for t in range(0,self.sim_length):		
			self.add_aggregate_count()
			self.add_particle_count() 

			for j in range(1,self.max_agg_size):			
				self.get_j_counts_at_step(j,t)

				self.aggregates[j] = self.init_no_particles * (t/self.tau)**(j - 1) * (1 + t/self.tau)**(-j-1)	
							
				self.count_total_aggregates(j)	
				self.count_total_particles(j)

			self.ttau.append(t/self.tau)
		
			sys.stdout.write(" Simulation progress: %.1f%%   \r" %(t*100/float(self.sim_length)))
			sys.stdout.flush()
			# print t/float(self.sim_length)
		print self.data
			
	def plot(self):
		for i in range(1,5):
			# print self.data[:,i]
			plt.plot(self.ttau,self.data[:,i], label = "N%i"%(i))
	
		# plt.plot(self.ttau, self.current_no_aggregates, label = "Total no. of aggregates") #no. aggregates as func time
		plt.plot(self.ttau, self.current_no_particles, label = "Total number of particles") #should be horizonta line at 1 but is decaying 
		plt.xlabel("t / tau")
		plt.ylabel("N/N_total")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,x2,y1,1.1))
		# plt.legend()
		plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		plt.show()	

	

class PlotExperimental():

	def __init__(self, file = "dist_aggs_sizes.dat"):
		self.file = str(file) 
		self.read_data()
		
	def read_data(self): #reads in data and uses numpy to make it array like objects
		
		data = np.loadtxt(self.file) 
		
		self.area_x = []
		self.agg_size_y = []
		
		for i in range(len(data)):
			self.area_x.append(data[i,0])
			self.agg_size_y.append(data[i,1])
		print data
		return (self.area_x,self.agg_size_y)
	
	def plot(self):
		plt.plot(self.area_x, self.agg_size_y, label = "experimantal ") #should be horizonta line at 1 but is decaying 
		plt.xlabel("Area")
		plt.ylabel("Aggregate Size")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,100,y1,y2))
		# plt.legend()
		plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		plt.show()	

class PlotAnalytical():
	
	def __init__(self, init_no_particles, sim_length, k = 1, max_agg_size = 5):
		self.init_no_particles = init_no_particles
		self.sim_length = sim_length #number of time steps for simulation
		self.k = float(k)
		
		self.max_agg_size = (max_agg_size + 1)#+1 so we can index from 1 not 0

	def sim(self):
		for t in range(0,self.sim_length):		
			self.add_aggregate_count()
			self.add_particle_count() 

			for k in range(1,self.max_agg_size):			
				# self.get_j_counts_at_step(j,t)

				# for i in range(1,self.max_agg_size):
				# 	self.decay(j,i)
				# 	if j != 1 and i < j:
				# 		self.produce(j,i)
				Pk = self.init_no_particles * (t/self.tau)^(j - 1) * (1 + t/self.tau)^(-j-1)	
							
				self.count_total_aggregates(j)	
				self.count_total_particles(j)

			self.ttau.append(t/self.tau)
		
			sys.stdout.write(" Simulation progress: %.1f%%   \r" %(t*100/float(self.sim_length)))
			sys.stdout.flush()
			# print t/float(self.sim_length)


def main():
	# sim = Smoluchowski(35, 5000, 1, 100)
	sim = Smoluchowski(35, 5000, 1, 10)
	exp = PlotExperimental().plot()
	# sim.run_sim()
	# sim.run_analytical()
	# sim.plot()
	# sim(100, 300, 0.00075, 4).update()
main()
