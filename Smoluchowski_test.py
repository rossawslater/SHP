import numpy as np
import matplotlib.pyplot as plt
import sys


class Smoluchowski():
	"""Simulator using Smoluchowski Eqn"""
	def __init__(self, init_no_particles, sim_length, max_agg_size = 5, k = 6E-4):
		self.init_no_particles = init_no_particles
		self.sim_length = sim_length #number of time steps for simulation
		# self.k = 3.167E-3 #Literature value of D for Pa
		self.k = k
		self.max_agg_size = (max_agg_size + 1)#+1 so we can index from 1 not 0
		# self.current_max_agg_size = 2 #need to adjust length of data array if using this?
		self.aggregates = np.zeros(self.max_agg_size)
		self.aggregates[1] = self.init_no_particles #set number of aggregate of size one to initial number of particles

		self.data = np.zeros((self.sim_length, self.max_agg_size)) #+1 for tau

		self.tau = 2/(self.k * self.init_no_particles) #t/tau needed to mimic plot from textbook
		self.time = [] #list for holding time

		self.current_no_aggregates = [] #empty lists for keeping track of aggregate and particle numbers over time
		self.current_no_particles = []

		self.aggregate_count = self.init_no_particles #count of aggregates and particles in system is initallised as initial no particles
		self.particle_count = self.init_no_particles

		# self.dt = 0.1
		self.dt = 0.05

		self.weighted_data = np.zeros((self.sim_length, self.max_agg_size))

		# self.gavins_times = [103,137,162,201,252,297]
		self.gavins_times = range(0,25,1)

	def produce(self,j,i):
		self.aggregates[j] += (self.k/2 * (self.aggregates[i]*self.aggregates[j-i]*self.dt))

	def decay(self,j,i):
		if i + j > self.max_agg_size:
			pass
		else:
			self.aggregates[j] -= (self.k*self.aggregates[j]*self.aggregates[i]*self.dt)

	def grow(self):
		self.aggregates[1] += 1.023*self.dt*self.aggregates[1]*10
		# self.aggregates[1] = self.init_no_particles

	def get_j_counts_at_step(self,j,t):
		# self.data[t,j]=self.aggregates[j]*j/self.init_no_particles
		self.data[t,j]=self.aggregates[j]/self.init_no_particles

	def get_weighted_j_counts_at_step(self,j,t): #dont ned this as already probablility?
		self.weighted_data[t,j] = (self.aggregates[j]*j)#/self.init_no_particles

	def get_aggregate_count(self):
		self.current_no_aggregates.append(self.aggregate_count/self.init_no_particles)
		self.aggregate_count = 0

	def get_particle_count(self):
		self.current_no_particles.append(self.particle_count/self.init_no_particles)
		self.particle_count = 0

	def	count_total_aggregates(self,j):
		self.aggregate_count += self.aggregates[j]

	def count_total_particles(self,j):
		self.particle_count += self.aggregates[j]*j

	def run_sim(self):
		# self.current_max_agg_size = self.max_agg_size
		for t in range(0,self.sim_length):
			# self.grow()
			self.get_aggregate_count()
			self.get_particle_count()

			for j in range(1,self.max_agg_size):
			# for j in range(1,self.current_max_agg_size):
				self.get_j_counts_at_step(j,t)
				self.get_weighted_j_counts_at_step(j,t)

				for i in range(1,self.max_agg_size):
				# for i in range(1,self.current_max_agg_size):
					self.decay(j,i)
					if j != 1 and i < j:
						self.produce(j,i)

				self.count_total_aggregates(j)
				self.count_total_particles(j)

			self.time.append(t*self.dt)
			self.data[t,0] = t*self.dt
			print t*self.dt

	def run_analytical(self):
		for t in range(0,self.sim_length):
			dt = self.dt * t
			for j in range(1,self.max_agg_size):

				self.aggregates[j] = self.init_no_particles * (dt/self.tau)**(j-1) * (1 + (dt/self.tau))**(-j-1)

				self.get_j_counts_at_step(j,t)

			self.time.append(t/self.tau)

		for i in self.data:
			print(i[1]/self.init_no_particles)

	def plot(self):
		for i in range(1,len(self.data[0])):#self.max_agg_size):
			# print(self.data[:,i])
			# print "\n", np.argmax(self.data[:,i])
			print len(self.time), len(self.data[:,i])
			# plt.plot(self.time,self.data[:,i], label = "N%i"%(i))

		plt.xlabel("t")
		plt.ylabel("N/N Initial")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,x2,y1,1.1))
		# plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		plt.show()

		for i in range(1,self.max_agg_size):
			plt.plot(self.weighted_data[:,0],self.weighted_data[:,i], label = "N%i"%(i))
		# plt.plot(self.time, self.current_no_aggregates, label = "Total no. of aggregates") #no. aggregates as func time
		# plt.plot(self.time, self.current_no_particles, label = "Total number of particles") #should be horizonta line at 1 but is decaying

		plt.xlabel("t")
		plt.ylabel("N/N Initial")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		# plt.axis((x1,x2,y1,y2))
		# plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		plt.show()

	def plot_hist(self):
		x = 0
		for i in self.gavins_times:
			# print i
			# print self.data[(i*10)]
			# print self.time[(i*10)]
			# plt.plot(self.data[(i*10),:], label = "N%i"%(i))
			plt.plot(self.weighted_data[(i*10),:], label = "N%i"%(i))
			x1,x2,y1,y2 = plt.axis()
			plt.axis((1,200,y1,2000))
			plt.title("Distribution of Aggregates at %i Minutes"%(i))
			plt.ylabel("Number of Bacteria in Agg of size N")
			plt.xlabel("Aggregate Size N")
			plt.savefig("Growth_Distribution_at_%i_mins.png"%(i))
			plt.show()
			print np.sum(self.data[(i*10),1:]) #total probablility not conserved
			x += np.sum(self.data[(i*10),1:])
		print x

	def load_data(self,file):
		self.data = np.loadtxt(file)

	def load_weighted_data(self,file):
		self.weighted_data = np.loadtxt(file)

def main():
	sim = Smoluchowski(20000, 6000, 10)
	sim.load_data("Reaggregation_data.txt")
	sim.load_weighted_data("Reaggregatoin_weighted_data.txt")
	# sim.run_sim()
	# sim.run_analytical()
	sim.plot()
	sim.plot_hist()
	# np.savetxt("Reaggregation_data.txt", sim.data)
	# np.savetxt("Reaggregatoin_weighted_data.txt", sim.weighted_data)
main()
