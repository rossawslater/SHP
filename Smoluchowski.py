import numpy as np
import matplotlib.pyplot as plt
import sys


class Smoluchowski():
	"""Simulator using Smoluchowski Eqn"""
	# def __init__(self, init_no_particles, sim_length, dt = 0.1, max_agg_size = 5, k = 6E-4, mu = 1.023):
	def __init__(self,input):
		#INITIAL PARAMETERS

		self.init_no_particles = int(input[0])
		self.len_mins = int(input[1])
		self.sim_length = int(int(input[1])/float(input[2])) #number of time steps for simulation,
		self.dt = float(input[2])
		self.max_agg_size = (int(input[3]) + 1)#+1 so we can index from 1 not 0
		self.k = float(input[4])# self.k = 3.167E-3 #Literature value of D for Pa
		self.mu = float(input[5])


		# self.current_max_agg_size = 2 #need to adjust length of data array if using this?

		#COUNTERS
		self.aggregates = np.zeros(self.max_agg_size)
		self.aggregates[1] = self.init_no_particles #set number of aggregate of size one to initial number of particles

		self.time = [] #list for holding time

		self.current_no_aggregates = [] #empty lists for keeping track of aggregate and particle numbers over time
		self.current_no_particles = []

		self.aggregate_count = self.init_no_particles #count of aggregates and particles in system is initallised as initial no particles
		self.particle_count = self.init_no_particles

		#DATA ARRAYS FOR OUTPUT
		self.data = np.zeros((self.sim_length, self.max_agg_size))
		self.weighted_data = np.zeros((self.sim_length, self.max_agg_size))
		self.init_ratio_data = np.zeros((self.sim_length, self.max_agg_size))

		# self.times = [103,137,162,201,252,297] #times in gavin's data
		self.times = range(0,100,10)

	def produce(self,j,i):
		self.aggregates[j] += (self.k/2 * (self.aggregates[i]*self.aggregates[j-i]*self.dt))

	def decay(self,j,i):
		if i + j > self.max_agg_size:
			pass
		else:
			self.aggregates[j] -= (self.k*self.aggregates[j]*self.aggregates[i]*self.dt)

	def grow(self):
		self.aggregates[1] += self.mu*self.dt*self.aggregates[1]


	def get_j_counts_at_step(self,j,t):
		# self.data[t,j]=self.aggregates[j]*j/self.init_no_particles
		self.data[t,j]=self.aggregates[j]#/self.init_no_particles

	def get_weighted_j_counts_at_step(self,j,t): #dont ned this as already probablility?
		self.weighted_data[t,j] = (self.aggregates[j]*j)#/self.init_no_particles

	def get_init_ratio_at_step(self,j,t):
		# self.data[t,j]=self.aggregates[j]*j/self.init_no_particles
		self.init_ratio_data[t,j]=self.aggregates[j]/self.init_no_particles

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
			self.grow()
			self.get_aggregate_count()
			self.get_particle_count()

			for j in range(1,self.max_agg_size):
			# for j in range(1,self.current_max_agg_size):
				self.get_j_counts_at_step(j,t)
				self.get_weighted_j_counts_at_step(j,t)
				self.get_init_ratio_at_step(j,t)

				for i in range(1,self.max_agg_size):
				# for i in range(1,self.current_max_agg_size):
					self.decay(j,i)
					if j != 1 and i < j:
						self.produce(j,i)

				self.count_total_aggregates(j)
				self.count_total_particles(j)

			self.time.append(t*self.dt)
			self.data[t,0] = t*self.dt
			self.weighted_data[t,0] = t*self.dt
			self.init_ratio_data[t,0] = t*self.dt

			print t*self.dt

	def run_analytical(self):
		self.tau = 2/(self.k * self.init_no_particles) #t/tau needed to mimic plot from textbook
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
			plt.plot(self.data[:,0],self.data[:,i], label = "N%i"%(i))

		plt.xlabel("t")
		plt.ylabel("N/N Initial")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,x2,y1,y2))
		# plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		plt.show()

		for i in range(1,len(self.init_ratio_data[0])):#self.max_agg_size):
			plt.plot(self.init_ratio_data[:,0],self.init_ratio_data[:,i], label = "N%i"%(i))

		plt.xlabel("t")
		plt.ylabel("N/N Initial")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,x2,y1,y2))
		# plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		plt.show()

		for i in range(1,self.max_agg_size):
			plt.plot(self.weighted_data[:,0],self.weighted_data[:,i], label = "N%i"%(i))

		plt.xlabel("t")
		plt.ylabel("N/N Initial")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		# plt.axis((x1,x2,y1,y2))
		# plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		plt.show()

	def plot_hists(self, data):
		x = 0
		for i in self.times:
			t = i/(self.dt)
			# print i, t, t*self.dt
			plt.plot(data[t,:])#, label = "N%i"%(i))
			x1,x2,y1,y2 = plt.axis()
			plt.axis((1,len(data[0,:]),y1,20000))
			plt.title("Distribution of Aggregates at %i Minutes"%(i))
			plt.ylabel("Number of Bacteria in Agg of size N")
			plt.xlabel("Aggregate Size N")
			# plt.savefig("Growth_Distribution_at_%i_mins.png"%(i))
			plt.show()
			# print np.sum(self.data[(i*10),1:]) #total probablility not conserved
			# x += np.sum(self.data[(i*10),1:])
		# print x

	def load_data(self,file):
		self.data = np.loadtxt(file, skiprows = 1)
		self.params = np.loadtxt(file)[0]

	def load_weighted_data(self,file):
		self.weighted_data = np.loadtxt(file, skiprows = 1)
		self.params = np.loadtxt(file)[0]

	def get_header(self):
		self.header = str(str(self.init_no_particles) + " " + str(self.len_mins) + " " + str(self.dt) + " " + str(self.max_agg_size-1) + " " + str(self.k) + " " + str(self.mu))

	def save(self):
		self.get_header()
		np.savetxt(str(self.header + " growth.txt"), self.data, header = self.header, comments='')
		np.savetxt(str(self.header + " growth weighted.txt"), self.weighted_data, header = self.header, comments='')

def main():
	sim = Smoluchowski(sys.argv[1:])
	# #LOAD PREVIOUS DATA
	sim.load_data("20000 300 0.05 100 0.0006 1.023.txt")
	sim.load_weighted_data("20000 300 0.05 100 0.0006 1.023 weighted.txt")

	#ANALYTICAL SOLUTION
	# sim.run_analytical()

	#RUN SIM AND OUTPUT
	# sim.run_sim()
	# sim.save()

	#PLOT
	# sim.plot()
	sim.plot_hists(sim.weighted_data)
main()
