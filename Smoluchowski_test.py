import numpy as np
import matplotlib.pyplot as plt
import sys


class Smoluchowski():
	"""Simulator using Smoluchowski Eqn"""
	def __init__(self, init_no_particles, sim_length, k = 1, max_agg_size = 5):
		self.init_no_particles = init_no_particles
		self.sim_length = sim_length #number of time steps for simulation
		self.k = float(k)
		self.k = 6E-4
		self.max_agg_size = (max_agg_size + 1)#+1 so we can index from 1 not 0
		self.current_max_agg_size = 2 #need to adjust length of data array if using this?
		self.aggregates = np.zeros(self.max_agg_size)
		self.aggregates[1] = self.init_no_particles #set number of aggregate of size one to initial number of particles

		self.data = np.zeros((self.sim_length, self.max_agg_size)) #+1 for tau

		self.tau = 2/(self.k * self.init_no_particles)
		self.ttau = [] #list for holding time/tau values

		self.current_no_aggregates = [] #empty lists for keeping track of aggregate and particle numbers over time
		self.current_no_particles = []

		self.aggregate_count = self.init_no_particles #count of aggregates and particles in system is initallised as initial no particles
		self.particle_count = self.init_no_particles

		self.dt = 0.0075

		self.weighted_data = np.zeros((self.sim_length, self.max_agg_size))
		print self.k

	def produce(self,j,i):
		self.aggregates[j] += (self.k/2 * (self.aggregates[i]*self.aggregates[j-i]*self.dt))

	def decay(self,j,i):
		if i + j > self.max_agg_size:
			pass
		else:
			self.aggregates[j] -= (self.k*self.aggregates[j]*self.aggregates[i]*self.dt)

	def grow(self):
		self.aggregates[1] += 1.023*self.dt*self.aggregates[1]
		# self.aggregates[1] = self.init_no_particles

	def get_j_counts_at_step(self,j,t):
		self.data[t,j]=self.aggregates[j]*j/self.init_no_particles
		# self.data[t,j]=self.aggregates[j]/self.init_no_particles

	def get_weighted_j_counts_at_step(self,j,t): #dont ned this as already probablility?
		self.weighted_data[t,j] = (self.aggregates[j]*j)/self.init_no_particles

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
		self.current_max_agg_size = self.max_agg_size
		for t in range(0,self.sim_length):
			self.grow()
			self.add_aggregate_count()
			self.add_particle_count()

			# for j in range(1,self.max_agg_size):
			for j in range(1,self.current_max_agg_size):
				self.get_j_counts_at_step(j,t)
				self.get_weighted_j_counts_at_step(j,t)

				# for i in range(1,self.max_agg_size):
				for i in range(1,self.current_max_agg_size):
					self.decay(j,i)
					if j != 1 and i < j:
						self.produce(j,i)


				self.count_total_aggregates(j)
				self.count_total_particles(j)


			self.ttau.append(t*self.dt)
			print t*self.dt
			# sys.stdout.write(" Simulation progress: %.1f%%   \r" %(t*100/float(self.sim_length)))
			# sys.stdout.flush()

	def run_analytical(self):
		for t in range(0,self.sim_length):
			dt = self.dt * t
			for j in range(1,self.max_agg_size):

				self.aggregates[j] = self.init_no_particles * (dt/self.tau)**(j-1) * (1 + (dt/self.tau))**(-j-1)

				self.get_j_counts_at_step(j,t)

			self.ttau.append(t/self.tau)

		for i in self.data:
			print(i[1]/self.init_no_particles)

	def plot(self):
		for i in range(1,self.max_agg_size):
			# print(self.data[:,i])
			print "\n", np.argmax(self.data[:,i])
			plt.plot(self.ttau,self.data[:,i], label = "N%i"%(i))
			# plt.plot(self.ttau,self.weighted_data[:,i], label = "N%i"%(i))
		# plt.plot(self.ttau, self.current_no_aggregates, label = "Total no. of aggregates") #no. aggregates as func time
		# plt.plot(self.ttau, self.current_no_particles, label = "Total number of particles") #should be horizonta line at 1 but is decaying

		plt.xlabel("t")
		plt.ylabel("N/N Initial")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		# plt.axis((x1,x2,y1,1.1))
		# plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		plt.show()

	def find_peaks(self):
		self.peaks = np.zeros(self.max_agg_size)
		for j in range(1,self.max_agg_size):
			self.peaks[j] = np.argmax(self.data[:,j])

		print self.peaks

	def plot_hist(self):

		# x = 0
		for i in range(1,self.max_agg_size):
			plt.plot(self.data[self.peaks[i],:], label = "N%i"%(i))
			x1,x2,y1,y2 = plt.axis()
			plt.axis((1,x2,y1,1))
			plt.title("Peak in aggregates of size %i"%(i))
			plt.ylabel("Probability of Bacteria being in Agg of size N")
			plt.xlabel("Aggregate Size N")
			# plt.savefig("Peak at N%i.pdf"%(i))
			plt.show()
			print np.sum(self.data[self.peaks[i],1:]) #total probablility not conserved
			# x += np.sum(self.data[self.peaks[i],1:])
		# print "x/N =", x/self.max_agg_size



class PlotExperimental():

	def __init__(self, file = "Gavins_data/weight_time_297_min.dat"):
		self.file = str(file)
		self.read_data()

	def read_data(self): #reads in data and uses numpy to make it array like objects

		data = np.loadtxt(self.file)

		self.area_x = []
		self.agg_size_y = []

		for i in range(len(data)):
			self.area_x.append(data[i,0])
			self.agg_size_y.append(data[i,1])
		print(data)

		self.total = 0
		for i in range(len(data)):
			self.total += self.agg_size_y[i]
		print self.total

		return (self.area_x,self.agg_size_y)

	def plot(self):
		plt.plot(self.area_x, (self.agg_size_y/self.total), label = "experimental") #should be horizonta line at 1 but is decaying
		# plt.plot(np.log(self.area_x), np.log((self.agg_size_y/self.total)), label = "experimental") #logged
		plt.ylabel("N/N_total")
		plt.title("Smoluchowski Aggregation Model")
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,200,y1,y2))
		# plt.legend()
		plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		plt.show()


def main():
	sim = Smoluchowski(200, 20000, 1, 10)
	# exp = PlotExperimental().plot()
	sim.run_sim()
	# sim.run_analytical()

	sim.plot()
	sim.find_peaks() #WHAT TIME? IS THIS PLOTTED FOR?
	sim.plot_hist()
main()
