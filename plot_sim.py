import numpy as np
import matplotlib.pyplot as plt
import sys

class Plotter():
    def __init__(self, input):#file, min, max, step):
        self.data = np.loadtxt(input[0], skiprows = 1)
        params = np.loadtxt(input[0])[0]
        self.init_no_particles = int(params[0])
        self.sim_length = int(int(params[1])/float(params[2])) #number of time steps for simulation,
        self.dt = float(params[2])
        self.max_agg_size = (int(params[3]) + 1)#+1 so we can index from 1 not 0
        self.k = float(params[4])# self.k = 3.167E-3 #Literature value of D for Pa
        self.mu = float(params[5])

        # self.times = [103,137,162,201,252,297] #times in gavin's data
        self.times = range(int(input[1]), int(input[2]), int(input[3]))

    def plot_hists(self):
        for i in self.times:
            t = i/(self.dt)
            plt.plot(self.data[t,:])#, label = "N%i"%(i))
            x1,x2,y1,y2 = plt.axis()
            plt.axis((1,len(self.data[0,:]),y1,2000))
            plt.title("Distribution of Aggregates at %i Minutes"%(i))
            plt.ylabel("Number of Bacteria in Agg of size N")
            plt.xlabel("Aggregate Size N")
            # plt.savefig("Growth_Distribution_at_%i_mins.png"%(i))
            plt.show()

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

def main():
    plotter = Plotter(sys.argv[1:])
    plotter.plot_hists()
main()
