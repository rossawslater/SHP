import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys
class PlotExperimental():

	def __init__(self, file, saveTrue): #= "Gavins_data/weight_time_297_min.dat"):
		self.file = str(file)
		self.read_data()
		self.plotname = self.file.split("/")[-1]
		self.plotname = self.plotname.split(".")[0]
		self.saveTrue = saveTrue

	def read_data(self): #reads in data and uses numpy to make it array like objects

		data = np.loadtxt(self.file)

		self.area_x = []
		self.agg_size_y = []

		for i in range(len(data)):
			self.area_x.append(data[i,0])
			self.agg_size_y.append(data[i,1])
		# print(data)

		self.total = 0
		for i in range(len(data)):
			self.total += self.agg_size_y[i]
		# print self.total

		return (self.area_x,self.agg_size_y)

	def plot(self):
		plt.plot(self.area_x, (self.agg_size_y)) #should be horizonta line at 1 but is decaying
		# plt.plot(np.log(self.area_x), np.log((self.agg_size_y/self.total)), label = "experimental") #logged
		plt.ylabel("Number of aggregates of size N")
		plt.xlabel("Aggregate area in microns squared")
		plt.title(self.plotname)
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,200,y1,50000))
		# plt.legend()
		# plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		if self.saveTrue:
			plt.savefig(str(os.getcwd() + "/Gavins_data/Plots/Unweighted/" + self.plotname + ".png"))
		# plt.show()
			plt.clf()
		else:
			plt.show()

	def plot_weighted(self):
		plt.plot(self.area_x, self.agg_size_y) #should be horizonta line at 1 but is decaying
		# plt.plot(np.log(self.area_x), np.log((self.agg_size_y/self.total)), label = "experimental") #logged
		plt.ylabel("Number of aggregates of size N, multiplied by N")
		plt.xlabel("Aggregate area in microns squared")
		plt.title(self.plotname)
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,200,y1,120000))
		# plt.legend()
		# plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		if self.saveTrue:
			plt.savefig(str(os.getcwd() + "/Gavins_data/Plots/Weighted/" + self.plotname + ".png"))
			plt.clf()
		else:
			plt.show()
	def plot_weighted_divide_total(self):
		plt.plot(self.area_x, (self.agg_size_y/self.total)) #should be horizonta line at 1 but is decaying
		# plt.plot(np.log(self.area_x), np.log((self.agg_size_y/self.total)), label = "experimental") #logged
		plt.ylabel("N/N_total")
		plt.xlabel("Aggregate area in microns squared")
		plt.title(self.plotname)
		x1,x2,y1,y2 = plt.axis()
		plt.axis((x1,200,y1,1))
		# plt.legend()
		# plt.legend(bbox_to_anchor=(0.5, 0.8), loc=2, borderaxespad=0.)
		if self.saveTrue:
			plt.savefig(str(os.getcwd() + "/Gavins_data/Plots/Weighted_Divide_Total/" + self.plotname + ".png"))
		# plt.show()
			plt.clf()
		else:
			plt.show()


def main():

	if len(sys.argv) == 2 and sys.argv[1] == "Save":
		saveTrue = True
	else:
		saveTrue = False

	for filename in glob.glob('Gavins_data/*.dat'):
		exp = PlotExperimental(str(filename), saveTrue)
		if "weight" in filename:
			exp.plot_weighted()
			exp.plot_weighted_divide_total()
		else:
			exp.plot()
main()
