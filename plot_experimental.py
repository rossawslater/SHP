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
