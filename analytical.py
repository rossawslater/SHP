def run_analytical(init_no_particles, sim_length ):
    tau = 2/init_no_particles
    for t in range(0,sim_length):


        for j in range(1,2):#,self.max_agg_size):

            self.aggregates[j] = self.init_no_particles * (1 + t/self.tau)**(-j-1)
            self.get_j_counts_at_step(j,t)

        # self.ttau.append(t/self.tau)

        sys.stdout.write(" Simulation progress: %.1f%%   \r" %(t*100/float(sim_length)))
        sys.stdout.flush()
        # print(t/float(self.sim_length))
    for i in self.data:
        print(i[1]/init_no_particles)
        # print(i[1])
        # print("ping")
