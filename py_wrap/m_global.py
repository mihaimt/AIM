from matplotlib import pyplot
import numpy
import m_analyse

def plot_energy(d):
	path = d['output_path']
	efile= d['output_e']
	file_name = m_analyse.remove_spaces(path+efile)
	out_name  = m_analyse.remove_spaces(path+"total_energy_cons.pdf")

	if m_analyse.remove_spaces(d['energy']) == 'yes': 
		data = numpy.genfromtxt(file_name, delimiter = ' ', skip_footer = 0)
		tdata = numpy.transpose(data)
		time =  tdata[0]
		energy = tdata[1]
		pyplot.figure(figsize = (10,10))
		pyplot.plot(time, energy/10**6, 'r.')
		pyplot.xlabel("time [c.u.]", fontsize = 20)
		pyplot.ylabel("Energy in the system [MJ]", fontsize = 20)
		pyplot.tick_params(axis='both', which='major', labelsize=15)
		pyplot.tick_params(axis='both', which='minor', labelsize=12)
		pyplot.savefig(out_name)
	
	

