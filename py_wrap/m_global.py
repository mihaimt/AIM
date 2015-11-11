from matplotlib import pyplot
import numpy

file_name = ""
def read_energy(file_name):
	data = numpy.genfromtxt(file_name, delimiter = ' ', skip_footer = 0)
	print data[0]
	

read_energy(file_name)
