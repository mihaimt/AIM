
##Author: Mihai Tomozeiu

import os
import numpy
from sys import argv
import re
from m_execute import read_parameter, execute_it
from matplotlib import pyplot


def read(d):

	path = ''
	for ch in d['output_path']:
		if ch not in [' ']:
			path = path +ch
	
	lista = sorted(os.listdir(path))
	sol_lista = []

	for e in lista:
		if e[0:3] == 'sol':
			sol_lista = sol_lista + [e]


	data_bunch= []

	for name in sol_lista:
		data = numpy.genfromtxt(path+name, delimiter = ' ', skip_footer = 0)
		data_bunch = data_bunch + [data]




	data_transpose= numpy.transpose(data_bunch)
	pos, den, vel, pre, tem = [], [], [], [], []
	steps = 0
	for e in data_bunch:
		steps = steps + 1
		pos = pos + [e[:,0]]
		den = den + [e[:,1]]
		vel = vel + [e[:,2]]
		pre = pre + [e[:,3]]
		tem = tem + [e[:,4]]
	
	size = len(e[:,0])
	sizes = [steps, size]

	return path, sizes, pos, den, vel, pre, tem


def plott_one(arr, sizes, name, path, timestep, cell_length):


	di = {}
	di['time'] = 'ms'
	di['temperature'] = 'K'
	di['pressure'] = 'Pa'
	di['velocity'] = r'$\rm m \cdot s^{-1}$'
	di['density'] = r'$\rm kg \cdot m^{-3}$'

	prename = path+"/"+name
	steps, size = sizes

	fs = [10,10]
	pyplot.figure(figsize = (fs[0],fs[1]))
	pyplot.imshow(arr, extent=[0,size*cell_length,0,steps*timestep*1000], aspect='auto', interpolation = 'none')
	pyplot.hot()
	pyplot.xlabel("position [m]", fontsize = 20)
	pyplot.ylabel("time [ms]", fontsize = 20)
	pyplot.tick_params(axis='both', which='major', labelsize=15)
        pyplot.tick_params(axis='both', which='minor', labelsize=12)


	cbar = pyplot.colorbar()

	cbar.set_label(name+' ['+di[name]+']', rotation=90, fontsize = 20)

	pyplot.savefig(prename+"_no_interp.pdf")


	pyplot.cla()
	pyplot.clf()

	fs = [10,10]
	pyplot.figure(figsize = (fs[0],fs[1]))
	pyplot.imshow(arr, extent=[0,size*cell_length,0,steps*timestep*1000], aspect='auto')
	pyplot.hot()

	pyplot.xlabel("position [m]", fontsize = 20)
	pyplot.ylabel("time [ms]", fontsize = 20)

	cbar = pyplot.colorbar()
	cbar.set_label(name+' ['+di[name]+']', rotation=90, fontsize = 20)

	pyplot.savefig(prename+"_interp.pdf")

          

def remove_spaces(string):
	ns = ''
	for i in string:
		if i not in [' ']:
			ns = ns + i
	return ns


def reverse(arr):
	narr = []
	for e in arr[::-1]:
		narr = narr + [e]
	return narr


def plott_all(d):
	path, sizes,  pos, den, vel, pre, tem = read(d)
	nd = {}
	nd['pressure'] = reverse(pre)
	nd['density']  = reverse(den)
	nd['velocity'] = reverse(vel)
	nd['temperature'] = reverse(tem)

	numpy.savez(path+"/data.npz", nd = nd, sizes = sizes)

	for key in d:
		if key in ['pressure', 'density', 'velocity', 'temperature']:
			if remove_spaces(d[key]) == 'yes':
				print key, 'done'
				plott_one(nd[key], sizes, key, path, d['timestep'], d['cell_length'])

	print "=Fin="




#print numpy.transpose(pos)
	
#Move this one back up
#execute_it(d)

#plot section


