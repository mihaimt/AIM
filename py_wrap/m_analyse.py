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
	pos, den, vel, pre = [], [], [], []
	steps = 0
	for e in data_bunch:
		steps = steps + 1
		pos = pos + [e[:,0]]
		den = den + [e[:,1]]
		vel = vel + [e[:,2]]
		pre = pre + [e[:,3]]
	
	size = len(e[:,0])
	sizes = [steps, size]

	return path, sizes, pos, den, vel, pre


def plott_one(arr, sizes, name, path):

	prename = path+"/"+name
	steps, size = sizes

	fs = [20,20]
	pyplot.figure(figsize = (fs[0],fs[1]))
	pyplot.matshow(arr, extent=[0,size,steps,0], aspect='auto')
	pyplot.hot()
	pyplot.xlabel("position [c.u.]")
	pyplot.ylabel("time [c.u.]")


	cbar = pyplot.colorbar()

	cbar.set_label(name+' [c.u.]', rotation=90)

	pyplot.savefig(prename+"_mat.pdf")


	pyplot.cla()
	pyplot.clf()

	fs = [20,20]
	pyplot.figure(figsize = (fs[0],fs[1]))
	pyplot.imshow(arr, extent=[0,size,steps,0], aspect='auto')
	pyplot.hot()

	pyplot.xlabel("position [c.u.]")
	pyplot.ylabel("time [c.u.]")

	cbar = pyplot.colorbar()
	cbar.set_label(name+' [c.u.]', rotation=90)

	pyplot.savefig(prename+"_ims.pdf")



def remove_spaces(string):
	ns = ''
	for i in string:
		if i not in [' ']:
			ns = ns + i
	return ns


def plott_all(d):
	path, sizes,  pos, den, vel, pre = read(d)
	nd = {}
	nd['pressure'] = pre
	nd['density']  = den
	nd['velocity'] = vel

	for key in d:
		if key in ['pressure', 'density', 'velocity']:
			if remove_spaces(d[key]) == 'yes':
				print key, 'done'
				plott_one(nd[key], sizes, key, path) 

	print "=Fin="




#print numpy.transpose(pos)
	
#Move this one back up
#execute_it(d)

#plot section


