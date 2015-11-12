
#Author: Mihai Tomozeiu#

import os
import numpy
from sys import argv
import re
from m_execute import read_parameter, execute_it, execution
from m_analyse import remove_spaces as rs
import m_analyse
import m_global


import os, errno

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: 
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def add_units(d):
	c_par      = open(rs(d['solver_input']), 'r')

        for linee in c_par.readlines():
		if linee[0:12] == "   output_dt":
                	timestep = float(linee[12:-1])
                        

                if linee[0:8] == "   x_max":
                        length = float(linee[8:-1])

                if linee[0:9] == "   n_cell":
                        n_cell = float(linee[9:-1])




        c_par.close()
        d['timestep']   = timestep
        d['cell_length']= length/n_cell
 


def write_parac(input_c_par, param_name, var_ar, main_directory, d):

	for e in var_ar:
		mkdir_p(main_directory+param_name+"/"+e)
		c_par      = open(name_c_par, 'r')
		print name_c_par
		path_new   = main_directory+param_name+ "/" + e + "/" 
		c_new      = open(path_new + "param.in", 'w')
		for linee in c_par.readlines():
			if linee == "   output_path  .\n":
				c_new.write("output_path "+ path_new[0:-1]+"\n")
#				print "   output_path "+ path_new[0:-1] 	
			elif (param_name in linee) and (len(linee)<50):
				c_new.write("   "+param_name+" "+e+ "\n")
#				print "   "+param_name+" "+e 

			else:
	       			c_new.write(linee)

			if linee[0:12] == "   output_dt":
				timestep = float(linee[12:-1])
				print timestep

			if linee[0:8] == "   x_max":
				length = float(linee[8:-1])
                        if linee[0:9] == "   n_cell":
                                n_cell = float(linee[9:-1])




		c_par.close()
		c_new.close()
		c = {}
		for key in d:
			c[key] = d[key]
		#print(c['
    	        c['solver_call'] = d['solver_call']
                c['solver_input']= path_new+"param.in"
                c['solver_dir']  = d['solver_dir']
		c['output_path'] = path_new
		c['timestep']   = timestep
		c['cell_length']= length/n_cell
		execute_it(c)	
		
	        m_analyse.plott_all(c)
	        m_global.plot_energy(c)


#def execute_them(list_input):
	
		



	
script, input = argv

d = read_parameter(input)


if rs(d['var_param']) == 'no':
	if rs(d['execute']) == 'yes':
		execute_it(d)

	add_units(d)
	
	m_analyse.plott_all(d)
	m_global.plot_energy(d)

else:
	varfile = rs(d['parameter_set_path'])+rs(d['parameter_set'])
	vf = open(varfile, 'r')
	set_var = []
	for line in vf.readlines():
		set_var = set_var + [line[:-1]]

	vf.close()
	for element in set_var:
		varpar = element

		param_file = open(rs(d['parameter_set_path'])+varpar, 'r')
		i = 0
		var_ar = []
		for line in param_file.readlines():
			if i == 0:
				param_name = rs(line[0:-1])
				print "="*50
				print "Exploring parameter: ", param_name
				print "="*50
			else:
				var_ar = var_ar + [rs(line[0:-1])]
			i = i + 1
		param_file.close()
		print var_ar


#Reading the c++ parameter file

		name_c_par = rs(d['solver_input'])
		output_main = rs(d['output_path'])
		write_parac(name_c_par, param_name, var_ar, output_main, d)
	

