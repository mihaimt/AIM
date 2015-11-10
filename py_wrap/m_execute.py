import os
import numpy
from sys import argv
import re

def execution(d):
	cal    = d['solver_call']
	input  = d['solver_input']
	path   = d['solver_dir']

	cals = cal.split(' ')
	cl = []
	for o in cals:
		if len(o)>1:
			cl = cl + [o]	
	string = path+"/"+cl[0]+" "+cl[1]+ " " +input
	return string 

def clean(string):
	sa, sb, sc = '', '', ''
        a = string.split(' ')
        for s in a:
                sa = sa + s
	if "\t" in sa:
		sc = sa.strip("\t")
		for s in sc:
			sa = sa + s
        return sa

def extract_info(s):
	if len(s)==2:
		key = s[0]
		value = s[1]
		return clean(key), str(value)
	else:
		return "error", "error" 


def read_parameter(input):
	in_file = open(input, 'r')
	d = {}
	for line in in_file:
       		s =  line[:-1].split("=")
        	key, value = extract_info(s)

        	if key <> "error":
                	d[str(key)] = str(value)
	in_file.close()

	return d

def execute_it(d):
	string = execution(d)
	print "="*20
	print "Executing: ", string
	print "="*20
	os.system(string)
	

