import os
import numpy
from sys import argv
import re
from m_execute import read_parameter, execute_it
import m_analyse
import m_global

	
script, input = argv

d = read_parameter(input)
if m_analyse.remove_spaces(d['execute']) == 'yes':
	execute_it(d)



m_analyse.plott_all(d)
m_global.plot_energy(d)

