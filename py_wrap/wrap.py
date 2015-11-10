import os
import numpy
from sys import argv
import re
from m_execute import read_parameter, execute_it
import m_analyse


	
script, input = argv

d = read_parameter(input)

execute_it(d)

m_analyse.plott_all(d)


