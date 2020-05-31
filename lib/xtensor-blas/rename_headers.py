import sys
import re
import os
from shutil import copyfile

args = sys.argv[1:]

expr = r"#include *[\<\"](?P<path>.+)[\>\"]"

def mkdir(path):
	dname = os.path.dirname(path)
	if not os.path.exists(dname):
		print("Making dir: ", dname)
		os.makedirs(dname)

for arg in args:
	print("Opening: ", arg)
	temp = ""
	if os.path.isdir(arg):
		continue

	arg2 = arg.replace('flens', '../xflens/')
	if arg.endswith(('.cxx', '.tcc', '.h', '.cc')):
		with open(arg) as f:
				for l in f.readlines():
					match = re.match(expr, l)
					if match:
						path = match.group('path')
						if path.startswith('cxxblas/') or path.startswith('cxxlapack/'):
							path = "xflens/" + path
							temp += '#include \"' + path + '"\n'  
						else:
							temp += l
					else:
						temp += l
		mkdir(arg2)
		with open(arg2, 'w+') as f:
			f.write(temp)
	else:
		mkdir(arg2)
		copyfile(arg, arg2)
