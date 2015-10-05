"""
Copyright (c) 2010, Matt Berger and Bigyan Ankur Mukherjee
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list
  of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list
  of conditions and the following disclaimer in the documentation and/or other
  materials provided with the distribution.
- Neither the name of the University of Utah nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import ConfigParser
import os
import subprocess
import sys


class RunError(Exception):
	def __init__(self, msg):
		self.value = msg
	def __str__(self):
		return 'Error: ' + self.value

def runcommand(args):
	try:
		print os.getcwd()
		print args
		subprocess.check_call(args)
	except OSError as e:
		raise RunError(args[0] + ': Fork falied with error \'' + e.strerror + '\'')
	except subprocess.CalledProcessError as e:
		raise RunError(args[0] + ': Execution failed with returncode = ' + repr(e.returncode))
		

def runUniform(config, pathdir, infile, outfile):
       if config.has_section("uniform"):
					args = []
					#args.append('bin/bash')
					#args.append('-c')
					args.append("./" + pathdir + "/" + config.get("uniform", "exec_name"))
					args.append(infile)
					args.append(outfile)

					# required

					args.append(config.get("uniform", "camera_res_x"))
					args.append(config.get("uniform", "camera_res_y"))
					args.append(config.get("uniform", "scan_res"))

					# optional

					if config.has_option("uniform", "min_range"):
						args.append("min_range")
						args.append(config.get("uniform", "min_range"))

					if config.has_option("uniform", "max_range"):
						args.append("max_range")
						args.append(config.get("uniform", "max_range"))

					if config.has_option("uniform", "num_stripes"):
						args.append("num_stripes")
						args.append(config.get("uniform", "num_stripes"))

					if config.has_option("uniform", "laser_fov"):
						args.append("laser_fov")
						args.append(config.get("uniform", "laser_fov"))

					if config.has_option("uniform", "peak_threshold"):
						args.append("peak_threshold")
						args.append(config.get("uniform", "peak_threshold"))

					if config.has_option("uniform", "std_threshold"):
						args.append("std_threshold")
						args.append(config.get("uniform", "std_threshold"))

					if config.has_option("uniform", "additive_noise"):
						args.append("additive_noise")
						args.append(config.get("uniform", "additive_noise"))

					if config.has_option("uniform", "laser_smoother"):
						args.append("laser_smoother")
						args.append(config.get("uniform", "laser_smoother"))

					if config.has_option("uniform", "registration_error"):
						args.append("registration_error")
						args.append(config.get("uniform", "registration_error"))

					if config.has_option("uniform", "normal_type"):
						args.append("normal_type")
						args.append(config.get("uniform", "normal_type"))

					if config.has_option("uniform", "pca_knn"):
						args.append("pca_knn")
						args.append(config.get("uniform", "pca_knn"))

					if config.has_option("uniform", "random_sample_rotation"):
						args.append("random_sample_rotation")
						args.append(config.get("uniform", "random_sample_rotation"))

					runcommand(args)

try:
	config = ConfigParser.ConfigParser()
	config.read(sys.argv[1])

	infile = config.get("sampler", "infile")
	outfile = config.get("sampler", "outfile")
	pathdir = config.get("sampler", "pathdir")
except Exception:
	raise RunError('Input and Output files must be specified')

	
if not os.access(infile, os.R_OK):
	raise RunError('Permission denied: ' + infile)

#os.chdir(pathdir)
runUniform(config, pathdir, infile, outfile)
