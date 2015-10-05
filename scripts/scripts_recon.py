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
#import OSError

class RunAlgoError(Exception):
	def __init__(self, msg):
		self.value = msg
	def __str__(self):
		return 'Error: ' + self.value

def runcommand(args):
	try:
		subprocess.check_call(args)
	except OSError as e:
		print args[0] + " : Fork failed with error : " + e.strerror
	except subprocess.CalledProcessError as e:
		print args[0] + " : Execution failed with returncode = : "  + repr(e.returncode)
		

def runPoisson(config, pathdir, infile, outdir):
	if config.has_section("poisson"):
		if not os.access(outdir + "/poisson", os.F_OK):
			os.mkdir(outdir + "/poisson")
		args = []

		args.append("./" + pathdir + "/" + config.get("poisson", "poisson_script"))
		#args.append("--in")
		args.append(indir + "/" + infile)
		#args.append("--out")
		args.append(outdir + "/poisson/" + infile.rstrip(".npts") + ".ply")

		if config.has_option("poisson", "depth"):
			args.append("--depth")
			args.append(config.get("poisson", "depth"))
		if config.has_option("poisson", "scale"):
			args.append("--scale")
			args.append(config.get("poisson", "scale"))
		if config.has_option("poisson", "solverDivide"):
			args.append("--solverDivide")
			args.append(config.get("poisson", "solverDivide"))
		if config.has_option("poisson", "isoDivide"):
			args.append("--isoDivide")
			args.append(config.get("poisson", "isoDivide"))
		if config.has_option("poisson", "samplesPerNode"):
			args.append("--samplesPerNode")
			args.append(config.get("poisson", "samplesPerNode"))
		runcommand(args)

# --- entry point --- #
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])

if config.has_option("dir_structs", "infile"):
	print 'WRONG SCRIPT: perhaps you want single_recon.py...'
	sys.exit(1)

pathdir = config.get("dir_structs", "pathdir")
indir = config.get("dir_structs", "indir")
outdir = config.get("dir_structs", "outdir")

file_num = 0
#os.chdir(pathdir)
for infile in os.listdir(indir):
	if not infile.endswith(".npts"):
		continue
	print 'run all on file ' + infile
	runPoisson(config, pathdir, infile, outdir)
	file_num=file_num+1;
