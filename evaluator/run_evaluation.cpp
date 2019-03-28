/*
Copyright (c) 2010, Matt Berger
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
*/

#include "recon_evaluation.h"
#include "../sampler/NPTSReader.h"
#include "../modeling/shape_loader.h"

#include <fstream>
#include <sstream>

void tokenize_line(vector<string>* tokens, const string& input, string sep)  {
	string comm;
	for (unsigned int i = 0; i < input.size(); i++)  {
		const char ch = input[i];
		bool added = false;
		for (unsigned int s = 0; s < sep.size(); s++)  {
			if (ch == sep[s])  {
				tokens->push_back(comm);
				comm = "";
				added = true;
				break;
			}
		}
		if (!added)
			comm += ch;
	}
	if (comm != "")
		tokens->push_back(comm);
}

int main(int argc, char **argv) {
	if(argc != 6)  {
		std::cerr << "usage: " << argv[0] << " mesh shape pc output_base write_correspondences" << std::endl;
		return 1;
	}
	int arg_num = 1;
	bool determine_inversion = false;

	// --- input --- //
	BasicTriMesh mesh;
   OpenMesh::IO::Options opt;
	if(!OpenMesh::IO::read_mesh(mesh, argv[arg_num++], opt))  {
		std::cerr << "unable to read mesh from file " << argv[arg_num] << std::endl;
		return 1;
	}

	string shape_filename = argv[arg_num++];
	ImplicitFunction* shape = ShapeLoader::load_shape(shape_filename);

	OrientedPointCloud pc;
	string pc_filename = argv[arg_num++];
	cout << "pc filename: " << pc_filename << endl;
	NPTSReader pc_reader(pc_filename);
	pc_reader.read_pc(&pc);

	// tease out computation time from mesh filename
	string mesh_filename = argv[1];
	string mesh_err = mesh_filename.substr(0, mesh_filename.size()-4) + ".err";

	char line_raw[256];
	double computation_time = 0;
	ifstream err_file(mesh_err.c_str());
	if (err_file.fail())
		cerr << "unable to parse out computation time: set to 0s by default." << endl;
	else  {
		for (;;)  {
			err_file.getline(line_raw, 256);
			if (err_file.eof())
				break;

			string line(line_raw);
			vector<string> tokens;
			tokenize_line(&tokens, line, " \t");

			if(tokens.size() == 3 && tokens[0].compare("execution") == 0 && tokens[1].compare("time") == 0)  {
				string str_time = tokens[2].substr(0, tokens[2].length()-1);
				computation_time += atof(str_time.c_str());
			}
		}
		err_file.close();
	}

	// --- output --- //
	string output_base = argv[arg_num++];
	int write_correspondences = atoi(argv[arg_num++]);
	bool to_write_correspondences = write_correspondences > 0;

	ReconEvaluation* evaluator = new ReconEvaluation(shape, &pc, &mesh);
	if(determine_inversion)
		evaluator->enableInversion();
	evaluator->evaluate_all();
	evaluator->write_evaluation(output_base, computation_time, to_write_correspondences);
	delete evaluator;
}
