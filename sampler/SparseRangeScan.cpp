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

#include "SparseRangeScan.h"

SparseRangeScan::SparseRangeScan(string _file)  {
	char line_raw[256];
	ifstream file(_file.c_str());
	if (file.fail())
		return;

	for (;;)  {
		file.getline(line_raw, 256);
		if (file.eof())
			break;

		string line(line_raw);
		vector<string> tokens;
		tokenize_line(&tokens, line, " \t");
		if(tokens.size() != 5)
			continue;

		PixelCoord new_pixel = PixelCoord(atoi(tokens[0].c_str()), atoi(tokens[1].c_str()));
		Vector3 new_position = Vector3(atof(tokens[2].c_str()), atof(tokens[3].c_str()), atof(tokens[4].c_str()));
		pixels.push_back(new_pixel);
		pts.push_back(new_position);
	}
	file.close();
}

SparseRangeScan::~SparseRangeScan()  {
}

void SparseRangeScan::tokenize_line(vector<string>* tokens, const string& input, string sep)  {
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
