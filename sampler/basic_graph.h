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

#ifndef BASICGRAPH_H
#define BASICGRAPH_H

#include <vector>
using namespace std;

class AdjacencyList  {
	public:
		AdjacencyList(int _src)  {
			src = _src;
		}
		~AdjacencyList()  {}

		void add_adjacency(int _dst)  {
			neighbors.push_back(_dst);
		}

		int get_src()  {
			return src;
		}

		int get_adjacency(int _i)  {
			return neighbors[_i];
		}

		int size()  {
			return neighbors.size();
		}

	private:
		int src;
		vector<int> neighbors;
};

class BasicGraph  {
	public:
		BasicGraph(int _n)  {
			num_nodes = _n;
			edges = new AdjacencyList*[num_nodes];
			for(int i = 0; i < num_nodes; i++)
				edges[i] = new AdjacencyList(i);
		}
		~BasicGraph()  {
			for(int i = 0; i < num_nodes; i++)
				delete edges[i];
			delete [] edges;
		}

		AdjacencyList* get_adjacencies(int _i)  {
			return edges[_i];
		}

		void add_edge(int _i, int _j)  {
			edges[_i]->add_adjacency(_j);
		}

	private:
		AdjacencyList** edges;
		int num_nodes;
};

#endif
