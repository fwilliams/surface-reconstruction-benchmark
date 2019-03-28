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

#ifndef UNIONFIND_H
#define UNIONFIND_H

template<typename Entity>
class UFEntity  {
	public:
		UFEntity(Entity _entity)  {
			entity = _entity;
			parent = this;
			rank = 0;
		}
		~UFEntity()  {}

		Entity entity;
		UFEntity* parent;
		int rank;
};

template<typename Entity>
class UFSet  {
	public:
		UFSet()  {}
		~UFSet()  {}

		bool perform_union(UFEntity<Entity>* entity1, UFEntity<Entity>* entity2)  {
			UFEntity<Entity>* entity1_root = this->find(entity1);
			UFEntity<Entity>* entity2_root = this->find(entity2);
			//??????????????????
			const int r1 = entity1_root->rank;
			const int r2 = entity2_root->rank;

			if(r1 > r2)  {
				entity2_root->parent = entity1_root;
				return true;
			}
			else if(r1 < r2)  {
				entity1_root->parent = entity2_root;
				return true;
			}
			else if(entity1_root != entity2_root)  {
				entity2_root->parent = entity1_root;
				entity1_root->rank = entity1_root->rank+1;
				return true;
			}

			return false;
		}

		UFEntity<Entity>* find(UFEntity<Entity>* entity)  {
			if(entity->parent == entity)
				return entity;
			// else
			entity->parent = this->find(entity->parent);
			return entity->parent;
		}
};

#endif
