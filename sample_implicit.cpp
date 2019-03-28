#include "shape_loader.h"

#include <iostream>
#include <vector>
#include <array>
#include <fstream>

int main(int argc, char const *argv[])
{
	if(argc != 4)
	{
		std::cerr<<"call with <implicit file> <cube size> <out name>"<<std::endl;
		return EXIT_FAILURE;
	}

	MPUSoup implicit(argv[1]);

	Vector3 ll_bound;
	Vector3 ur_bound;
	implicit.bounds(ll_bound, ur_bound);

	double cube_size = atof(argv[2]);

	std::cout<<"dx: "<<(ur_bound.entry(0) - ll_bound.entry(0))<<" dy: "<<(ur_bound.entry(1) - ll_bound.entry(1))<<" dz: "<<(ur_bound.entry(2) - ll_bound.entry(2))<<std::endl;

	const int nx = (ur_bound.entry(0) - ll_bound.entry(0))/cube_size;
	const int ny = (ur_bound.entry(1) - ll_bound.entry(1))/cube_size;
	const int nz = (ur_bound.entry(2) - ll_bound.entry(2))/cube_size;

	std::vector<std::array<double, 3> > points;
	std::array<double, 3> tmp;

	for(int i = -1; i <= nx; ++i)
	{
		for(int j = -1; j <= ny; ++j)
		{
			for(int k = -1; k <= nz; ++k)
			{
				const Vector3 pt1(ll_bound.entry(0) + i*cube_size,     ll_bound.entry(1) + j*cube_size,   ll_bound.entry(2) + k*cube_size);
				const double f1 = implicit.function(pt1);

				{
					const Vector3 pt2(ll_bound.entry(0) + (i+1)*cube_size, ll_bound.entry(1) + (j)*cube_size, ll_bound.entry(2) + (k)*cube_size);
					const double f2 = implicit.function(pt2);

					if( f1*f2 < 0)
					{
						tmp[0] = ll_bound.entry(0) + (i+0.5)*cube_size;
						tmp[1] = ll_bound.entry(1) + (j)*cube_size;
						tmp[2] = ll_bound.entry(2) + (k)*cube_size;
						points.push_back(tmp);
					}
				}

				{
					const Vector3 pt2(ll_bound.entry(0) + (i)*cube_size, ll_bound.entry(1) + (j+1)*cube_size, ll_bound.entry(2) + (k)*cube_size);
					const double f2 = implicit.function(pt2);

					if( f1*f2 < 0)
					{
						tmp[0] = ll_bound.entry(0) + (i)*cube_size;
						tmp[1] = ll_bound.entry(1) + (j+0.5)*cube_size;
						tmp[2] = ll_bound.entry(2) + (k)*cube_size;
						points.push_back(tmp);
					}
				}

				{
					const Vector3 pt2(ll_bound.entry(0) + (i)*cube_size, ll_bound.entry(1) + (j)*cube_size, ll_bound.entry(2) + (k+1)*cube_size);
					const double f2 = implicit.function(pt2);

					if( f1*f2 < 0)
					{
						tmp[0] = ll_bound.entry(0) + (i)*cube_size;
						tmp[1] = ll_bound.entry(1) + (j)*cube_size;
						tmp[2] = ll_bound.entry(2) + (k+0.5)*cube_size;
						points.push_back(tmp);
					}
				}
			}
		}


		std::cout<<i<<"/"<<nx<<std::endl;

	}

	std::ofstream file(argv[3]);
	for(size_t i = 0; i < points.size(); ++i)
	{
		file << points[i][0] << " " << points[i][1] << " " << points[i][2] << "\n";
	}
	file.close();

	return EXIT_SUCCESS;
}