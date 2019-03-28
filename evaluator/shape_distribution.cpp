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

#include "shape_distribution.h"

ShapeDistribution::ShapeDistribution(ShortestDistanceMap* _i2m, ShortestDistanceMap* _m2i)  {
	i2m_map = _i2m;
	m2i_map = _m2i;
}

ShapeDistribution::~ShapeDistribution()  {}

void ShapeDistribution::quartiles(vector<double>* _distribution, double& min, double& lq, double& median, double& uq, double& max)  {
	int size = _distribution->size();
	median = 0.5*(_distribution->at((size-1)/2)+_distribution->at(size/2));

	double lq_perc = 0.25*(size-1);
	int lq_low = lq_perc;
	int lq_high = lq_low+1;
	lq = (1.0-(lq_perc-lq_low))*_distribution->at(lq_low) + (lq_perc-lq_low)*_distribution->at(lq_high);

	double uq_perc = 0.75*(size-1);
	int uq_low = uq_perc;
	int uq_high = uq_low+1;
	uq = (1.0-(uq_perc-uq_low))*_distribution->at(uq_low) + (uq_perc-uq_low)*_distribution->at(uq_high);

	min = _distribution->at(0);
	max = _distribution->at(size-1);
}

bool ShapeDistribution::write_to_file(string _distFilename)  {
	string dist_filename = _distFilename;

	vector<double> distances;
	vector<double> angles;

	double mean_m2i_dist = 0, mean_m2i_angle = 0, mean_dist = 0, mean_angle = 0;
	int num_m2i = 0, num_m2i_nan = 0;
	for(int i = 0; i < m2i_map->num_correspondences(); i++)  {
		PointPair point_pair = m2i_map->getPointCorrespondence(i);
		NormalPair normal_pair = m2i_map->getNormalCorrespondence(i);
		if(point_pair.has_nan() || normal_pair.has_nan())  {
			num_m2i_nan++;
			continue;
		}
		num_m2i++;

		double next_dist = point_pair.distance();
		distances.push_back(next_dist);

		double next_angle = normal_pair.angle();
		angles.push_back(next_angle);

		if(isnan(next_angle))
			cout << "bad angle: " << normal_pair.normal1 << " : " << normal_pair.normal2 << endl;

		mean_m2i_dist += next_dist;
		mean_dist += next_dist;
		mean_m2i_angle += next_angle;
		mean_angle += next_angle;
	}

	double mean_i2m_dist = 0, mean_i2m_angle = 0;
	int num_i2m = 0, num_i2m_nan = 0;
	for(int i = 0; i < i2m_map->num_correspondences(); i++)  {
		PointPair point_pair = i2m_map->getPointCorrespondence(i);
		NormalPair normal_pair = i2m_map->getNormalCorrespondence(i);
		if(point_pair.has_nan() || normal_pair.has_nan())  {
			num_i2m_nan++;
			continue;
		}
		num_i2m++;

		double next_dist = point_pair.distance();
		distances.push_back(next_dist);

		double next_angle = normal_pair.angle();
		angles.push_back(next_angle);

		if(isnan(next_angle))
			cout << "bad angle: " << normal_pair.normal1 << " : " << normal_pair.normal2 << endl;

		mean_i2m_dist += next_dist;
		mean_dist += next_dist;
		mean_i2m_angle += next_angle;
		mean_angle += next_angle;
	}

	std::sort(distances.begin(),distances.end());
	std::sort(angles.begin(),angles.end());

	mean_m2i_dist /= ((double)num_m2i);
	mean_m2i_angle /= ((double)num_m2i);
	mean_i2m_dist /= ((double)num_i2m);
	mean_i2m_angle /= ((double)num_i2m);
	mean_dist /= (double(num_m2i+num_i2m));
	mean_angle /= (double(num_m2i+num_i2m));

	cout << "m2i nan: " << num_m2i_nan << " i2m nan: " << num_i2m_nan << endl;
	if((num_m2i+num_i2m) != distances.size() || (num_m2i+num_i2m) != angles.size())
		cout << "ERROR: sizes don't match up!" << endl;

	double min_dist, lq_dist, median_dist, uq_dist, max_dist;
	quartiles(&distances, min_dist, lq_dist, median_dist, uq_dist, max_dist);
	double min_angle, lq_angle, median_angle, uq_angle, max_angle;
	quartiles(&angles, min_angle, lq_angle, median_angle, uq_angle, max_angle);

	FILE* dist_file = fopen(dist_filename.c_str(), "wb");
	double dist_distribution[] = { min_dist, lq_dist, median_dist, uq_dist, max_dist };
	fwrite(dist_distribution, sizeof(double), 5, dist_file);
	fwrite(&mean_m2i_dist, sizeof(double), 1, dist_file);
	fwrite(&mean_i2m_dist, sizeof(double), 1, dist_file);
	fwrite(&mean_dist, sizeof(double), 1, dist_file);

	double angle_distribution[] = { min_angle, lq_angle, median_angle, uq_angle, max_angle };
	fwrite(angle_distribution, sizeof(double), 5, dist_file);
	fwrite(&mean_m2i_angle, sizeof(double), 1, dist_file);
	fwrite(&mean_i2m_angle, sizeof(double), 1, dist_file);
	fwrite(&mean_angle, sizeof(double), 1, dist_file);

	fclose(dist_file);

	return true;
}
