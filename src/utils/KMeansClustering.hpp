//
//
// GFLIP - Geometrical FLIRT Phrases for Large Scale Place Recognition
// Copyright (C) 2012-2013 Gian Diego Tipaldi and Luciano Spinello and Wolfram
// Burgard
//
// This file is part of GFLIP.
//
// GFLIP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GFLIP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GFLIP.  If not, see <http://www.gnu.org/licenses/>.
//

// #include <KMeansClustering.h>

#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>

template <typename ClusterType>
KMeansClustering<ClusterType>::KMeansClustering(unsigned int maxIterations, double minError):
    m_maxIterations(maxIterations)
    , m_minError(minError)
{
    
}

template <typename ClusterType>
void KMeansClustering<ClusterType>::clusterPoints(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds) const
{
	std::vector< std::vector<unsigned int> > assignment;
	clusterPoints(points, seeds, assignment);
}

template <typename ClusterType>
void KMeansClustering<ClusterType>::clusterPoints(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds, std::vector< std::vector<unsigned int> >& assignment) const
{
	if(points.size() < seeds.size()) {
		seeds = points;
		return;
	}
	double oldError = 0;
	for(unsigned int i = 0; i < m_maxIterations; i++) {
// 	std::cout << "Cluster centers: ";
// 	for(unsigned int c = 0; c < seeds.size(); c++) {
// 	    std::cout << seeds[c].getMean() << ", ";
// 	}
// 	std::cout << std::endl;
		double error = 0;
		assignment.clear();
		assignment.resize(seeds.size());
		for(unsigned int p = 0; p < points.size(); p++) {
			unsigned int bestCluster = 0;
			double maxSim = 0;
			for(unsigned int c = 0; c < seeds.size(); c++) {
				double sim = seeds[c].sim(&points[p]);
				if(sim > maxSim) {
					maxSim = sim;
					bestCluster = c;
				}
			}
			error += maxSim;
			assignment[bestCluster].push_back(p);
// 	    std::cout << "Point " << points[p].getMean() << " assgined to " << bestCluster << " with center " << seeds[bestCluster].getMean() << std::endl;
		}
		for(unsigned int c = 0; c < seeds.size(); c++) {
	    if(assignment[c].size()) seeds[c] = points[assignment[c].front()];
	    for(unsigned int p = 1; p < assignment[c].size(); p++) {
				seeds[c].merge(&points[assignment[c][p]]);
	    }
		}
// 		std::cout << "Iteration " << i << ", error difference = " << fabs(error - oldError) << ", min error = " << m_minError << std::endl;
// 		std::cout << "             , error = " << error << ", old error = " <<  oldError << std::endl;
		if(fabs(error - oldError) < m_minError) {break;}
		oldError = error;
	}
}


template <typename ClusterType>
void ForgyKmeansInitialization<ClusterType>::operator()(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds)
{
	if(points.size() < seeds.size()) {
		unsigned int seedSize = seeds.size();
		seeds = points;
		seeds.resize(seedSize);
		return;
	}
	std::mt19937 rng;	
	std::uniform_int_distribution<int> distribution(0, points.size() - 1);
	auto generator = std::bind ( distribution, rng );
// 	std::variate_generator<std::mt19937&, std::uniform_int<int> > generator(rng, distribution);
	for(unsigned int i = 0; i < seeds.size(); i++) {
		seeds[i] = points[generator()];
	}
}

template <typename ClusterType>
void RandomKmeansInitialization<ClusterType>::operator()(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds)
{
	if(points.size() < seeds.size()) {
		unsigned int seedSize = seeds.size();
		seeds = points;
		seeds.resize(seedSize);
		return;
	}
	std::mt19937 rng;
	std::uniform_int_distribution<int> distribution(0, seeds.size() - 1);
	auto generator = std::bind ( distribution, rng );
// 	std::uniform_int<int> distribution(0, seeds.size() - 1);
// 	std::variate_generator<std::mt19937&, std::uniform_int<int> > generator(rng, distribution);
	std::vector< std::vector<unsigned int> > assignment(seeds.size());
	for(unsigned int p = 0; p < points.size(); p++) {
		unsigned int bestCluster = generator();
		assignment[bestCluster].push_back(p);
	}
	for(unsigned int c = 0; c < seeds.size(); c++) {
		if(assignment[c].size()) seeds[c] = points[assignment[c].front()];
		for(unsigned int p = 0; p < assignment[c].size(); p++) {
			seeds[c].merge(&points[assignment[c][p]]);
		}
	}
}

template <typename ClusterType>
void PlusPlusKmeansInitialization<ClusterType>::operator()(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds)
{
	if(points.size() < seeds.size()) {
		unsigned int seedSize = seeds.size();
		seeds = points;
		seeds.resize(seedSize);
		return;
	}
	std::mt19937 rng;
	std::uniform_real_distribution<double> distribution(0, 1);
	std::uniform_int_distribution<int> distribution2(0, seeds.size() - 1);
	auto generator2 = std::bind(distribution2, rng);
// 	std::variate_generator<std::mt19937&, std::uniform_int<int> > generator2(rng, distribution2);
	auto generator = std::bind(distribution, rng);
// 	std::variate_generator<std::mt19937&, std::uniform_real<double> > generator(rng, distribution);
	double maxCumulative = 0;
	std::vector<double> cumulative(points.size());
	
	// Pick a first at random
	seeds[0] = points[generator2()];
	// Initialize the distances
	for(unsigned int j = 0; j < points.size(); j++) {
		double distance = 1. - seeds[0].sim(&points[j]);
		cumulative[j] = maxCumulative + distance * distance;
		maxCumulative = cumulative[j];
	}
	
	// Loop over the cluster seeds
	for(unsigned int i = 1; i < seeds.size(); i++) {
		// Sample the new seed
		std::vector<double>::iterator it = std::lower_bound(cumulative.begin(), cumulative.end(), generator() * maxCumulative);
		unsigned int index = std::distance(cumulative.begin(), it);
		seeds[i] = points[index];
		
		// Update the distances
		double distance = 1. - seeds[i].sim(&points[0]);
		distance = distance * distance;
		double correction = 0;
		if (distance < cumulative[0]) {
			correction = cumulative[0] - distance;
			cumulative[0] = fabs(cumulative[0] - correction);
		}		
		for(unsigned int j = 1; j < points.size(); j++) {
			cumulative[j] = fabs(cumulative[j] - correction);
			double distance = 1. - seeds[i].sim(&points[j]);
			distance = distance * distance;
			double oldDistance = cumulative[j] - cumulative[j-1];
			if (distance < oldDistance) {
				correction += oldDistance - distance;
				cumulative[j] = fabs(cumulative[j] - (oldDistance - distance));
			}
		}
		maxCumulative = cumulative.back();
	}
	
// 	std::binary_search<>();
}

