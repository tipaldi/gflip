/* *
 * GFLIP - Geometrical FLIRT Phrases for Large Scale Place Recognition
 * Copyright (C) 2012-2013 Gian Diego Tipaldi and Luciano Spinello and Wolfram
 * Burgard
 *
 * This file is part of GFLIP.
 *
 * GFLIP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GFLIP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with GFLIP.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HIERARCHICALKMEANSLUSTERING_H_
#define HIERARCHICALKMEANSLUSTERING_H_

#include <utils/KMeansClustering.h>

#include <vector>
#include <cmath>

/** 
 * Implement the Hierarchical K-Means clustering algorithm.
 *
 * @author Gian Diego Tipaldi
 *
 */

template <typename ClusterType>
class HierarchicalKMeansClustering {
	public:
	/** 
	 * Default constructor. It set the maximum iterations for the clustering, the minimum cluster difference and the fanout.
	 *
	 */
	HierarchicalKMeansClustering(unsigned int maxIterations, double minError, unsigned int fanout = 10);
	
	/** 
	 * Cluster the @param points into clusters. It initialize the centroids with the provided Initialization class.
	 * The number of clusters is the size of @param seeds. The @param seeds are modified to hold the clusters.
	 * 
	 */
	template<template <typename Type > class Strategy>
	void clusterPoints(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds) const
	{ 
		unsigned int levels = double(seeds.size())/double(m_fanout); 
		clusterPoints<Strategy>(points, seeds, levels, m_fanout);
	}
	
	/** 
	 * Cluster the @param points into clusters. It initialize the centroids with the @param seeds.
	 * The number of clusters is the size of @param seeds. The @param seeds are modified to hold the clusters.
	 * 
	 */
	template<template <typename Type > class Strategy>
	void clusterPoints(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds, unsigned int levels, unsigned int fanout = 0) const
	{
		fanout = fanout + bool(!fanout) * m_fanout;
		std::vector<ClusterType> localSeeds(fanout);
		std::vector< std::vector<unsigned int> > assignment;
		Strategy<ClusterType> initializer;
		initializer(points, localSeeds);
		m_clustering.clusterPoints(points, localSeeds, assignment);
		if(levels == 0 || localSeeds.size() < fanout) {
			seeds.swap(localSeeds);
			return;
		}
		seeds.clear();
		seeds.reserve(pow(fanout, levels + 1));
		for(unsigned int i = 0; i < fanout; i++) {
			std::vector<ClusterType> localPoints(assignment[i].size());
			for(unsigned int j = 0; j < assignment[i].size(); j++){
				localPoints[j] = points[assignment[i][j]];
			}
			clusterPoints< Strategy >(localPoints, localSeeds, levels - 1, fanout);
			seeds.insert(seeds.end(), localSeeds.begin(), localSeeds.end());
		}
	}
	
	protected:
	unsigned int m_maxIterations; /**< The maximum number of iterations. */
	double m_minError; /**< The maximum number of iterations. */
	unsigned int m_fanout; /**< The number of clusters per level. */
	KMeansClustering< ClusterType > m_clustering; /**< The low level KMeans class for clustering each node. */
    
};


#include "HierarchicalKMeansClustering.hpp"

#endif
