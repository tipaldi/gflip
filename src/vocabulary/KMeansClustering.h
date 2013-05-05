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

#ifndef KMEANSLUSTERING_H_
#define KMEANSLUSTERING_H_

#include <vector>

/** 
 * Implement the K-Means clustering algorithm.
 *
 * The general implementation requires a ClusterType to implement the sim function, 
 * which returns the similarity between two clusters in the [0,1] range, and the merge function,
 * which merges two clusters into one.
 * 
 * @author Gian Diego Tipaldi
 *
 */

template <typename ClusterType>
class KMeansClustering {
	public:
		
	/** 
	 * Default constructor. It set the maximum iterations for the clustering and the minimum error difference for convergence.
	 *
	 */
	KMeansClustering(unsigned int maxIterations, double minError);
	
	/** 
	 * Cluster the @param points into clusters. It initialize the centroids with the @param seeds.
	 * The number of clusters is the size of @param seeds. The @param seeds are modified to hold the clusters.
	 * 
	 */
	void clusterPoints(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds) const;
	
	/** 
	 * Cluster the @param points into clusters. It initialize the centroids with the @param seeds.
	 * The number of clusters is the size of @param seeds. The @param seeds are modified to hold the clusters.
	 * This overloaded version returns also the point assignments.
	 * 
	 */
	void clusterPoints(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds, std::vector< std::vector<unsigned int> >& assignment) const;
	
	/** 
	 * Initialize the @param points into the centroids @param seeds. 
	 * The number of clusters is the size of @param seeds.
	 * The strategy to initialize them is defined by the functor @param strategy.
	 * 
	 */
	template<template <typename Type > class Strategy>
	inline void initializeClusters(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds) const
		{Strategy<ClusterType> initializer; initializer(points, seeds);}
	
	protected:
	
	unsigned int m_maxIterations; /**< The maximum number of iterations. */
	double m_minError; /**< The minimum error difference. */
    
};


/** 
 * Implement the Forgy initialization for the K-Means clustering algorithm.
 *
 * @author Gian Diego Tipaldi
 *
 */

template <typename ClusterType>
class ForgyKmeansInitialization {
	public:
		void operator()(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds);
};

/** 
 * Implement the Random Partition initialization for the K-Means clustering algorithm.
 *
 * @author Gian Diego Tipaldi
 *
 */

template <typename ClusterType>
class RandomKmeansInitialization {
	public:
		void operator()(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds);
};

/** 
 * Implement the KMeans++ initialization for the K-Means clustering algorithm.
 *
 * @author Gian Diego Tipaldi
 *
 */

template <typename ClusterType>
class PlusPlusKmeansInitialization {
	public:
		void operator()(std::vector<ClusterType>& points, std::vector<ClusterType>& seeds);
};

#include "KMeansClustering.hpp"

#endif
