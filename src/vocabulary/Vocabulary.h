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

#ifndef VOCABULARY_H_
#define VOCABULARY_H_

#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/export.hpp>
#include <utils/HistogramDistances.h>

/**
 * Representation for a generic feature vector to be used as a word in a dictionary.
 * The class represents a generic feature vector with weights for each dimension. It is used as a word in the Bag of Words sense. 
 * It defines the interface for computing a similiarity values between words and estimate a mean vector given multiple vectors.
 *
 * @author Gian Diego Tipaldi
 */


class HistogramFeatureWord {
    public:
	typedef std::list< std::vector<double> > ElementList;
	typedef std::list< std::vector<double> > WeightList;
	typedef std::list< std::vector<double> >::iterator ElementListIt;
	typedef std::list< std::vector<double> >::iterator WeightListIt;
	typedef std::list< std::vector<double> >::const_iterator ElementListCIt;
	typedef std::list< std::vector<double> >::const_iterator WeightListCIt;
	
	/** 
	 * Constructor. It creates the feature word providing the feature vector, a distance function and the weights' vector.
	 *
	 * @param histogram The feature vector representing this word.
	 * @param distance The distance function used to compute the similarity between feature vectors.
	 * @param weights The weights' vector to scale the individual dimension of the feature vector.
	 */
	HistogramFeatureWord(const std::vector<double>& histogram = std::vector<double>(), const HistogramDistance<double>* distance = NULL, const std::vector<double>& weights = std::vector<double>());
	
	/** Returns the histogram used to represent the feature vector. */
	inline const std::vector<double>& getHistogram() const
	    {return m_histogram;}
	
	/** Returns the vector used to represent the mean of different words. It is used during the KMeans algorithm. */    
	inline const std::vector<double>& getMean() const
	    {return m_mean;}
	
	/** Returns the histogram used to represent the weights' vector. */
	inline const std::vector<double>& getWeights() const
	    {return m_weights;}
	    
	/** Returns the similarity between feature words. Mainly used during KMeans. */
	double sim(const HistogramFeatureWord* other) const;
	
	/** Returns the similarity between the feature word and a feature vector. */
	double sim(const std::vector<double>& histogram) const;

	/** Returns the similarity between the feature word and a feature vector, considering the weights for each dimension. */
	double sim(const std::vector<double>& histogram, const std::vector<double>& weights) const;
	
	/** Merges the current feature vector with @param other. Mainly used during KMeans. */
	void merge(HistogramFeatureWord* other);
	
	/** Returns the elements belonging to this word during KMeans. */
	inline const std::list< std::vector<double> >& getElements() const
	    {return m_elements;}
	    
	/** Sets the distance function to be used for copmuting the similarity. */
	inline void setDistance(const HistogramDistance<double>* distance)
	    {m_distance = distance;}
	
    protected:
	std::vector<double> m_histogram; /**< The feature vector as histogram. */
	std::vector<double> m_mean; /**< The feature vector as the mean of the elements of this word. Mainly used during KMeans. */
	std::vector<double> m_weights; /**< The weights' vector. */
	unsigned int m_number; /**< The number of elements for this word. Mainly used for KMeans. */
	const HistogramDistance<double>* m_distance; /**< The distance function. */
	std::list< std::vector<double> > m_elements; /**< The elements of this word. Mainly used during KMeans. */
	std::list< std::vector<double> > m_elementsWeights; /** The weights' vector for the elements of this word. Mainly used during KMeans. */
	
	friend class boost::serialization::access;
	
	/** Serializes the class using boost::serialization. */ 
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);
	
};

/** Representation of a vocabulary as a vector of words. */
typedef std::vector< HistogramFeatureWord > HistogramVocabulary;

// template<typename HistogramDistance>
// class HistogramFeatureVocabulary: std::vector< HistogramFeatureWord<HistogramFeatureWord> > {
//     
// };


template<class Archive>
void HistogramFeatureWord::serialize(Archive& ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(m_histogram);
    ar & BOOST_SERIALIZATION_NVP(m_mean);
    ar & BOOST_SERIALIZATION_NVP(m_weights);
    ar & BOOST_SERIALIZATION_NVP(m_number);
    ar & BOOST_SERIALIZATION_NVP(m_elements);
}

#endif
