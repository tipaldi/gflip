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

#include <vocabulary/Vocabulary.h>

HistogramFeatureWord::HistogramFeatureWord(const std::vector<double>& histogram, const HistogramDistance<double>* distance, const std::vector<double>& weights):
    m_histogram(histogram),
    m_mean(histogram),
    m_number(1),
    m_distance(distance),
    m_weights(weights)
{
    if(!m_distance){
	m_distance = &standardEuclideanDistance;
    }
    if(m_weights.size() != m_mean.size()){
	m_weights.resize(m_mean.size(), 1.);
    }
    m_elements.push_back(m_histogram);
    m_elementsWeights.push_back(m_weights);
}

double HistogramFeatureWord::sim(const HistogramFeatureWord* other) const
{
    double similarity = 0.;
    WeightListCIt firstW = m_elementsWeights.begin();
    for(ElementListCIt first = m_elements.begin(); first != m_elements.end(); first++, firstW++){
	WeightListCIt secondW = other->m_elementsWeights.begin();
	for(ElementListCIt second = other->m_elements.begin(); second != other->m_elements.end(); second++, secondW++){
	    similarity += exp(-m_distance->distance(*first, *firstW, *second, *secondW));
	}
    }
	///FIXME test the commented stuff
//     return similarity/double(m_elements.size() * other->m_elements.size());
    return other? exp(-m_distance->distance(m_mean, m_weights, other->m_mean, other->m_weights)) : 0.;
}

double HistogramFeatureWord::sim(const std::vector<double>& histogram) const{
    return exp(-m_distance->distance(m_mean, histogram));
}

double HistogramFeatureWord::sim(const std::vector<double>& histogram, const std::vector<double>& weights) const{
    return exp(-m_distance->distance(m_mean, m_weights, histogram, weights));
}

	
void HistogramFeatureWord::merge(HistogramFeatureWord* other)
{
    if(!other || m_mean.size() != other->m_mean.size()) return;
    double normalizer = 1./double(m_number + other->m_number);
    for(unsigned int i = 0; i < m_mean.size(); i++){
// 	m_mean[i] = normalizer * (double(m_number) * m_mean[i] + double(other->m_number) * other->m_mean[i]);
// 	m_weights[i] = (double(m_number) * m_weights[i] + double(other->m_number) * other->m_weights[i]);
	m_mean[i] = (m_weights[i] * m_mean[i] + other->m_weights[i] * other->m_mean[i])/(m_weights[i] + other->m_weights[i]);
	m_weights[i] = m_weights[i] + other->m_weights[i];
    }
    m_number += other->m_number;
    m_elements.splice(m_elements.end(), other->m_elements, other->m_elements.begin(), other->m_elements.end());
    m_elementsWeights.splice(m_elementsWeights.end(), other->m_elementsWeights, other->m_elementsWeights.begin(), other->m_elementsWeights.end());
		
		if(!other->m_elements.size()){
			other->m_elements.push_back(other->m_histogram);
			other->m_elementsWeights.push_back(other->m_weights);
		}
}
