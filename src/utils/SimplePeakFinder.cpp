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

#include "SimplePeakFinder.h"

SimplePeakFinder::SimplePeakFinder(double minValue, double minDifference):
    PeakFinder(),
    m_minValue(minValue),
    m_minDifference(minDifference)
{

}

void SimplePeakFinder::findPeaks(const std::vector<double>& signal, std::vector<unsigned int>& indexes) const	{
    indexes.clear();
    for (unsigned int i = 1; i < signal.size() - 1; i++){
	if(isPeak(signal, i)){
	    indexes.push_back(i);
	}
    }
}

bool SimplePeakFinder::isPeak(const std::vector<double>& signal, unsigned int index) const {
    return signal[index] > m_minValue && 
	   signal[index] - signal[index - 1] > m_minDifference &&
	   signal[index] - signal[index + 1] > m_minDifference;
}

void SimplePeakFinder::findPeaks(const std::vector< std::vector<double> >& signal, std::vector< std::vector<unsigned int> >& indexes) const{
    indexes.resize(signal.size());
    for (unsigned int i = 0; i < signal.size(); i++){
	findPeaks(signal[i], indexes[i]);
    }
}

bool SimplePeakFinder::isPeak(const std::vector< std::vector<double> >& signal, unsigned int index1, unsigned int index2) const{
    return isPeak(signal[index1],index2);
}
