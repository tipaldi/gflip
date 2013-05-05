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

#include <feature/Descriptor.h>
#include <feature/ShapeContext.h>
#include <feature/BetaGrid.h>
#include <feature/InterestPoint.h>
#include <vocabulary/Vocabulary.h>
#include <utils/HistogramDistances.h>
#include <vocabulary/KMeansClustering.h>
#include <vocabulary/HierarchicalKMeansClustering.h>
#include <geometry/point.h>

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <utility>

#include <sys/time.h>

#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <kdtree++/kdtree.hpp>
#include <Eigen/Core>

#include <limits>


typedef HierarchicalKMeansClustering<HistogramFeatureWord> VocabularyClustering;
typedef std::vector< std::vector< InterestPoint *> > InterestPointLog;
typedef std::vector< std::vector< unsigned int> > LabelLog;

struct VMeasure {
    double homogeneity;
    double completeness;
};

struct LabelledPoint {
    typedef double value_type;
    double point[2];
    unsigned int label;
    double& operator[](unsigned int i) { return point[i];}
    double operator[](unsigned int i) const { return point[i];}
};

typedef KDTree::KDTree<2,LabelledPoint> SearchTree;
typedef std::pair<KDTree::KDTree<2,LabelledPoint>::const_iterator, KDTree::KDTree<2,LabelledPoint>::distance_type> SearchResult;

std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec){
	out << "[ ";
	for(std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); it++){
		out << *it << " ";
	}
	out << "]" << std::endl;
	return out;
}

VMeasure testCluster(const HistogramVocabulary& clusters, const LabelLog& tree, const InterestPointLog& points, unsigned int label, double beta){
    // Checking the  assignment
    Eigen::MatrixXd assignment = Eigen::MatrixXd::Constant(clusters.size(), label, std::numeric_limits< double >::min());
    for(unsigned int i = 0; i < points.size(); i++) {
	for(unsigned int j = 0; j < points[i].size(); j++) {
	    LabelledPoint point, closestPoint;
	    point[0]=points[i][j]->getPosition().x;
	    point[1]=points[i][j]->getPosition().y;
// 	    SearchResult closest = tree.find_nearest(point);
	    unsigned int cluster = 0, realCluster = tree[i][j];
	    double maxSim = 0.;
	    std::vector<double> description;
	    std::vector<double> weights;
	    points[i][j]->getDescriptor()->getWeightedFlatDescription(description, weights);
	    for(unsigned int w = 0; w < clusters.size(); w++) {
		double sim = clusters[w].sim(description, weights);
		if(sim > maxSim) {
		    maxSim = sim;
		    cluster = w;
		}
	    }
	    assignment(cluster, realCluster) = assignment(cluster, realCluster) + 1;
// 	    assignment(0, realCluster) = assignment(0, realCluster) + 1;
// 	    assignment(cluster, 0) = assignment(cluster, 0) + 1;
	}
    }
    
    // Computing the entropy based measure
//     std::cout << "assignment = " << assignment << std::endl;
    double number = assignment.sum();
//     std::cout << "number = " << number << std::endl;
    assignment *= 1./number;
    Eigen::MatrixXd assignmentC = assignment.colwise().sum();
    Eigen::MatrixXd assignmentK = assignment.rowwise().sum();
//     std::cout << "assignment = " << assignment << std::endl;
//     std::cout << "assignmentK = " << assignmentK << std::endl;
//     std::cout << "assignmentC = " << assignmentC << std::endl;
    double entropyKandC = - (assignment.array() * (assignment.array().log())).sum();
    double entropyK = - (assignmentK.array() * (assignmentK.array().log())).sum();
    double entropyC = - (assignmentC.array() * (assignmentC.array().log())).sum();
/*    std::cout << std::endl;
    std::cout << "entropyKC = " << entropyKandC << std::endl;
    std::cout << "entropyK = " << entropyK << std::endl;
    std::cout << "entropyC = " << entropyC << std::endl;
    std::cout << "entropyKC - entropyK= " << entropyKandC - entropyK << std::endl;
    std::cout << "entropyKC - entropyC= " << entropyKandC - entropyC << std::endl;
    std::cout << "entropyCgivenK = " << (entropyKandC - entropyK)/entropyC << std::endl;
    std::cout << "entropyKgivenC = " << (entropyKandC - entropyC)/entropyK << std::endl;*/
    VMeasure result;
    result.homogeneity = entropyC ? 1. - (entropyKandC - entropyK)/entropyC : 1.;
    result.completeness = entropyK ? 1. - (entropyKandC - entropyC)/entropyK : 1.;
    return result;
}


int main(int argc, char **argv){

    
    std::vector< std::string > filenames;
    std::string outfile("Vocabulary");
    unsigned int maxFeatures = 1000, distanceType = 2, samplingType = 1;
    double beta = 1.;
		
		unsigned int vocabularyLevels = 3, vocabularySize = 5;
    
    int i = 1;
    while(i < argc){
	if(strncmp("-maxFeatures", argv[i], sizeof("-maxFeatures")) == 0 ){
	    maxFeatures = atoi(argv[++i]);
	    i++;
	} else if(strncmp("-outfile", argv[i], sizeof("-outfile")) == 0 ){
	    outfile = argv[++i];
	    i++;
	} else if(strncmp("-levels", argv[i], sizeof("-levels")) == 0 ){
	    vocabularyLevels = atoi(argv[++i]);
	    i++;
	} else if(strncmp("-size", argv[i], sizeof("-size")) == 0 ){
	    vocabularySize = atoi(argv[++i]);
	    i++;
	} else if(strncmp("-distanceType", argv[i], sizeof("-distanceType")) == 0 ){
	    distanceType = atoi(argv[++i]);
	    i++;
	} else if(strncmp("-samplingType", argv[i], sizeof("-samplingType")) == 0 ){
		samplingType = atoi(argv[++i]);
		i++;
	} else if(strncmp("-beta", argv[i], sizeof("-beta")) == 0 ){
		beta = atof(argv[++i]);
		i++;
	} else {
	    filenames.push_back(argv[i++]);
	}
    }
    
    
    HistogramDistance<double> *dist = NULL;
    
    std::string distance("");
    switch(distanceType){
	case 0:
	    dist = new EuclideanDistance<double>();
	    distance = "euclid";
	    break;
	case 1:
	    dist = new Chi2Distance<double>();
	    distance = "chi2";
	    break;
	case 2:
	    dist = new SymmetricChi2Distance<double>();
	    distance = "symchi2";
	    break;
	case 3:
	    dist = new BatthacharyyaDistance<double>();
	    distance = "batt";
	    break;
	case 4:
	    dist = new KullbackLeiblerDistance<double>();
	    distance = "kld";
	    break;
	case 5:
	    dist = new JensenShannonDistance<double>();
	    distance = "jsd";
	    break;
	default:
	    std::cerr << "Wrong distance type" << std::endl;
	    exit(-1);
    }
	
	switch (samplingType){
	case 0:
		std::cerr << "low variance sampling" << std::endl;
		break;
	case 1:
		std::cerr << "random uniform sampling" << std::endl;
		break;
	default:
		std::cerr << "Wrong sampling type" << std::endl;
		exit(-1);
	}
    
    std::cerr << "Processing files:\t";
    for(std::vector<std::string>::const_iterator it = filenames.begin(); it != filenames.end(); it++) {
	std::cerr << *it << " ";
    }
    
    std::cerr << "\nDistance:\t\t" << distance 
	      << "\nSize:\t\t\t" << maxFeatures 
	      << "\nLevels:\t\t\t" << vocabularyLevels
	      << "\nSize:\t\t\t" << vocabularySize
	      << std::endl;
    
    unsigned int fileFeatures = ceil(double(maxFeatures)/double(filenames.size()));

		
    VocabularyClustering clustering(30, 0.001, vocabularySize);
    boost::mt19937 rng;
    
		HistogramVocabulary currentVocabulary;
		currentVocabulary.reserve(maxFeatures);

		std::cout << "Reading files: " << std::endl; 
		for(std::vector<std::string>::const_iterator it = filenames.begin(); it != filenames.end(); it++) {
			std::ifstream inputStream(it->c_str());
			boost::archive::binary_iarchive inputArchive(inputStream);
			std::vector< std::vector< InterestPoint *> > points;
			inputArchive >> BOOST_SERIALIZATION_NVP(points);

			if (samplingType == 0){
				///low variance sampling
				double step = double(points.size() -1) / double(fileFeatures);
				boost::uniform_real<double> generator(0, step);
				std::cerr << "initial seed in " << step << " step= " << step << std::endl;
				unsigned int seed = generator(rng);
				for (unsigned int p = 0; p < fileFeatures; p++){
					unsigned int index1 = seed + double(p)*step;
					while(points[index1].size() < 1){
						index1 -= 1;
					}
					unsigned int max = points[index1].size() - 1;
					boost::uniform_int<int> generator2(0, max);
					unsigned int index2 = generator2(rng);
					std::vector<double> description;
					std::vector<double> weights;
// 			 		std::cout << "p= " << p << " i1= "<< index1 << " i2= " << index2 << " " << max << std::endl;
					points[index1][index2]->getDescriptor()->getWeightedFlatDescription(description, weights);
	//		 		std::cout << description << std::endl;
					HistogramFeatureWord word(description, dist, weights);
	// 				HistogramFeatureWord word(description, dist);
					currentVocabulary.push_back(word);
				}
			}
	
			if (samplingType == 1){
				///uniform sampling
				boost::uniform_int<int> generator(0, points.size() - 1);
				for(unsigned int p = 0; p < fileFeatures; p++){
						unsigned int index1 = generator(rng); ///here are samples generated (this is on feature files, not on log files!)
						unsigned int max = points[index1].size() - 1;
						if(points[index1].size() < 1){
					p--;
					continue;
						}
						boost::uniform_int<int> generator2(0, max);
						unsigned int index2 = generator2(rng);
						std::vector<double> description;
						std::vector<double> weights;
	//			    std::cout << "p= " << p << " i1= "<< index1 << " i2= " << index2 << " " << max << std::endl;
						points[index1][index2]->getDescriptor()->getWeightedFlatDescription(description, weights);
	//			    std::cout << description << std::endl;
						HistogramFeatureWord word(description, dist, weights);
						currentVocabulary.push_back(word);
				}
			}
	
	
		
			for(unsigned int j = 0; j < points.size(); j++){
				for(unsigned int k = 0; k < points[j].size(); k++){
						delete points[j][k];
				}
			}
						
			inputStream.close();
			std::cout << *it << " " << points.size() << " | " << currentVocabulary.size() << std::endl;
		}
    
		KMeansClustering<HistogramFeatureWord> KMeans(30, 0.001);
		HistogramVocabulary clusters2(vocabularySize);
		KMeans.initializeClusters<PlusPlusKmeansInitialization>(currentVocabulary, clusters2);
		KMeans.clusterPoints(currentVocabulary, clusters2);
		std::ostringstream outfileK;
		outfileK << outfile << "_0_" << vocabularySize << "KMEANS.voc";
		std::ofstream outputStreamK(outfileK.str().c_str());
		boost::archive::binary_oarchive outputArchiveK(outputStreamK);
		outputArchiveK << BOOST_SERIALIZATION_NVP(clusters2);
		std::cout << "Writing Vocabulary: " << outfileK.str() << " Size = " << clusters2.size() << std::endl;
		
		
		
    double bestScore = -1;
    HistogramVocabulary bestVocabulary;
    for(unsigned int i = 0; i < vocabularyLevels; i++){
			HistogramVocabulary clusters(currentVocabulary);
			
			clustering.clusterPoints< PlusPlusKmeansInitialization >(currentVocabulary, clusters, i, vocabularySize);
		
		
			// Test the filenames
			std::ifstream inputStream(filenames[0].c_str());
			boost::archive::binary_iarchive inputArchive(inputStream);
			std::vector< std::vector< InterestPoint *> > points;
			inputArchive >> BOOST_SERIALIZATION_NVP(points);
			KDTree::KDTree<2,LabelledPoint> tree;
			unsigned int label = 0, features = 0;
			
			// Creating ground truth labels
			LabelLog labelLog(points.size());
			for(unsigned int j = 0; j < points.size(); j++) {
				labelLog[j].resize(points[j].size());
				for(unsigned int k = 0; k < points[j].size(); k++) {
					LabelledPoint point, closestPoint;
					point[0]=points[j][k]->getPosition().x;
					point[1]=points[j][k]->getPosition().y;
					SearchResult closest = tree.find_nearest(point, 0.04);
					features++;
					if(closest.first == tree.end()) {
						point.label = label++;
						tree.insert(point);
						labelLog[j][k] = point.label;
					} else {
						labelLog[j][k] = closest.first->label;
					}
				}
			}
		
			std::cout << "Testing vocabulary on " << features << " interest points" << std::endl;
		
			std::string bar(50,' ');
			bar[0] = '#';
			unsigned int progress = 0;
				
		
		
			VMeasure result;// = testCluster(clusters, labelLog, points, label, beta);
			double score = 0;//(1.+ beta) * result.homogeneity * result.completeness / (beta * result.homogeneity + result.completeness);
// 			unsigned int vocabularySize = std::numeric_limits< unsigned int >::max();
		/*	if(score > bestScore){
					std::cout << "Vocabulary with " << clusters.size() << " words, tested on " << label << " labels, scored " << score << " c = " << result.completeness << " h = " << result.homogeneity << std::endl;
					bestScore = score;
					bestVocabulary = clusters;
			}*/
// 			currentVocabulary.swap(clusters);
		
		// 	std::cout << "Vocabulary with " << clusters.size() << " words, tested on " << label << " labels, scored " << score << std::endl;
/*			unsigned int dendit = 1;
			for(HistogramVocabulary::iterator it = clusters.begin(); it != clusters.end(); it++, dendit++) {
		// 	    std::cout << "Merging cluster " << it->second.first << " with " << it->second.second << ", distance = " << it->first << std::endl;
				unsigned int currentProgress = (dendit*100)/(clusters.size());
				if (progress < currentProgress){
					progress = currentProgress;
					bar[progress/2] = '#';
					std::cout << "\rTesting vocabulary  [" << bar << "] " << progress << "%" << std::flush;
				}*/
				result = testCluster(clusters, labelLog, points, label, beta);
				score = (1.+ beta) * result.homogeneity * result.completeness / (beta * result.homogeneity + result.completeness);
				std::cout << "Vocabulary with " << clusters.size() << " words, tested on " << label << " labels, scored " << score << " c = " << result.completeness << " h = " << result.homogeneity << std::endl;
				if(score > bestScore){
					std::cout << "Improved score from " << bestScore << " to " << score << std::endl;
					bestScore = score;
					bestVocabulary = clusters;
// 					vocabularySize = clusters.size();
				}/* else if(fabs(score - bestScore) < 0.001 && clusters.size() < vocabularySize){
					std::cout << "Improved size from " << bestScore << " to " << score << std::endl;
					bestVocabulary = clusters;
// 					vocabularySize = clusters.size();
				}*/
/*			}
			std::cout << "\rTesting vocabulary  [" << bar << "] 100% done." << std::endl;*/
		
	// 	    HistogramVocabulary clusters(20);
	// 	    boost::mt19937 rng;
	// 	    boost::uniform_smallint<int> generator(0, maxFeatures);
	// 	    for(int j = 0; j < 20; j++){
	// 		clusters[j] = currentVocabulary[generator(rng)];
	// 	    }

	// 	    clustering.clusterPoints(currentVocabulary, clusters);

			std::ostringstream outfileL;
			outfileL << outfile << "_" << i << "_" << vocabularySize << ".voc";

			std::ofstream outputStreamL(outfileL.str().c_str());
		// 	    boost::archive::xml_oarchive archive(outputStream);
			boost::archive::binary_oarchive outputArchiveL(outputStreamL);
			
			outputArchiveL << BOOST_SERIALIZATION_NVP(clusters);
		
			std::cout << "Writing Vocabulary: " << outfileL.str() << std::endl;
		}
		std::ostringstream outfileB;
		outfileB << outfile << "_" << maxFeatures << "_B.voc";

		std::ofstream outputStream(outfileB.str().c_str());
	// 	    boost::archive::xml_oarchive archive(outputStream);
		boost::archive::binary_oarchive outputArchive(outputStream);
		
		outputArchive << BOOST_SERIALIZATION_NVP(bestVocabulary);
	
		std::cout << "Writing Vocabulary: " << outfileB.str() << std::endl;
}

