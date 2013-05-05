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

#include <feature/Detector.h>
#include <feature/ShapeContext.h>
#include <feature/BetaGrid.h>
#include <feature/RangeDetector.h>
#include <feature/CurvatureDetector.h>
#include <feature/NormalBlobDetector.h>
#include <feature/NormalEdgeDetector.h>
#include <sensorstream/CarmenLog.h>
#include <sensorstream/LogSensorStream.h>
#include <sensorstream/SensorStream.h>
#include <utils/SimpleMinMaxPeakFinder.h>
#include <utils/HistogramDistances.h>
#include <vocabulary/Vocabulary.h>

#include <gflip/gflip_engine.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include <iostream>
#include <string>
#include <string.h>
#include <sstream>
#include <utility>
#include <map>

#include <sys/time.h>

LogSensorStream m_sensorReference(NULL,NULL);

CurvatureDetector *m_detectorCurvature = NULL;
NormalBlobDetector *m_detectorNormalBlob = NULL;
NormalEdgeDetector *m_detectorNormalEdge = NULL;
RangeDetector *m_detectorRange = NULL;
Detector* m_detector = NULL;

BetaGridGenerator *m_betaGenerator = NULL;
ShapeContextGenerator *m_shapeGenerator = NULL;
DescriptorGenerator *m_descriptor = NULL;

std::vector< std::vector<InterestPoint *> > m_pointsReference;
std::vector< OrientedPoint2D > m_posesReference;

struct timeval detectTime, describeTime, vocabularyTime;

HistogramVocabulary histogramVocabulary;

gflip_engine *m_gfpMatcher = NULL;
int m_type = 2, m_neighborood = 50;
std::vector< std::vector< int > >m_bowReference;

void help(){
    std::cerr << "Geometrical FLIRT Phrases for Large Scale Place Recognition in 2D Range Data" << std::endl
    << "G. D. Tipaldi, L. Spinello, W. Burgard -- Int. Conf. Robotics and Automation (ICRA) 2013" << std::endl
    << std::endl
    << "Usage: generateNN -filename <logfile> [options] " << std::endl
    << "Options:" << std::endl
    << " -filename          \t The logfile in CARMEN format to process (mandatory)." << std::endl
    << " -vocabulary        \t The vocabulary file to use (default=Vocabolary.voc)." << std::endl
    << " -scale             \t The number of scales to consider (default=5)." << std::endl
    << " -dmst              \t The number of spanning tree for the curvature detector (deafult=2)." << std::endl
    << " -window            \t The size of the local window for estimating the normal signal (default=3)." << std::endl
    << " -detector          \t The type of detector to use. Available options are:" << std::endl
    << "                    \t     0 - Curvature detector [DEFAULT];" << std::endl
    << "                    \t     1 - Normal edge detector;" << std::endl
    << "                    \t     2 - Normal blob detector;" << std::endl
    << "                    \t     3 - Range detector." << std::endl
    << " -descriptor        \t The type of descriptor to use. Available options are:" << std::endl
    << "                    \t     0 - Beta - grid;" << std::endl
    << "                    \t     1 - Shape context [DEFAULT]." << std::endl
    << " -distance          \t The distance function to compare the descriptors. Available options are:" << std::endl
    << "                    \t     0 - Euclidean distance;" << std::endl
    << "                    \t     1 - Chi square distance;" << std::endl
    << "                    \t     2 - Symmetric Chi square distance [DEFAULT];" << std::endl
    << "                    \t     3 - Bhattacharyya distance;" << std::endl
    << "                    \t     4 - Kullback Leibler divergence;" << std::endl
    << "                    \t     5 - Jensen Shannon divergence." << std::endl
    << " -baseSigma         \t The initial standard deviation for the smoothing operator (default=0.2)." << std::endl
    << " -sigmaStep         \t The incremental step for the scales of the smoothing operator. (default=1.4)." << std::endl
    << " -minPeak           \t The minimum value for a peak to be detected (default=0.34)." << std::endl
    << " -minPeakDistance   \t The minimum difference with the neighbors for a peak to be detected (default=0.001)." << std::endl
    
    << " -type              \t The type of matching to perform. Available options are:" << std::endl
    << "                    \t     0 - Standard BoW;" << std::endl
    << "                    \t     1 - Geometrical FLIRT Phrases (GFP) [DEFAULT]." << std::endl
    << " -kernel            \t The kernel size to use for GFP. It must be > 2 (default = 2)." << std::endl
    << " -kbest             \t The number of nearest neighbors to be returned for each query (default = 50)." << std::endl
    << " -bow_subtype       \t The type of TF-IDF scaling to use. Available options are:" << std::endl
    << "                    \t     0 - Standard TF-IDF scaling [DEFAULT];" << std::endl
    << "                    \t     1 - Sublinear TF-IDF scaling;" << std::endl
    << "                    \t     2 - Lenght smoothing TF-IDF scaling." << std::endl
    << " -alpha_vss         \t The value of the aplha smoother to use for the lenght smoothing TF-IDF scaling (default = 0.4)." << std::endl
    << std::endl
    << "The program generates the list of nearest neighbour, according to the strategy and kernel used, for each scan present in the logfile." << std::endl
    << "The program writes the NN list in a file named <logfile>.nn " << std::endl
    << "Each line of the file corresponds to the list of NN for each laser reading of the logfile." << std::endl
    << "The data is written in a space seprated value, following the format:" << std::endl
    << "    <#NN> <#SSET> <ScanID> <ScanID> ... <ScanID>" << std::endl
    << "where #NN are teh number of nearest neigbors, #SSET the size of the scoreset and ScanID the IDs of the matched scans." << std::endl
    << std::endl;
}


bool pairComparator ( const std::pair<double, unsigned int>& l, const std::pair<double, unsigned int>& r)
   { return l.first > r.first; }


void writePoses(){
    unsigned int i = 0;
    unsigned int position = m_sensorReference.tell();
    m_sensorReference.seek(0,END);
    unsigned int last = m_sensorReference.tell();
    m_sensorReference.seek(0);
    
    std::string bar(50, ' ');
    bar[0] = '#';
    unsigned int progress = 0;
    
    while(!m_sensorReference.end()){
	unsigned int currentProgress = (m_sensorReference.tell()*100)/last;
	if (progress < currentProgress){
	    progress = currentProgress;
	    bar[progress/2] = '#';
	    std::cout << "\rWriting poses     [" << bar << "] " << (m_sensorReference.tell()*100)/last << "%" << std::flush;
	}
	const LaserReading* lreadReference = dynamic_cast<const LaserReading*>(m_sensorReference.next());
	if (lreadReference){
	    m_posesReference[i] = lreadReference->getLaserPose();
	    i++;
	}
    }
    m_sensorReference.seek(position);
    std::cout << " done." << std::endl;
}

void detectLog(){
    unsigned int i = 0;
    unsigned int position = m_sensorReference.tell();
    m_sensorReference.seek(0,END);
    unsigned int last = m_sensorReference.tell();
    m_sensorReference.seek(0);
    
    std::string bar(50, ' ');
    bar[0] = '#';
    unsigned int progress = 0;
    
    struct timeval start, end;
    gettimeofday(&start, NULL);
    while(!m_sensorReference.end()){
	unsigned int currentProgress = (m_sensorReference.tell()*100)/last;
	if (progress < currentProgress){
	    progress = currentProgress;
	    bar[progress/2] = '#';
	    std::cout << "\rDetecting points  [" << bar << "] " << (m_sensorReference.tell()*100)/last << "%" << std::flush;
	}
	const LaserReading* lreadReference = dynamic_cast<const LaserReading*>(m_sensorReference.next());
	if (lreadReference){
	    m_detector->detect(*lreadReference, m_pointsReference[i]);
	    m_posesReference[i] = lreadReference->getLaserPose();
	    i++;
	}
    }
    gettimeofday(&end,NULL);
    timersub(&end,&start,&detectTime);
    m_sensorReference.seek(position);
    std::cout << " done." << std::endl;
}

void countLog(){
    double flirtNum = 0.;
    uint count = 0;
    
    std::string bar(50, ' ');
    bar[0] = '#';

    unsigned int progress = 0;
    for(unsigned int i = 0; i < m_pointsReference.size(); i++){
	unsigned int currentProgress = (i*100)/(m_pointsReference.size() - 1);
	if (progress < currentProgress){
	    progress = currentProgress;
	    bar[progress/2] = '#';
	    std::cout << "\rCounting points   [" << bar << "] " << progress << "%" << std::flush;
	}
	if(m_pointsReference[i].size()){
	    flirtNum += m_pointsReference[i].size();
	    count++;
	}
    }
    flirtNum=count?flirtNum/double(count):0.;
    std::cout << " done.\nFound " << flirtNum << " FLIRT features per scan." << std::endl;
}

void describeLog(){
    unsigned int i = 0;
    unsigned int position = m_sensorReference.tell();
    m_sensorReference.seek(0,END);
    unsigned int last = m_sensorReference.tell();
    m_sensorReference.seek(0);

    std::string bar(50, ' ');
    bar[0] = '#';
    struct timeval start, end;
    gettimeofday(&start, NULL);
    unsigned int progress = 0;
    
    while(!m_sensorReference.end()){
	unsigned int currentProgress = (m_sensorReference.tell()*100)/last;
	if (progress < currentProgress){
	    progress = currentProgress;
	    bar[progress/2] = '#';
	    std::cout << "\rDescribing points [" << bar << "] " << progress << "%" << std::flush;
	}
	const LaserReading* lreadReference = dynamic_cast<const LaserReading*>(m_sensorReference.next());
	if (lreadReference){
	    for(unsigned int j = 0; j < m_pointsReference[i].size(); j++){
		m_pointsReference[i][j]->setDescriptor(m_descriptor->describe(*m_pointsReference[i][j], *lreadReference));
	    }
	    i++;
	}
    }
    gettimeofday(&end,NULL);
    timersub(&end,&start,&describeTime);
    m_sensorReference.seek(position);
    std::cout << " done." << std::endl;
}

struct WordResult {
  unsigned int word;
  OrientedPoint2D pose;
};

void generateBoWDescription(const OrientedPoint2D& pose, const std::vector<InterestPoint *>& pointsVector, std::multimap<double,WordResult>& signature) {
    for(unsigned int j = 0; j < pointsVector.size(); j++){
	InterestPoint * point = pointsVector[j];
	OrientedPoint2D localpose = pose.ominus(point->getPosition());
	double angle = atan2(localpose.y, localpose.x);
	unsigned int bestWord = 0;
	double bestMatch = 0.;
	std::vector<double> descriptor;
	std::vector<double> weights;
	point->getDescriptor()->getWeightedFlatDescription(descriptor, weights);
	HistogramFeatureWord word(descriptor, NULL, weights);
	for(unsigned int w = 0; w < histogramVocabulary.size(); w++) {
	    double score = histogramVocabulary[w].sim(&word);
	    if(score > bestMatch) {
		bestMatch = score;
		bestWord = w;
	    }
	}
	WordResult best; best.pose = localpose; best.word = bestWord;
	signature.insert(std::make_pair(angle,best));
    }
}

void insertBoWDescriptor(const std::multimap<double,WordResult>& signature, std::vector< int>& wordscan) {
    wordscan.resize(signature.size());
    std::vector< double > xpos(signature.size());
    std::vector< double > ypos(signature.size());
    unsigned int w = 0;
    for(std::multimap<double,WordResult>::const_iterator it = signature.begin(); it != signature.end(); ++it, ++w){
	wordscan[w] = it->second.word;
	xpos[w] = it->second.pose.x;
	ypos[w] = it->second.pose.y;
    }
    m_gfpMatcher->insert_wordscan(wordscan, xpos, ypos);
}

void buildBoWDatabase(){
    std::string bar(50, ' ');
    bar[0] = '#';

    timerclear(&vocabularyTime);
    unsigned int progress = 0;
    struct timeval start, end, diff, sum;
    for(unsigned int i = 0; i < m_pointsReference.size(); i++){
	unsigned int currentProgress = (i*100)/(m_pointsReference.size() - 1);
	if (progress < currentProgress){
	    progress = currentProgress;
	    bar[progress/2] = '#';
	    std::cout << "\rBuilding database [" << bar << "] " << progress << "%" << std::flush;
	}
	gettimeofday(&start, NULL);
	std::multimap<double,WordResult> signature;
	generateBoWDescription(m_posesReference[i], m_pointsReference[i], signature);
	insertBoWDescriptor(signature, m_bowReference[i]);
	gettimeofday(&end,NULL);
	timersub(&end,&start,&diff);
	timeradd(&vocabularyTime, &diff, &sum);
	vocabularyTime = sum;
    }
    gettimeofday(&start, NULL);
    m_gfpMatcher->prepare();
    gettimeofday(&end,NULL);
    timersub(&end,&start,&diff);
    timeradd(&vocabularyTime, &diff, &sum);
    vocabularyTime = sum;
    std::cout << " done." << std::endl;
}

void writeNN(std::ofstream& out) {
    struct timeval start, end, queryTime;
    std::string bar(50, ' ');
    bar[0] = '#';
    unsigned int progress = 0;
    for(unsigned int i =0; i < m_pointsReference.size(); i++){
	unsigned int currentProgress = (i*100)/(m_pointsReference.size() - 1);
	if (progress < currentProgress){
	    progress = currentProgress;
	    bar[progress/2] = '#';
	    std::cout << "\rComputing NN      [" << bar << "] " << progress << "%" << std::flush;
	}
	std::vector < std::pair <double, int> > *scorequery = NULL;
	gettimeofday(&start, NULL);
	m_gfpMatcher->query(m_type, m_bowReference[i],  &scorequery);
	gettimeofday(&end, NULL);  
	timersub(&end,&start,&queryTime);
	out << m_neighborood << " " << scorequery->size() << " " << double(queryTime.tv_sec) + 1e-06 * double(queryTime.tv_usec);
	int maxNN = std::min((unsigned int)m_neighborood, (unsigned int)scorequery->size());
	for(int ii=0;ii<maxNN;ii++)
	    out << " " << (*scorequery)[ii].second;
	out << std::endl;
    }
    std::cout << " done." << std::endl;
}


int main(int argc, char **argv) {

    std::string filename(""), vocabulary("Vocabulary.voc");
    unsigned int scale = 5, dmst = 2, window = 3, detectorType = 0, descriptorType = 0, distanceType = 2;
    double baseSigma = 0.2, sigmaStep = 1.4, minPeak = 0.34, minPeakDistance = 0.001;
    bool useMaxRange = false;
    
    // GFP parameters
    double alpha_vss = 0.4;
    int kernel = 2, bow_subtype = 0;
    char bag = 0;

    
    int i = 1;
    while(i < argc){
		if(strncmp("-filename", argv[i], sizeof("-filename")) == 0 ){
			filename = argv[++i];
			i++;
		} else if(strncmp("-vocabulary", argv[i], sizeof("-vocabulary")) == 0 ){
			vocabulary = argv[++i];
			i++;
		} else if(strncmp("-detector", argv[i], sizeof("-detector")) == 0 ){
			detectorType = atoi(argv[++i]);
			i++;
		} else if(strncmp("-descriptor", argv[i], sizeof("-descriptor")) == 0 ){
			descriptorType = atoi(argv[++i]);
			i++;
		} else if(strncmp("-distance", argv[i], sizeof("-distance")) == 0 ){
			distanceType = atoi(argv[++i]);
			i++;
		} else if(strncmp("-baseSigma", argv[i], sizeof("-baseSigma")) == 0 ){
			baseSigma = strtod(argv[++i], NULL);
			i++;
		} else if(strncmp("-sigmaStep", argv[i], sizeof("-sigmaStep")) == 0 ){
			sigmaStep = strtod(argv[++i], NULL);
			i++;
		} else if(strncmp("-minPeak", argv[i], sizeof("-minPeak")) == 0 ){
			minPeak = strtod(argv[++i], NULL);
			i++;
		} else if(strncmp("-minPeakDistance", argv[i], sizeof("-minPeakDistance")) == 0 ){
			minPeakDistance = strtod(argv[++i], NULL);
			i++;
		} else if(strncmp("-scale", argv[i], sizeof("-scale")) == 0 ){
			scale = atoi(argv[++i]);
			i++;
		} else if(strncmp("-dmst", argv[i], sizeof("-dmst")) == 0 ){
			dmst = atoi(argv[++i]);
			i++;
		} else if(strncmp("-window", argv[i], sizeof("-window")) == 0 ){
			window = atoi(argv[++i]);
			i++;
		} else if(strncmp("-type", argv[i], sizeof("-type")) == 0 ){
			m_type = atoi(argv[++i]);
			i++;
		} else if(strncmp("-kernel", argv[i], sizeof("-kernel")) == 0 ){
			kernel = atoi(argv[++i]);
			i++;
		} else if(strncmp("-kbest", argv[i], sizeof("-kbest")) == 0 ){
			m_neighborood = atoi(argv[++i]);
			i++;
		} else if(strncmp("-bow_subtype", argv[i], sizeof("-bow_subtype")) == 0 ){
			bow_subtype = atoi(argv[++i]);
			i++;
		} else if(strncmp("-bag", argv[i], sizeof("-bag")) == 0 ){
			bag = atoi(argv[++i]);
			i++;
		} else if(strncmp("-alpha_vss", argv[i], sizeof("-alpha_vss")) == 0 ){
			alpha_vss = strtod(argv[++i], NULL);
			i++;
		} else if(strncmp("-help", argv[i], sizeof("-localSkip")) == 0 ){
			help();
			exit(0);
		} else {
			i++;
		}
    }
    
    if(!filename.size()){
		help();
		exit(-1);
    }
    

    
    CarmenLogWriter writer;
    CarmenLogReader reader;
    
    m_sensorReference = LogSensorStream(&reader, &writer);
    
    m_sensorReference.load(filename);
    
    SimpleMinMaxPeakFinder *m_peakMinMax = new SimpleMinMaxPeakFinder(minPeak, minPeakDistance);
    
    
    std::string detector("");
    switch(detectorType){
	case 0:
	    m_detectorCurvature = new CurvatureDetector(m_peakMinMax, scale, baseSigma, sigmaStep, dmst);
	    m_detectorCurvature->setUseMaxRange(useMaxRange);
	    m_detector = m_detectorCurvature;
	    detector = "curvature";
	    break;
	case 1:
	    m_detectorNormalEdge = new NormalEdgeDetector(m_peakMinMax, scale, baseSigma, sigmaStep, window);
	    m_detectorNormalEdge->setUseMaxRange(useMaxRange);
	    m_detector = m_detectorNormalEdge;
	    detector = "edge";
	    break;
	case 2:
	    m_detectorNormalBlob = new NormalBlobDetector(m_peakMinMax, scale, baseSigma, sigmaStep, window);
	    m_detectorNormalBlob->setUseMaxRange(useMaxRange);
	    m_detector = m_detectorNormalBlob;
	    detector = "blob";
	    break;
	case 3:
	    m_detectorRange = new RangeDetector(m_peakMinMax, scale, baseSigma, sigmaStep);
	    m_detectorRange->setUseMaxRange(useMaxRange);
	    m_detector = m_detectorRange;
	    detector = "range";
	    break;
	default:
	    std::cerr << "Wrong detector type" << std::endl;
	    exit(-1);
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
    
    std::string descriptor("");
    switch(descriptorType){
	case 0:
	    m_betaGenerator = new BetaGridGenerator(0.02, 0.5, 4, 12);
	    m_betaGenerator->setDistanceFunction(dist);
	    m_descriptor = m_betaGenerator;
	    descriptor = "beta";
	    break;
	case 1:
	    m_shapeGenerator = new ShapeContextGenerator(0.02, 0.5, 4, 12);
	    m_shapeGenerator->setDistanceFunction(dist);
	    m_descriptor = m_shapeGenerator;
	    descriptor = "shape";
	    break;
	default:
	    std::cerr << "Wrong descriptor type" << std::endl;
	    exit(-1);
    }
    
    std::cerr << "Processing file:\t" << filename << 
		 "\nDetector:\t\t" << detector << 
		 "\nDescriptor:\t\t" << descriptor << 
		 "\nDistance:\t\t" << distance << 
		 "\nVocabulary:\t\t" << vocabulary << 
		 "\nNeighborood:\t\t" << m_neighborood << std::endl;
    
    std::ifstream vocabularyStream(vocabulary.c_str());
    boost::archive::binary_iarchive vocabularyArchive(vocabularyStream);
    vocabularyArchive >> histogramVocabulary;
    
    m_gfpMatcher = new gflip_engine(kernel, m_neighborood, bag, bow_subtype, alpha_vss);
    
    m_sensorReference.seek(0,END);
    unsigned int end = m_sensorReference.tell();
    m_sensorReference.seek(0,BEGIN);

    m_pointsReference.resize(end + 1);
    m_posesReference.resize(end + 1);
    m_bowReference.resize(end + 1);
    
    writePoses();
    try {
      std::string featureFile = filename.substr(0,filename.find_last_of('.')) + ".flt";
      std::ifstream featureStream(featureFile.c_str());
      boost::archive::binary_iarchive featureArchive(featureStream);
      std::cout << "Loading feature file " << featureFile << " ...";
      featureArchive >> m_pointsReference;
      std::cout << " done." << std::endl;
    } catch(boost::archive::archive_exception& exc) {
      detectLog();
      countLog();
      describeLog();
    }
    
    std::string nnFile = filename.substr(0,filename.find_last_of('.')) + ".nn";
    std::ofstream nnStream(nnFile.c_str());
    
    buildBoWDatabase();
    
    writeNN(nnStream);
    
}

