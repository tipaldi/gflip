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

#include <iostream>
#include <string>
#include <string.h>
#include <sstream>
#include <utility>

#include <sys/time.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/multi_array.hpp>

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

struct timeval detectTime, describeTime;

void detectLog(){
    unsigned int i = 0;
    unsigned int position = m_sensorReference.tell();
    m_sensorReference.seek(0,END);
    unsigned int last = m_sensorReference.tell();
    m_sensorReference.seek(0);
    
    std::string bar(50, ' ');
    
    struct timeval start, end;
    gettimeofday(&start, NULL);
    while(!m_sensorReference.end()){
	bar[(m_sensorReference.tell()*50)/last] = '#';
	std::cout << "\rDetecting points  [" << bar << "] " << (m_sensorReference.tell()*100)/last << "%" << std::flush;
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
    std::cout << "done. Time elapsed " << double(detectTime.tv_sec) + 1e-06 * double(detectTime.tv_usec) << std::endl;
}

void countLog(){
    double flirtNum = 0.;
    uint count = 0;
    
    std::string bar(50, ' ');

    for(unsigned int i = 0; i < m_pointsReference.size(); i++){
	bar[(i*50)/m_pointsReference.size()] = '#';
	std::cout << "\rCounting points  [" << bar << "] " << (i*100)/m_pointsReference.size() << "%" << std::flush;
	if(m_pointsReference[i].size()){
	    flirtNum += m_pointsReference[i].size();
	    count++;
	}
    }
    flirtNum=count?flirtNum/double(count):0.;
    std::cout << "done: Total: " << flirtNum*double(count) << ", Avg: " << flirtNum << std::endl;
}

void describeLog(){
    unsigned int i = 0;
    unsigned int position = m_sensorReference.tell();
    m_sensorReference.seek(0,END);
    unsigned int last = m_sensorReference.tell();
    m_sensorReference.seek(0);

    std::string bar(50, ' ');

    struct timeval start, end;
    gettimeofday(&start, NULL);
    while(!m_sensorReference.end()){
	bar[(m_sensorReference.tell()*50)/last] = '#';
	std::cout << "\rDescribing points  [" << bar << "] " << (m_sensorReference.tell()*100)/last << "%" << std::flush;
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
    std::cout << "done. Time elapsed " << double(describeTime.tv_sec) + 1e-06 * double(describeTime.tv_usec) << std::endl;
}

void writeGspan(std::ostream& out, const std::vector <std::vector <InterestPoint* > >& points){
    for(unsigned int i = 0; i < points.size(); i++) {
	out << "t # " << i << std::endl;
	for(unsigned int v = 0; v < points[i].size(); v++) {
	    out << "v " << v << " a " << points[i][v]->getPosition().x << " " << points[i][v]->getPosition().y << " 0.0" << std::endl;
	}
	if(points[i].size() < 2) continue;
	std::vector<bool> usedVertexes(points[i].size(), false);
	usedVertexes[0] = true;
	unsigned int vertexCount = 1;
	while(vertexCount != points[i].size()){
	    unsigned int minIndex1 = 0;
	    unsigned int minIndex2 = 0;
	    double minDistance = 1e10;
	    for(unsigned int j = 0; j < points[i].size(); j++) {
		if(!usedVertexes[j]) continue;
		for(unsigned int k = 0; k < points[i].size(); k++) {
		    if(usedVertexes[k]) continue;
		    Point2D delta = points[i][j]->getPosition() - points[i][k]->getPosition();
		    double currentDistance = delta * delta;
		    if(currentDistance < minDistance){
			minDistance = currentDistance;
			minIndex1 = j;
			minIndex2 = k;
		    }
		}
	    }
	    out << "e " << minIndex1 << " " << minIndex2 << " sparse" << std::endl;
	    usedVertexes[minIndex2] = true;
	    vertexCount++;
	}
    }
}

int main(int argc, char **argv){

    std::string filename("");
    unsigned int scale = 5, dmst = 2, window = 3, detectorType = 0, descriptorType = 0, distanceType = 2;
    double baseSigma = 0.2, sigmaStep = 1.4, minPeak = 0.34, minPeakDistance = 0.001;
    bool useMaxRange = false, gspan = false;
    
    int i = 1;
    while(i < argc){
	if(strncmp("-filename", argv[i], sizeof("-filename")) == 0 ){
	    filename = argv[++i];
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
	} else if(strncmp("-gspan", argv[i], sizeof("-gspan")) == 0 ){
	    gspan = true;
	    i++;
	} else {
	    i++;
	}
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
	    m_detector = m_detectorNormalEdge;
	    detector = "edge";
	    break;
	case 2:
	    m_detectorNormalBlob = new NormalBlobDetector(m_peakMinMax, scale, baseSigma, sigmaStep, window);
	    m_detector = m_detectorNormalBlob;
	    detector = "blob";
	    break;
	case 3:
	    m_detectorRange = new RangeDetector(m_peakMinMax, scale, baseSigma, sigmaStep);
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
    
    std::cerr << "Processing file:\t" << filename << "\nDetector:\t\t" << detector << "\nDescriptor:\t\t" << descriptor << "\nDistance:\t\t" << distance << std::endl;
    
    m_sensorReference.seek(0,END);
    unsigned int end = m_sensorReference.tell();
    m_sensorReference.seek(0,BEGIN);

    m_pointsReference.resize(end + 1);
    m_posesReference.resize(end + 1);
    
    detectLog();
    
    countLog();
//     return 1;
    
    describeLog();
    
    std::string outfile = filename.substr(0, filename.find_last_of("."));
    outfile.append(".flt");
    
    if(gspan){
	std::string gspanfile = filename.substr(0, filename.find_last_of("."));
	gspanfile.append(".gspan");
	std::ofstream gspanStream(gspanfile.c_str());
	writeGspan(gspanStream, m_pointsReference);
    }

    std::ofstream outputStream(outfile.c_str());
//     boost::archive::xml_oarchive archive(outputStream);
    boost::archive::binary_oarchive archive(outputStream);    
    
    
    std::cout << "Saving the points reference " << BOOST_PP_STRINGIZE(m_pointsReference) << std::endl;
    archive << BOOST_SERIALIZATION_NVP(m_pointsReference);
    std::cout << "Saving the points generators " << BOOST_PP_STRINGIZE(m_betaGenerator->getPhiEdges()) << std::endl;
    if(m_betaGenerator) {
			
			archive << boost::serialization::make_nvp("phiEdges", m_betaGenerator->getPhiEdges());
			archive << boost::serialization::make_nvp("rhoEdges", m_betaGenerator->getRhoEdges());
    } else if(m_shapeGenerator) {
			archive << boost::serialization::make_nvp("phiEdges", m_shapeGenerator->getPhiEdges());
			archive << boost::serialization::make_nvp("rhoEdges", m_shapeGenerator->getRhoEdges());
    }
    std::cout << "done" << std::endl;
    
    outputStream.close();

    // Uncomment for Serialization test
//     const InterestPoint* testPoint = m_pointsReference.back().back();
//     unsigned int numPoints = 0;
//     for(unsigned int j = 0; j < m_pointsReference.size(); numPoints += m_pointsReference[j++].size()){
// 	std::cerr << numPoints << " " << m_pointsReference[j].size() << std::endl;
// 	
//     }
//     
//     std::cout << "Filename: " << outfile 
// 	      << ", Scans: " << m_pointsReference.size() 
// 	      << ", Points: " << numPoints 
// 	      << ", Avg: " << double(numPoints)/double(m_pointsReference.size()) << std::endl << std::endl;
//     std::cout << "Test Point: end, end" << std::endl
// 	      << testPoint->getPosition() << std::endl
// 	      << testPoint->getScale() << std::endl
// 	      << testPoint->getScaleLevel() << std::endl
// 	      << typeid(*testPoint->getDescriptor()).name() << std::endl
// 	      << testPoint->getSupport().size() << std::endl;


}

