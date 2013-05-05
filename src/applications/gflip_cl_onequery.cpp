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

#include <gflip/gflip_engine.hpp>
#include <iostream>


typedef struct
{
	double meters,angle,alpha_vss;
	int type,kernel, kbest,bow_subtype;
	char bag;
	std::string filein ;
	std::string outdir;
}sw_param_str;

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

void program_info(void)
{
	std::cout << "-i BOW input file " << std::endl;
 	std::cout << "-t kind of matching {1: standard bow, 2: GFP geometrical flirt phrases [DEFAULT]} " << std::endl;
 	std::cout << "-st TFIDF flavour {0: TFIDF, 1: Sublinear TFIDF scaling, 2: lenght smoothing TFIDF} (used only with -t 1)" << std::endl;
 	std::cout << "-alpha [0,1] smoother for length (used only with -t 1 -st 2) [0.4 DEFAULT]" << std::endl;
	std::cout << "-k [2..N] GFP kernel size [2 DEFAULT] (used only with -t 2) " << std::endl;
	std::cout << "-b bag of distance words (histograms of pairwise distances) [NO DEFAULT]" << std::endl;
	std::cout << "-kbest [0..N] returns best k results for each query [50 DEFAULT]" << std::endl;
}

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

int parse_command_line(int argc, char **argv, sw_param_str *sw_param)
{
	int i;
	
	sw_param -> type = 2;
	sw_param -> kernel = 2;	
 	sw_param -> bag = 0;
	sw_param -> outdir = "./";
	sw_param -> kbest = 50;
	sw_param ->  bow_subtype = 0;
	sw_param ->  alpha_vss = 0.4;

	for(i=0; i<argc; i++)
	{

		if(!strcmp(argv[i], "-i"))
				sw_param -> filein = argv[i+1];

		if(!strcmp(argv[i], "-b"))		
				sw_param -> bag = 1;
	 
		if(!strcmp(argv[i], "-t"))		
				sw_param -> type = atoi(argv[i+1]);

		if(!strcmp(argv[i], "-st"))		
				sw_param -> bow_subtype = atoi(argv[i+1]);

		if(!strcmp(argv[i], "-alpha"))		
				sw_param -> alpha_vss = atof(argv[i+1]);	
	
		if(!strcmp(argv[i], "-k"))		
				sw_param -> kernel = atoi(argv[i+1]);

		if(!strcmp(argv[i], "--help"))		
		{
			program_info();
			exit(1);
		}
				

		if(!strcmp(argv[i], "-kbest"))		
				sw_param -> kbest = atoi(argv[i+1]);

 	}
			
			
	if(sw_param -> kbest == 1)
	{
		std::cout << "Kbest must be > 1" << std::endl;
		exit(1); 
	}
			
	if(!sw_param -> filein.size())
	{
		printf("Input filename missing\n");
		exit(1);
	}

	if(sw_param -> kernel < 2)
	{
		printf("Standard BOW selected. Readjusting params\n");
		sw_param -> kernel = 2;
		sw_param -> type = 1;
	}

	
	std::cout << "[PAR] BOW input filename: "  << sw_param -> filein  << std::endl;	
	if(sw_param -> type == 1)
	{
		std::cout << "[PAR] Standard Bag-of-words -- NO GFP"  << std::endl;	
		std::cout << "[PAR] TF-IDF flavour: " << sw_param -> bow_subtype  << std::endl;	
	}
	else
	{
		std::cout << "[PAR] Geometrical FLIRT Phrases matching "  << std::endl;	
		std::cout << "[PAR] Kernel size : "  << sw_param -> kernel  << std::endl;	
	}
		
	if(sw_param -> bag == 1)
		std::cout << "[PAR] Bag of distances " << sw_param -> alpha_vss;
	if(sw_param -> type == 1 && sw_param -> bow_subtype == 2)
	std::cout << "[PAR] Type of matching: "   << sw_param -> type  << std::endl;	
	if(sw_param -> type == 1 && sw_param -> bow_subtype == 2)
		std::cout << "[PAR] alpha smoothing " << sw_param -> alpha_vss;

	std::cout << "[PAR] Kbest results : "  << sw_param -> kbest  << std::endl;	
	return(1);
}

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

int main (int argc, char **argv)
{
	sw_param_str sw_param;
	std::cout << std::endl;
	std::cout << ">> Geometrical FLIRT Phrases for Large Scale Place Recognition in 2D Range Data" << std::endl;
	std::cout << ">> G.D. Tipaldi, L.Spinello, W. Burgard -- Int. Conf. Robotics and Automation (ICRA) 2013" << std::endl;
	std::cout << ">> coded: L. Spinello 2013" << std::endl  << std::endl;
	std::cout << "** One query scan example **" << std::endl  << std::endl;

	if(!parse_command_line(argc, argv, &sw_param))
		exit(0);
	
	class gflip_engine gfp (sw_param.kernel, sw_param.kbest, sw_param.bag, sw_param.bow_subtype, sw_param.alpha_vss);

	int ret2 = gfp.read_wordscan_file(sw_param.filein);
	std::cout << "Read FLIRT word scans: " << ret2 << std::endl;
	 	
	std::cout << "Preparing inverted file index and TF-IDF" << std::endl;
	gfp.prepare( );
	
	std::cout << "Example of querying a scan " << std::endl;
	std::vector <int> query_v (26);
	int qry[] = {157, 153, 157, 160, 156,
				 158, 84,  162, 168, 162,
				 168, 224, 187, 134, 157, 
				 130 ,130, 153 ,187, 130, 
				 167, 194, 160, 130, 130,181};
	query_v.assign(&qry[0], &qry[0]+26);
	
	std::vector < std::pair <double, int> > *scorequery = NULL;
	gfp.query(sw_param.type, query_v,  &scorequery);
	for(int ii=0;ii<sw_param.kbest;ii++)
		std::cout << "{scan id: " << (*scorequery)[ii].second << " score: "<< (*scorequery)[ii].first<<"}" << std::endl;
	std::cout << "done." << std::endl; 

}

