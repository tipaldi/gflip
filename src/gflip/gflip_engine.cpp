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


// ---------------------------------------------------------

void gflip_engine::prepare(void)
{	
 	cache_binomial_coeff();

	if(bow_type==1)
		reformulate_to_bagofdistances();
	
	build_tfidf();	
	

	//~ norms
	normgfp_rc_idf_sum.resize(max_bow_len);
	normgfp_rc_weak_match.resize(max_bow_len);
	for(uint i=0;i<laserscan_bow.size();i++)
		laserscan_bow[i].norm_wgv = norm_gfp(laserscan_bow[i].w);

	//~ prepare for matching	
	mtchgfp_rc_weak_match.resize(laserscan_bow.size() * max_bow_len);
	mtchgfp_rc_idf_sum.resize(laserscan_bow.size() * max_bow_len);
	
	mtchgfp_used_doc_idx.resize(laserscan_bow.size());
	mtchgfp_max_det_idx = std::vector <int> (laserscan_bow.size());
	mtchgfp_min_det_idx = std::vector <int> (laserscan_bow.size());
}

 
// ---------------------------------------------------------
 
void gflip_engine::cache_binomial_coeff(void)
{
	cached_binomial_coeff.resize (DEFAULT_CACHEBINOMIAL,0);
	for(uint i=wgv_kernel_size-1;i<DEFAULT_CACHEBINOMIAL;i++)
		cached_binomial_coeff[i] = boost::math::binomial_coefficient <double>(i, (double)wgv_kernel_size-1);
}
 
// ---------------------------------------------------------

void gflip_engine:: reformulate_to_bagofdistances(void)
{
	int sz = (bow_dst_end-bow_dst_start)/bow_dst_interval+1;

	for(uint i=0;i<laserscan_bow.size();i++)
	{
		std::vector <int> w_tmp;
		laserscan_bow[i].w.clear();
		for(int j=0;j<(int)laserscan_bow[i].w_x.size()-1;j++)
		{
			double dx = laserscan_bow[i].w_x[j] - laserscan_bow[i].w_x[j+1];
			double dy = laserscan_bow[i].w_y[j] - laserscan_bow[i].w_y[j+1];
			double d = sqrt(dx*dx + dy*dy);
			int idx = (d-bow_dst_start)/bow_dst_interval;
			if(idx > sz)
				idx = sz +1;
			w_tmp.push_back(idx);
			laserscan_bow[i].w = w_tmp;
			//~ std::cout << d << " " << dx << " "<< dy << " "<< idx << " " << sz << std::endl;
		}
	}
		
}

// ---------------------------------------------------------

double gflip_engine::norm_gfp(std::vector <int> & query_v)
{
	double norm2 = 0,query_v_norm=1;
	std::fill(normgfp_rc_idf_sum.begin(), normgfp_rc_idf_sum.end(), 0);
 	std::fill(normgfp_rc_weak_match.begin(), normgfp_rc_weak_match.end(), 0);

	int min_det_idx_qry=INT_MAX;
	int max_det_idx_qry=-INT_MAX;
	int middleidx = max_bow_len/2;

	for(uint j=0;j<query_v.size();j++)
	{
		int word_id = query_v[j];
		for(uint h=0;h<query_v.size();h++)
		{
			if(word_id == query_v[h])
			{
				int w_order_dif=j-h;
				normgfp_rc_weak_match[middleidx+w_order_dif]++;
				normgfp_rc_idf_sum[middleidx+w_order_dif] += tf_idf [word_id].idf;
				
				if(middleidx+w_order_dif < min_det_idx_qry)
					min_det_idx_qry = middleidx+w_order_dif;
				if(middleidx+w_order_dif > max_det_idx_qry)
					max_det_idx_qry = middleidx+w_order_dif;		
			}
		}
	}

	for(int b=min_det_idx_qry;b<=max_det_idx_qry;b++)
	{
		double combo = 0;
		if(normgfp_rc_weak_match[b] >= wgv_kernel_size )
			combo = cached_binomial_coeff[normgfp_rc_weak_match[b]-1];
		norm2 +=normgfp_rc_idf_sum[b] * combo;
	}
	//~ rounding errs
	if(norm2 > 0)
		query_v_norm = sqrt(norm2); 

	return(query_v_norm);
}

 
 
// ---------------------------------------------------------

void gflip_engine::matching_gfp(std::vector <int> &query_v)
{
	int middleidx = max_bow_len/2;
	std::fill(mtchgfp_used_doc_idx.begin(), mtchgfp_used_doc_idx.end(), 0);
	std::fill(mtchgfp_min_det_idx.begin(), mtchgfp_min_det_idx.end(), +INT_MAX);
	std::fill(mtchgfp_max_det_idx.begin(), mtchgfp_max_det_idx.end(), -INT_MAX);
 
 	std::fill(mtchgfp_rc_weak_match.begin(), mtchgfp_rc_weak_match.end(), 0);
 	std::fill(mtchgfp_rc_idf_sum.begin(), mtchgfp_rc_idf_sum.end(), 0);

	//~ query norm
	double query_v_norm = norm_gfp(query_v);
 	//~ every word of the query
	for(uint j=0;j<query_v.size();j++)
	{
		int word_id = query_v[j];
		for(uint a=0;a<tf_idf[word_id].doc_id.size();a++)
		{
			int doc_idx = tf_idf[word_id].doc_id[a];
			mtchgfp_used_doc_idx[doc_idx] = 1;

			//~ match positions
			for(uint h=0;h<tf_idf[word_id].word_order[a].pos.size();h++)
			{
				int w_order_dif=j-tf_idf[word_id].word_order[a].pos[h];
				int rcidx = (max_bow_len * doc_idx) + (middleidx+w_order_dif);
							
				mtchgfp_rc_weak_match[rcidx]++;
				mtchgfp_rc_idf_sum[rcidx] += tf_idf [word_id].idf;
				
				if(middleidx+w_order_dif < mtchgfp_min_det_idx[doc_idx])
					mtchgfp_min_det_idx[doc_idx] = middleidx+w_order_dif;
				if(middleidx+w_order_dif > mtchgfp_max_det_idx[doc_idx])
					mtchgfp_max_det_idx[doc_idx] = middleidx+w_order_dif;								
			}
		}
	}

	int num_used_doc_idx = 0;
	for(uint j=0;j<mtchgfp_used_doc_idx.size();j++)
		if(mtchgfp_used_doc_idx[j] == 1)
			num_used_doc_idx++;
	scoreset.resize(num_used_doc_idx);
	for(uint j=0, u_idx=0;j<mtchgfp_used_doc_idx.size();j++)
	{
		if(mtchgfp_used_doc_idx[j] == 0)
			continue;

		//~ default values
		int doc_idx = j;
		double score = 0;
		scoreset[u_idx].first = 1.0;
		scoreset[u_idx].second = doc_idx;
		
		//~ compute score
		for(int b=mtchgfp_min_det_idx[doc_idx];b<=mtchgfp_max_det_idx[doc_idx];b++)
		{
			double combo = 0;
			int rcidx = (max_bow_len * doc_idx) + b;
			 
			if(mtchgfp_rc_weak_match[rcidx] >= wgv_kernel_size )
				combo = cached_binomial_coeff[ mtchgfp_rc_weak_match[rcidx]-1 ];
			score +=mtchgfp_rc_idf_sum[rcidx] * combo;
		}		
		//~ normed istance
		score = score / (laserscan_bow[doc_idx].norm_wgv * query_v_norm);
 		
		//~ avoids no go zone
		if( doc_idx <= start_l || doc_idx >= stop_l)
			scoreset[u_idx].first = 1.0 - score;
		
		u_idx++;	
  	}
	
	//~ sort the results
 	sort(scoreset.begin(), scoreset.end(), isBettermatched);
	 
}

// ---------------------------------------------------------

bool isBettermatched(std::pair <double, int> x, std::pair <double, int> y) 
{
    return x.first < y.first;
}

// ---------------------------------------------------------

void gflip_engine::matching_bow(std::vector <int> &query_v )
{
	std::vector <double> image_db_scores(laserscan_bow.size(),0);
	std::set<int> used_doc_idx;
	
	double query_v_norm = 1, qsum = 0;
	std::vector <double> query_bow = std::vector <double> (dictionary_dimensions,0);
	for(uint j=0;j<query_v.size();j++)
		query_bow[query_v[j]]++;

	for(uint j=0;j<query_bow.size();j++)
		qsum+=query_bow[j]*query_bow[j];
	if(qsum != 0)
		query_v_norm = sqrt(qsum);

	//~ query tdidf voting
	for(uint j=0;j<query_v.size();j++)
	{
		int word_id = query_v[j];
		for(uint a=0;a<tf_idf[word_id].doc_id.size();a++)
		{
			int img_idx = tf_idf[word_id].doc_id[a];
			
			if(bow_subtype == 0)
				image_db_scores[img_idx] += tf_idf[word_id].tf_idf_doc_normed[a];
			if(bow_subtype == 1)
				image_db_scores[img_idx] += tf_idf[word_id].wf_idf_doc_normed[a];
			if(bow_subtype == 2)
				image_db_scores[img_idx] += tf_idf[word_id].ntf_idf_doc_normed[a];

			used_doc_idx.insert(img_idx);
		}
	}

	scoreset.resize(used_doc_idx.size());
	int u_idx=0;
	for (std::set<int>::iterator it=used_doc_idx.begin(); it!=used_doc_idx.end(); it++)
	{
		double score = image_db_scores[*it]/query_v_norm;
		
		scoreset[u_idx].first = 1.0;
		scoreset[u_idx].second = *it;

		//~ avoids no-go zone
		if( *it <= start_l || *it >= stop_l)
			scoreset[u_idx].first = 1.0 - score;
		u_idx++;
	}
 

	sort(scoreset.begin(), scoreset.end(), isBettermatched);
 	 
 }


// ---------------------------------------------------------

 
void gflip_engine::query(int dtype, std::vector <int> &query_v, std::vector < std::pair <double, int> > **scoreoutput)
{
	//~ avoids any skip 
	start_l = 0; 
	stop_l = 0;
	//~ does the search
	if(dtype ==1)
		matching_bow(query_v);
	if(dtype ==2)
		matching_gfp(query_v);
		
	*scoreoutput = &scoreset;	
}

// ---------------------------------------------------------



void gflip_engine::run_evaluation(int dtype)
{
	struct timeval tim_st,tim_ed;  
	char nosave = 0;

	char buff[2000];
	sprintf(buff,"./%s.nn", fileoutput_rootname.c_str());
	FILE *f;
	
	if(!nosave)
		f=fopen(buff, "wt");
		
	double dtime_avg=0;
	int countv=0;
	for(uint i=0;i<number_of_scans;i++)
	{
		//~ match this query,just the seq found
		std::vector <int> query_v  = laserscan_bow[i].w;

		//~ if 0 len
		if(!query_v.size())
		{
			if(!nosave)
				fprintf (f, "\n");
			continue;
		}
		//~ exclude the i scan
		start_l = i-1; 
		stop_l = i+1;

		//~ match with several techs
    	gettimeofday(&tim_st, NULL);  
		if(dtype ==1)
			matching_bow(query_v);
		if(dtype ==2)
			matching_gfp(query_v);
		gettimeofday(&tim_ed, NULL);  
		double dtimeqry = (tim_ed.tv_sec-tim_st.tv_sec) + (tim_ed.tv_usec-tim_st.tv_usec)/1000000.0;
		dtime_avg += dtimeqry;
		if(!nosave)
		{
			fprintf (f, "%d %d %7.7f ", kbest, (int)scoreset.size(), dtimeqry);
			uint minscoresetsize = std::min((int)kbest,(int)scoreset.size());
			for(uint ii=0;ii<minscoresetsize;ii++)
				fprintf (f, "%d ", scoreset[ii].second);
			fprintf (f, "\n");
		}
		countv++;
	}
	std::cout << "Averge query time: "<< (double)dtime_avg/(double)countv << " total  time: " << dtime_avg << " # scans: " << countv << std::endl;
	if(!nosave)
		fclose(f);
}

 
 

// ---------------------------------------------------------


void gflip_engine::build_tfidf(void)
{
	//~ find id size, maxlen
	int maxid= -1;
	int maxid_idx = -1;
	max_bow_len = -INT_MAX;

	for(uint i=0;i<laserscan_bow.size();i++)
	{
		for(uint j=0;j<laserscan_bow[i].w.size();j++)
		{
			if(laserscan_bow[i].w[j] > maxid)
			{
				maxid = laserscan_bow[i].w[j];
				maxid_idx = i;
			}
		
		}
			
		if((int)laserscan_bow[i].w.size() > max_bow_len)
			max_bow_len=laserscan_bow[i].w.size();			
	}
	//~ include last number
	maxid+=1;
	dictionary_dimensions = maxid;
	//~ do it large
	max_bow_len = (max_bow_len+1)*2;
	std::cout << "Detected dictionary dimension: "<< maxid << " @ "<< maxid_idx << std::endl;
	std::cout << "Detected max bow len : "<< max_bow_len << std::endl;

	tf_idf = std::vector <tf_idf_db> (maxid);
	for(int word_id=0; word_id<maxid; word_id++)
	{
		//~ doc id
		for(uint i=0;i<laserscan_bow.size();i++)
		{
			int term_count_unnormalized=0;
			tf_idf_db_ordercache w_order;
			for(uint j=0;j<laserscan_bow[i].w.size();j++)
			{
				if(laserscan_bow[i].w[j] == word_id)
				{
					term_count_unnormalized++;
					w_order.pos.push_back(j);
				}
			}
			//~ found a word match
			if(term_count_unnormalized)
			{
				tf_idf[word_id].term_count_unnormalized.push_back(term_count_unnormalized);
				tf_idf[word_id].num_words.push_back(laserscan_bow[i].w.size());
				tf_idf[word_id].doc_id.push_back(i);
				tf_idf[word_id].term_count.push_back((double)term_count_unnormalized / (double)laserscan_bow[i].w.size());
				tf_idf[word_id].word_order.push_back(w_order);
				tf_idf[word_id].tf_idf_doc_normed.push_back(-1);
				tf_idf[word_id].wf_idf_doc_normed.push_back(-1);
				tf_idf[word_id].ntf_idf_doc_normed.push_back(-1);
			}
		}
		tf_idf[word_id].num_doc_containing_the_word = tf_idf[word_id].doc_id.size();
		tf_idf[word_id].corpus_size = laserscan_bow.size();
		tf_idf[word_id].idf = log( (double)tf_idf[word_id].corpus_size / (double) tf_idf[word_id].num_doc_containing_the_word );
	}
	
	//~ tfsmoothing
	std::vector < double > mxtf_val (laserscan_bow.size(),-DBL_MAX);
	for(int doc_id=0;doc_id<(int)laserscan_bow.size();doc_id++)
	{
		std::set<int> used_idx;
		for(uint j=0;j<laserscan_bow[doc_id].w.size();j++)
		{
			int word_id = laserscan_bow[doc_id].w[j];
			for(uint h=0; h<tf_idf[word_id].doc_id.size(); h++)
				if(tf_idf[word_id].doc_id[h] == doc_id)
					if( tf_idf[word_id].term_count_unnormalized[h] > mxtf_val [doc_id] )
						mxtf_val[doc_id] = tf_idf[word_id].term_count_unnormalized[h];
		}
		//~ printf("%d %g %d\n",doc_id,mxtf_val,mxtf_idx);
	}
	
	//~ normedtfidf, wtfidf
	for(int doc_id=0;doc_id<(int)laserscan_bow.size();doc_id++)
	{
		std::set<int> used_idx;
		for(uint j=0;j<laserscan_bow[doc_id].w.size();j++)
			used_idx.insert(laserscan_bow[doc_id].w[j]);
		//~ sum
		double sum=0,sum_wf=0,sum_vss=0;
		for (std::set<int>::iterator word_id_iter=used_idx.begin(); word_id_iter!=used_idx.end(); word_id_iter++)
			for(uint h=0;h<tf_idf[*word_id_iter].doc_id.size();h++)
				if(tf_idf[*word_id_iter].doc_id[h] == doc_id)
				{
					double val = (tf_idf[*word_id_iter].term_count[h] * tf_idf[*word_id_iter].idf);
					double val_wf = (1 + log(tf_idf[*word_id_iter].term_count_unnormalized[h]) ) * tf_idf[*word_id_iter].idf;
					double val_vss = alpha_vss + ( (1.0 - alpha_vss) * tf_idf[*word_id_iter].term_count_unnormalized[h] ) / mxtf_val[doc_id];
					sum+=val*val;
					sum_wf+=val_wf*val_wf;
					sum_vss+=val_vss*val_vss;
				}
					
		//~ norm
		double norm=sqrt(sum);
		double norm_wf=sqrt(sum_wf);
		double norm_vss=sqrt(sum_vss);
		double versum=0, versum_wf=0, versum_vss=0;
		for (std::set<int>::iterator word_id_iter=used_idx.begin(); word_id_iter!=used_idx.end(); word_id_iter++)
			for(uint h=0;h<tf_idf[*word_id_iter].doc_id.size();h++)
				if(tf_idf[*word_id_iter].doc_id[h] == doc_id)
				{
					tf_idf[*word_id_iter].tf_idf_doc_normed[h]=(tf_idf[*word_id_iter].term_count[h] * tf_idf[*word_id_iter].idf)/norm;
					tf_idf[*word_id_iter].wf_idf_doc_normed[h]=((1 + log(tf_idf[*word_id_iter].term_count_unnormalized[h]) )  * tf_idf[*word_id_iter].idf)/norm_wf;
					tf_idf[*word_id_iter].ntf_idf_doc_normed[h]=( alpha_vss + ( (1.0 - alpha_vss) * tf_idf[*word_id_iter].term_count_unnormalized[h] ) / mxtf_val[doc_id]) / norm_vss;
					
					versum += tf_idf[*word_id_iter].tf_idf_doc_normed[h]*tf_idf[*word_id_iter].tf_idf_doc_normed[h];
					versum_wf += tf_idf[*word_id_iter].wf_idf_doc_normed[h]*tf_idf[*word_id_iter].wf_idf_doc_normed[h];
					versum_vss += tf_idf[*word_id_iter].ntf_idf_doc_normed[h]*tf_idf[*word_id_iter].ntf_idf_doc_normed[h];
				}
	
		//~ verification			
		if( fabs(sqrt(versum) -1 ) > 0.00001 && sum > 0.00001 )
		{
			std::cout << "ERROR NORMALIZ FAIL "<<sqrt(versum)<< " "<< doc_id<< " "<< laserscan_bow[doc_id].w.size() << std::endl;
			exit(1);
		}
		
		if( fabs(sqrt(versum_wf) -1 ) > 0.00001 && sum_wf > 0.00001 )
		{
			std::cout << "ERROR WFIDF NORMALIZ FAIL "<<sqrt(versum_wf)<< " "<< doc_id<< " "<< laserscan_bow[doc_id].w.size() << std::endl;
			exit(1);
		}

		if( fabs(sqrt(versum_vss) -1 ) > 0.00001 && sum_vss > 0.00001 )
		{
			std::cout << "ERROR VSSIDF NORMALIZ FAIL "<<sqrt(versum_vss)<< " "<< doc_id<< " "<< laserscan_bow[doc_id].w.size() << std::endl;
			exit(1);
		}

	}

	//~ verification
	for(uint i=0;i<tf_idf.size();i++)
		for(uint h=0;h<tf_idf[i].doc_id.size();h++)
		{
			if(tf_idf[i].tf_idf_doc_normed[h] < 0)
			{
				std::cout << "ERROR NORMALIZ NO INT"<<std::endl;
				exit(1);
			}
			if(tf_idf[i].wf_idf_doc_normed[h] < 0)
			{
				std::cout << "ERROR wf_idf NO INT"<< i << " "<< h << " "<< tf_idf[i].wf_idf_doc_normed[h] << " "<<  std::endl;
				exit(1);
			}
			if(tf_idf[i].ntf_idf_doc_normed[h] < 0)
			{
				std::cout << "ERROR vss NO INT"<< i << " "<< h << " "<< tf_idf[i].ntf_idf_doc_normed[h] << " "<<  std::endl;
				exit(1);
			}
			
		}
}

// ---------------------------------------------------------


void gflip_engine::insert_wordscan(std::vector <int> scanbow, std::vector <double> xpos, std::vector <double> ypos)
{
	uint numwords = scanbow.size();
	scan_bow tmpbow(scanbow.size());
 	
 	for(uint i=0;i<numwords;i++)
	{
		tmpbow.w[i] = scanbow[i];
		tmpbow.w_x[i] = xpos[i];
		tmpbow.w_y[i] = ypos[i];
	}
	laserscan_bow.push_back(tmpbow);
	number_of_scans = laserscan_bow.size();

}
// ---------------------------------------------------------

int gflip_engine::read_wordscan_file(std::string filename)
{
	std::ifstream ifs(filename.c_str());
	std::string line;

	//~ getfilename
	std::vector<std::string> ftokens;
	LSL_stringtoken(filename, ftokens, "/");
	fileoutput_rootname = ftokens[ftokens.size()-1];
	//~ std::cout<< fileoutput_rootname <<std::endl;

    int count=0;
	// index it
	while( std::getline( ifs, line ) )
	{
		std::vector <std::string> tokens;
		LSL_stringtoken(line, tokens, " ");  
	 
		if(tokens.size())
		{
			uint numwords=atoi(tokens[0].c_str());
			scan_bow tmpbow(numwords);			
			
			if(tokens.size() != 1+numwords*3)
			{
				std::cout << "Error: input file must contain coordinates" << std::endl;
				exit(1);
			}
			
 			for(uint i=0, a=0;i<numwords*3;i+=3,a++)
			{
				tmpbow.w[a] = atoi(tokens[i+1].c_str());
				tmpbow.w_x[a] = atof(tokens[i+1+1].c_str());
				tmpbow.w_y[a] = atof(tokens[i+2+1].c_str());
				//~ printf("%d ",tmpbow.w[a]);
 			}
			laserscan_bow.push_back(tmpbow);
			count++;
 
		}
	}

	//~ set num scans
	number_of_scans = laserscan_bow.size();

	return(count);
}
  
// ---------------------------------------------------------

void LSL_stringtoken(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters)
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
 
