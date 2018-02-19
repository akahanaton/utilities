
#include <set>
#include <unistd.h>
#include <sys/stat.h>

#include "../myHeader/DNA.h"
#include "../myHeader/blastparser.h"

using namespace std;

typedef map<string, vector<single_gap> > map_GAPS;

void usage();
void DealSoapResualt( char* in_file, map_GAPS& gaps, map_DNA& reads_a_, map_DNA& reads_b_);
int binary_pos(vector<single_gap>& gap_regions, const int& regions_num, int& key_pos );
void GetGapInScaf( map_GAPS &gap_regions, map_DNA &dna, const int &insert_size);
void OutputGapinfo(map_GAPS& gaps, map_DNA &dna);

void ConnectContigAccording2PE();

int main(int argc, char *argv[])
{
	if(argc == 1)
		usage();
	else
	{
      clock_t start,finish;
      double totaltime;
      start=clock();

      int opt, in_size = 400;
      char *scaffold_file = 0,*reads_a_file = 0, *reads_b_file = 0, *soap_file = 0, *out_folder = 0, *parameter_file = 0;
	  while((opt = getopt(argc, argv, "f:c:a:b:2:i:o:p:")) != -1 )
	  {
		 switch(opt)
		  {
			  case 'f':
				scaffold_file = optarg;
				break;
			  case 'a':
				reads_a_file = optarg;
				break;
			  case 'b':
				reads_b_file = optarg;
				break;
			  case '2':
				soap_file = optarg;
				break;
			  case 'o':
				out_folder = optarg;
				break;
			  case 'p':
				parameter_file = optarg;
				break;
			  case 'i':
				in_size = boost::lexical_cast<int>(optarg);
				break;
			  case '?':
				 printf("unknown option: %c\n", optopt);
				 return EXIT_SUCCESS;
				 break;
		  }
	  }

		map_DNA scaffold;
		map_DNA  reads_a, reads_b;
		//--------------------------------------------------
		// scaffold.readFasta(scaffold_file);
		// cout << scaffold.getChainNum() << " scaffold read in" << endl;
		//--------------------------------------------------
		readFastaMap(scaffold_file, scaffold);
		cout << scaffold.size() << " scaffold read in" << endl;

		map_GAPS gaps_in_scaffold;
		GetGapInScaf(gaps_in_scaffold, scaffold, in_size);

		if(parameter_file)
		{
			FILE * pFile;
			char buf[4096];
			char *pch;

			pFile = fopen (parameter_file, "r");
			if (pFile == NULL) perror ("Error opening file");
			else {
				while( fgets(buf,4096,pFile )!= NULL)
				{
					pch = strtok (buf," ");
					pch = strtok (NULL, " ");
					reads_a_file = pch;
					pch = strtok (NULL, " ");
					pch = strtok (NULL, " ");
					reads_b_file = pch;
					pch = strtok (NULL, " ");
					pch = strtok (NULL, " \n");
					soap_file = pch;
					cout << reads_a_file << "\t" << reads_b_file << "\t" << soap_file << endl;
					reads_a.clear();
					reads_b.clear();
					readFastqMap(reads_a_file, reads_a);
					cout << reads_a.size() << " reads read in" << endl;
					readFastqMap(reads_b_file, reads_b);
					cout << reads_b.size() << " reads read in" << endl;
					DealSoapResualt(soap_file, gaps_in_scaffold, reads_a, reads_b);
				}

				fclose (pFile);
			}

		}
		else if(reads_a_file && reads_b_file)
		{
			readFastqMap(reads_a_file, reads_a);
			cout << reads_a.size() << " reads read in" << endl;
			readFastqMap(reads_b_file, reads_b);
			cout << reads_b.size() << " reads read in" << endl;
			DealSoapResualt(soap_file, gaps_in_scaffold, reads_a, reads_b);
		}

		OutputGapinfo(gaps_in_scaffold, scaffold);


	}// end else

  return EXIT_SUCCESS;
}

void usage()
{
  cout << "Get the solexa reads info in the gap region" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f STR: scaffold sequences, fasta file"<<endl;
  cout << "\t\t-a STR: reads a sequences, fastq file"<<endl;
  cout << "\t\t-b STR: reads b sequences, fastq file"<<endl;
  cout << "\t\t-i INT: paired end insert size ,default 400"<<endl;
  cout << "\t\t-2 STR: soap result file, single aligned hits file"<<endl;
  cout << "\t\t-o STR: output folder name" << endl;
  cout << "\t\t-p STR: parameter file" << endl;
  cout << "\t\t\t: -a fq1.file -b fq2.file -2 soap.single.out" << endl;
}

void OutputGapinfo(map_GAPS& gaps ,map_DNA &dna)
{
	ofstream Fout("gap.info");

	if(mkdir("reads_in_gap", S_IRWXU | S_IRWXG | S_IRWXO ) < 0)
		cout << "mkdir error " << "reads_in_gap"  << endl;
	if(chdir("reads_in_gap") < 0)
		cout << "chdir error " << "reads_in_gap "<< endl;
	else
		cout << "chdir ok" << endl;

	string out_seq_file;
	map<string, vector< single_gap > >::iterator map_it = gaps.begin();
	string tmp_str;
	while (map_it != gaps.end())
	{
		Fout <<">"<< map_it->first << "\t" << dna[map_it->first].length() << endl;
		int gap_num = map_it->second.size();
		for(int i = 0; i < gap_num; ++i)
		{
			Fout << "# " <<  i << "\t"
				 << map_it->second.at(i).real_gap_size << "\t"
				 << map_it->second.at(i).gap1_start_pos << "\t"
				 << map_it->second.at(i).gap1_end_pos << "\t"
				 << map_it->second.at(i).gap2_start_pos << "\t"
				 << map_it->second.at(i).gap2_end_pos << "\t"
				 << map_it->second.at(i).long_reads.getChainNum() << "\t"
				 << map_it->second.at(i).short_reads.getChainNum() << endl;

			tmp_str = map_it->first + "_" +  boost::lexical_cast<string>(i) + "_" + boost::lexical_cast<string>( map_it->second.at(i).real_gap_size );
			out_seq_file = tmp_str + "_l.seq";
			map_it->second.at(i).long_reads.writeFasta(const_cast<char*>(out_seq_file.c_str()));
			out_seq_file = tmp_str + "_s.seq";
			map_it->second.at(i).short_reads.writeFasta(const_cast<char*>(out_seq_file.c_str()));
			int short_seq_num = map_it->second.at(i).short_reads.getChainNum();
			for ( int j = 0; j < short_seq_num; j++)
				Fout <<"\t" << map_it->second.at(i).short_reads.chain.at(j).description << endl;
		}
		++map_it;
	}
	Fout.close();
}

void DealSoapResualt( char* in_file, map_GAPS& gaps, map_DNA& reads_a_, map_DNA& reads_b_)
{
	try
	{
		ifstream Fin(in_file);
		if(!Fin) throw strcat( "Cannot open input result file", in_file);

		fastaDNA tmp_chain;
		DNA tmp_dna;

		string buf, scaf_name, read_name, read_index, tmp_read_id;
		int tmp_pos;
		string tmp_pos_str;
		int pos_1 = 0, pos_2 = 0;
		short int rcom = 1;
		int read_len = 0;
		int hit_num = 0;

		map<string, string> read_index_swither;
		read_index_swither.insert(pair<string, string>("1", "2"));
		read_index_swither.insert(pair<string, string>("2", "1"));

		vector<single_gap>::iterator it1, it2;

		for(;;)
		{
			getline(Fin, buf);
			if(Fin.eof())
				break;
			else
			{
				//soap format: 44 \t +/- \t scaffold
				pos_1 = buf.find("caffold");
				if( pos_1!= string::npos )
				{
					if(buf.substr(pos_1-2, 1) == "-")
						rcom = 1;
					else
						rcom = 0;
					//--------------------------------------------------
					// read_len = boost::lexical_cast<int>( buf.substr( pos_1-6,2 ) )
					//--------------------------------------------------
					hit_num = boost::lexical_cast<int>( buf.substr( pos_1-10,1 ) );
					pos_2 = buf.find_first_of("\t", pos_1);
					scaf_name = buf.substr(pos_1 -1 , pos_2 - pos_1 + 1);
					boost::to_lower(scaf_name);
					pos_1 = buf.find_first_of("\t", pos_2 + 1);
					tmp_pos_str =  buf.substr(pos_2+1, pos_1 - pos_2);
					tmp_pos = boost::lexical_cast<int>( tmp_pos_str );

					pos_1 = buf.find_first_of("\t");
					read_name = buf.substr(0, pos_1 - 1);
					read_index = buf.substr(pos_1 - 1, 1);
					//--------------------------------------------------
					// cout << scaf_name <<"\t" <<tmp_pos  << "\t" << rcom << endl;
					//--------------------------------------------------
					if (hit_num == 1)
					{
						if((int)gaps.count(scaf_name) > 0)
						{
							it1 = gaps[scaf_name].begin();
							it2 = gaps[scaf_name].end();
							for ( ; it1 != it2; ++it1)
							{
								if((tmp_pos>=(*it1).gap1_start_pos) && (tmp_pos <= (*it1).gap1_end_pos) && rcom == 1 ||
								   (tmp_pos>=(*it1).gap2_start_pos) && (tmp_pos <= (*it1).gap2_end_pos) && rcom == 0 )
								{
									tmp_read_id = read_name + read_index_swither[read_index];
									if(read_index == "1")
										tmp_chain.seq = reads_b_[tmp_read_id];
									if(read_index == "2")
										tmp_chain.seq = reads_a_[tmp_read_id];
									tmp_chain.description = tmp_read_id + " " + tmp_pos_str;
									(*it1).short_reads.addChain(tmp_chain);
									break;
								}
							} // end for
						}
					}
				 }
			}// end else
		}//end for(;;)
		Fin.close();
	}//end try
	catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl;}
}


void GetGapInScaf( map_GAPS& gap_regions, map_DNA& dna, const int& insert_size)
{

	map<string, string> :: iterator map_dna_it = dna.begin();

	vector<single_gap> vet_gaps;

	for(; map_dna_it != dna.end(); ++map_dna_it)
	{
		int  cur_gap_begin = 0, next_gap_begin = 0,  gap_end = 0, real_gap_size = 0;
		fastaDNA tmp_chain;

		string cur_scaffold_id = map_dna_it->first;
		int cur_scaffold_len = map_dna_it->second.length();
		vet_gaps.clear();
		next_gap_begin = map_dna_it->second.find_first_of("Nn", cur_gap_begin);
		while (next_gap_begin != string::npos)
		{
			cur_gap_begin = next_gap_begin;
			single_gap tmp_single_gap;
			if((cur_gap_begin - gap_end) < 100)
			{
				tmp_chain.description = cur_scaffold_id + "_" + boost::lexical_cast<string>(gap_end) + "_" + boost::lexical_cast<string>(cur_gap_begin);
				tmp_chain.seq = map_dna_it->second.substr(gap_end+1, cur_gap_begin - gap_end -1);
				tmp_single_gap.long_reads.addChain(tmp_chain);
			}
			else
			{
				// 100bp befor cur_gap_begin
				tmp_chain.description = cur_scaffold_id + "_" + boost::lexical_cast<string>(cur_gap_begin - 100) + "_" + boost::lexical_cast<string>(cur_gap_begin);
				tmp_chain.seq = map_dna_it->second.substr(cur_gap_begin - 100, 100);
				tmp_single_gap.long_reads.addChain(tmp_chain);
			}

			gap_end = map_dna_it->second.find_first_not_of("Nn",cur_gap_begin);
			next_gap_begin = map_dna_it->second.find_first_of("Nn", gap_end+1);
			if( next_gap_begin != string::npos )
			{
				if( gap_end + 100 < next_gap_begin)
				{
					tmp_chain.description = cur_scaffold_id + "_" + boost::lexical_cast<string>(gap_end) + "_" + boost::lexical_cast<string>(gap_end + 100);
					tmp_chain.seq = map_dna_it->second.substr(gap_end, 100);
					tmp_single_gap.long_reads.addChain(tmp_chain);
				}
				else
				{
					tmp_chain.description = cur_scaffold_id + "_" + boost::lexical_cast<string>(gap_end) + "_" + boost::lexical_cast<string>(next_gap_begin);
					tmp_chain.seq = map_dna_it->second.substr(gap_end, next_gap_begin  - gap_end -1);
					tmp_single_gap.long_reads.addChain(tmp_chain);
				}
			}
			else
			{
				if( gap_end + 100 < cur_scaffold_len )
				{
					tmp_chain.description = cur_scaffold_id + "_" + boost::lexical_cast<string>(gap_end) + "_" + boost::lexical_cast<string>(gap_end + 100);
					tmp_chain.seq = map_dna_it->second.substr(gap_end, 100);
					tmp_single_gap.long_reads.addChain(tmp_chain);
				}
				else
				{
					tmp_chain.description = cur_scaffold_id + "_" + boost::lexical_cast<string>(gap_end) + "_" + boost::lexical_cast<string>(cur_scaffold_len);
					tmp_chain.seq = map_dna_it->second.substr(gap_end, cur_scaffold_len  - gap_end -1);
					tmp_single_gap.long_reads.addChain(tmp_chain);
				}
			}

			if(cur_gap_begin - insert_size -100 > 0)
			{
				tmp_single_gap.gap1_start_pos = cur_gap_begin - insert_size -100;
				tmp_single_gap.gap1_end_pos = gap_end - insert_size +100;

				if( gap_end + insert_size + 100 < cur_scaffold_len )
				{
					tmp_single_gap.gap2_start_pos = cur_gap_begin + insert_size - 100;
					tmp_single_gap.gap2_end_pos = gap_end + insert_size + 100;
				}
			}
			else
			{
				if( gap_end + insert_size + 100 < cur_scaffold_len )
				{
					tmp_single_gap.gap2_start_pos = cur_gap_begin + insert_size - 100;
					tmp_single_gap.gap2_end_pos = gap_end + insert_size + 100;
				}
			}
			tmp_single_gap.real_gap_size = gap_end - cur_gap_begin;
			vet_gaps.push_back(tmp_single_gap);
		}// endl while
		if (vet_gaps.size() > 0)
			gap_regions.insert(pair<string, vector<single_gap> >( cur_scaffold_id, vet_gaps ) );
	} // end for index
}
