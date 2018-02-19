
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include "DNA.h"

using namespace std;

int main(int argc, char *argv[])
{
if(argc == 1)
{
  cout << "Get the length information whitin in fasta file" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f fasta format: input file name" << endl;
}
else
{
      clock_t start,finish;
      double totaltime;

      int opt;
      int longer_len = 0, shorter_len = 0, equal_len = 0;
      char *in_genbank = 0, *out_fasta=0, *in_fasta= 0;

  while((opt = getopt(argc, argv, "f:")) != -1 )
  {
      switch(opt)
      {
      case 'f':
	  in_fasta = optarg;
	  break;
      case '?':
	  printf("unknown option: %c\n", optopt);
	  return EXIT_SUCCESS;
	  break;
      }
   }

  DNA mydna;

  if (in_fasta)
   mydna.readFasta(in_fasta);

    int dna_num = mydna.getChainNum();
    unsigned long long tmp_len, total_len = 0, N50_Total = 0, N50Size = 0;
    bool N50_set = false;
    cout << "Read " << dna_num << " sequence in."<< endl;
    cout << endl;


  /// sequence length statistics
     map<int,int> reads_len_counter;
     map<int,int>::const_iterator map_it;
     vector<int> reads_len;

    for(int i = 0; i < dna_num; i++)
	{
 	 tmp_len =  mydna.chain.at(i).length();
	 cout <<"\t" << mydna.chain.at(i).id << "\t" << tmp_len << endl;
	 reads_len.push_back(tmp_len);
	 total_len += tmp_len;
	if (reads_len_counter.count(tmp_len))
		++reads_len_counter[tmp_len];
	  else
		reads_len_counter.insert(map<int,int>::value_type(tmp_len,1));
	}

    map_it = reads_len_counter.begin();
    while ( map_it != reads_len_counter.end())
      {
      cout <<"\t"<< map_it->first << "\t" << map_it->second << endl;
      ++map_it;

      } // end while

      sort(reads_len.begin(),reads_len.end(),greater<int>());
	for (int i = 0; i < reads_len.size(); i++)
		{
		N50_Total += reads_len.at(i);
		if(N50_Total >= total_len/2)
			{
				N50Size = reads_len.at(i);
				break;
			}
		}
      cout << "\tTotal sequences " <<"\t"<< dna_num << endl;
      cout << "\tTotal Base " << "\t" << total_len << endl;
      cout << "\tLongest length" <<"\t"<< reads_len.at(0) <<endl;
      cout << "\tShortest length" <<"\t"<< reads_len.at(reads_len.size() - 1 ) <<endl;
      cout << "\tavg. length" <<"\t"<<total_len / dna_num <<endl;
      cout << "\tN50. Size" << "\t" << N50Size << endl;
}
  return EXIT_SUCCESS;
}
