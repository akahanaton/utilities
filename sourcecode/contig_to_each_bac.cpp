/***************************************************************************
 *   Copyright (C) 2008 by wenming   *
 *   wenming@genomics.org.cn   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../myHeader/DNA.h"

using namespace std;

void usage();
void DealEachBAC( char* id_file, map_DNA& seq, map<string, string>& id_info);

int main(int argc, char *argv[])
{
	if(argc == 1)
		usage();
	else
 	{
      int opt;
      char *in_seq = 0, *sanger_id_file = 0;
	  while((opt = getopt(argc, argv, "f:i:")) != -1 )
	  {
		 switch(opt)
		  {
			  case 'f':
				in_seq = optarg;
				break;
			  case 'i':
				sanger_id_file = optarg;
				break;
			  case '?':
				 printf("unknown option: %c\n", optopt);
				 return EXIT_SUCCESS;
				 break;
		  }
	   }
	 map_DNA ext_bac;
	 readFastaMap(in_seq, ext_bac);
	 map<string, string> index_id;
	 map< string, string>::iterator map_it;
	 for(map_it = ext_bac.begin(); map_it != ext_bac.end(); ++map_it)
	 {
		 int pos = map_it->first.find("le");
		 string k = map_it->first.substr(0, pos);
		 index_id.insert(pair<string, string>(k, map_it->first));
	 }
	 DealEachBAC(sanger_id_file, ext_bac, index_id);
	}// end else

  return EXIT_SUCCESS;
}


void usage()
{
  cout << "Get the solexa reads info in the gap region" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f *.fasta:  contigs sequences"<<endl;
  cout << "\t\t-b tab string: blast result file, m 8 format"<<endl;
  cout << "\t\t-s tab string: soap result file, paired alignmet hits file"<<endl;
  cout << "\t\t-o tab string: result file"<<endl;
}

void DealEachBAC( char* id_file, map_DNA& seq, map<string, string>& id_info)
{
	try
	{
		ifstream Fin(id_file);
		if(!Fin) throw strcat( "Cannot open input result file", id_file);
		string out_bac_file;
		map_DNA out_bac_seq;

		string buf, pre_bac_id, ext_index;
		int index = 1;
		int same_index_b = 1, same_index_e = 0;

		getline(Fin, pre_bac_id);

		for(;;)
		{
			getline(Fin, buf);
			if(Fin.eof())
			{
				out_bac_file="phrap_contig/" + pre_bac_id + "/extent.fa";
				same_index_e = index;
				out_bac_seq.clear();
				for(int j = same_index_b; j <= same_index_e; ++j)
				{
					string str_index = boost::lexical_cast<string>(j);
					ext_index = "contigs_" + str_index + "_";
					string k = id_info[ext_index];
					out_bac_seq.insert(pair<string, string>(k, seq[k]));
				}
				writeFastaMap(const_cast<char *>(out_bac_file.c_str()), out_bac_seq);
				break;
			}
			else
			{
				if( buf == pre_bac_id)
					++index;
				else
				{
					// process the found bac first
					out_bac_file="phrap_contig/" + pre_bac_id + "/extent.fa";
					same_index_e = index;
					out_bac_seq.clear();
					for(int j = same_index_b; j <= same_index_e; ++j)
					{
						string str_index = boost::lexical_cast<string>(j);
						ext_index = "contigs_" + str_index + "_";
						string k = id_info[ext_index];
						out_bac_seq.insert(pair<string, string>(k, seq[k]));
					}
					writeFastaMap(const_cast<char *>(out_bac_file.c_str()), out_bac_seq);
					// set the begin info for a new bac
					++index;
					same_index_b = index;
					pre_bac_id = buf;
				}
			}
		}//end for
		Fin.close();
	}//end try
	catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl;}
}
