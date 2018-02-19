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

#include "/share/raid5/wenming/myProgram/myHeader/DNA.h"
using namespace std;

struct contig
{
	string 	contig_id;
	int 	contig_pos;
	string	recomplement;
};

typedef map<string, string> FastaSeq;

void 	connectcontig(vector<contig>& input_scaffold_info, const string& input_scaffold_id, DNA& all_scaffold, FastaSeq& input_contig );
void 	ReadScafInfo(char* filename, DNA& out_scaffold, FastaSeq& contigs);



int main(int argc, char *argv[])
{
if(argc == 1)
{
  cout << "Use The contigs to build scaffold accoding the scaf info" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f fasta format: input file name, contains the contig" << endl;
  cout << "\t\t-s scaf info, according the map2contig resualt" << endl;
  cout << "\t\t-o fasta format: output file name" << endl;
}
else
{
      clock_t start,finish;
      double totaltime;
      start=clock();

      int opt;
      int longer_len = 0, shorter_len = 0, equal_len = 0;
      char *in_fasta = 0, *out_fasta=0, *in_scaf= 0;

  while((opt = getopt(argc, argv, "f:s:o:")) != -1 )
  {
      switch(opt)
      {
      case 'f':
	 	in_fasta = optarg;
	  break;
      case 's':
	 	in_scaf = optarg;
	  break;
      case 'o':
	 	out_fasta = optarg;
	  break;
      case '?':
	 	 printf("unknown option: %c\n", optopt);
	  return EXIT_SUCCESS;
	  break;
      }
   }

  DNA mydna, scaffold;
    mydna.readFasta(in_fasta);
	FastaSeq all_contigs;
   size_t seq_num = mydna.getChainNum();
   for( size_t i = 0; i < seq_num; ++i)
   {
	all_contigs.insert( map<string, string>::value_type( mydna.chain.at(i).id, mydna.chain.at(i).seq ) );
   }
	cout << "test 1" << endl;
   ReadScafInfo( in_scaf, scaffold, all_contigs );

	cout << "test 2" << endl;
   scaffold.writeFasta( out_fasta);

	finish=clock();
   totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
   cout<<"\nProgram takes "<<totaltime<<" seconds"<<endl;
}

  return EXIT_SUCCESS;
}


void connectcontig(vector<contig>& input_scaffold_info, const string& input_scaffold_id, DNA& all_scaffold,  FastaSeq& input_contig )
{
	vector<contig>::iterator contig_it( input_scaffold_info.begin() );
	 cout << "test 5" << endl;
	fastaDNA tmp_fasta_dna;
	tmp_fasta_dna.set_id(input_scaffold_id);
	tmp_fasta_dna.set_description( input_scaffold_id );

	string tmp_seq;

	for( ; contig_it != input_scaffold_info.end(); ++contig_it)
	{
		string tmp_gap(contig_it->contig_pos - tmp_fasta_dna.seq.length(), 'X');
		tmp_fasta_dna.seq.append(tmp_gap);
		if( contig_it->recomplement == "-" )
		{
			reverse_copy( input_contig[contig_it->contig_id].begin(), input_contig[contig_it->contig_id].end(), tmp_seq.begin() );
			transform( tmp_seq.begin(), tmp_seq.end(), input_contig[contig_it->contig_id].begin(), compliment_nucleotide );
		}
		tmp_fasta_dna.seq.append( input_contig[contig_it->contig_id] );
	}
	cout <<" test 4" << endl;
	all_scaffold.addChain(tmp_fasta_dna);
}

void ReadScafInfo(char* filename,  DNA& out_scaffold, FastaSeq& input_contig)
{
    try
    {
	ifstream Fin(filename);
		if(!Fin) throw strcat("Cannot open input scaf file ", filename);

		int 		newLine;
		size_t 		tmp_pos;
		string		buf;
		string 		tmp_scaffold_id("xxxxxxx");
 		contig 		tmp_contig;
		vector<contig> tmp_scaffold_info;
				boost ::regex expression( "\\s+" );
				typedef vector<string> split_vector_type;
				split_vector_type SplitVec;
		for(;;)
		{

			getline(Fin, buf);

			if(buf[0] == '>' ) 		newLine = 1;  // find Scofold
			else if(Fin.eof()) 		newLine = 3;
			else 					newLine = 2;

			if (newLine == 1)
			{
				if(!tmp_scaffold_info.empty())
				{
					cout << tmp_scaffold_id << endl;
					connectcontig(tmp_scaffold_info, tmp_scaffold_id, out_scaffold, input_contig);
					tmp_scaffold_info.clear();
				}
				cout << buf << endl;
				tmp_pos =  buf.find_first_of(" ");
			//	if(tmp_pos != string::npos)
			//		tmp_scaffold_id = buf.substr(1, tmp_pos);
			}
			if (newLine == 2)
			{
				string test = boost::regex_replace( buf, expression,"\t");
				SplitVec.clear();
				boost::iter_split( SplitVec, test, boost::first_finder("\t") );
				cout << SplitVec.at(0)<< "\t" << SplitVec.at(1)<< "\t" << SplitVec.at(2)<< "\t"<< endl;
				tmp_contig.contig_id = SplitVec.at(0);
				tmp_contig.contig_pos = boost::lexical_cast<int>(SplitVec.at(1));
				tmp_contig.recomplement = SplitVec.at(2);
				tmp_scaffold_info.push_back( tmp_contig );
			}
			if(newLine == 3)
			{
				// process the last scaffold;
				if(!tmp_scaffold_info.empty())
				{
					connectcontig(tmp_scaffold_info, tmp_scaffold_id, out_scaffold, input_contig);
					tmp_scaffold_info.clear();
				}
				break;
			}

		} // end for
	}// end try
	catch (char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}


