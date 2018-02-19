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
#include "../myHeader/blastparser.h"

using namespace std;

void usage();

struct PairedReads
{
	string Seq_1;
	string Seq_2;
	int    Pos_1;
	int    Pos_2;
};


typedef map<string, PairedReads> PairedInfo;

//void MatePair( PairedInfo& out_pe, const SoapParser& in_sp );

int binary_pos(vector<single_gap>& gap_regions, const int& regions_num, int& key_pos );
void DealSoapResual( char* in_file, char *out_file, vector<single_gap>& gaps);

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

      int opt;
      int align_len_cutoff = 0;
      float identity = 90.0;
      char *in_seq = 0, *blast_file = 0, *soap_file = 0, *out_file = 0;
	  while((opt = getopt(argc, argv, "f:b:s:o:")) != -1 )
	  {
		 switch(opt)
		  {
			  case 'f':
				in_seq = optarg;
				break;
			  case 'b':
				blast_file = optarg;
				break;
			  case 's':
				soap_file = optarg;
				break;
			  case 'o':
				out_file = optarg;
				break;
			  case '?':
				 printf("unknown option: %c\n", optopt);
				 return EXIT_SUCCESS;
				 break;
		  }
	   }

	  ///blast result filter
	   BlastParser bp;
	   if(in_seq)
		   bp.SetReads(in_seq);
	   if(blast_file)
	   {
		   bp.ReadBlast8Result(blast_file,identity);
		   bp.FilterBlast8Result(true, true, align_len_cutoff, 0);
	   }
		vector<single_gap> gr;
		bp.GetGapRegion(gr);

		bp.BlastResults.clear();
	    bp.Reads.clear();

		cout << gr.size()<< endl;
		finish=clock();
		totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
		cout<<"\nProgram takes "<<totaltime<<" seconds"<<endl;

		DealSoapResual( soap_file, out_file, gr);

		finish=clock();
		totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
		cout<<"\nProgram takes "<<totaltime<<" seconds"<<endl;
	}// end else

  return EXIT_SUCCESS;
}


void ConnectContigAccording2PE()
{
  // pre: 1 Reads must have been seted;
  // 	  2 blast results must have been readed in, and hvae filtered useing the function FilterBlast8Result(false,150);
  // 	  3 soap paired hits results must have been readed in;
  // post: According to the alignment results between contigs and reference (blast results)
  //       every two contigs which have overlap or the gaps between them is shorter than 100bp is connected as one;
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

void DealSoapResual( char* in_file, char *out_file, vector<single_gap>& gaps)
{
	try
	{
		ifstream Fin(in_file);
		if(!Fin) throw strcat( "Cannot open input result file", in_file);

		ofstream Fout(out_file);

		string buf, tmp_str;
		int gaps_num = gaps.size();
		int tmp_pos, in_which_gap;
		int tab_1;
		for(;;)
		{
			getline(Fin, buf);
			if(Fin.eof())
				break;
			else
			{
				tab_1 = buf.find("->");
				if( tab_1 != string::npos )
						tmp_str = buf.substr(tab_1 -20, 16);
				else
						tmp_str = buf.substr( buf.size()-15, 13);

				tab_1 = tmp_str.find("\t");
				tmp_pos = boost::lexical_cast<int>(tmp_str.substr(tab_1+1));
				in_which_gap = binary_pos( gaps, gaps_num, tmp_pos);
				if(in_which_gap > 0)
					Fout << buf << endl;
			 }
		}//end for
	//	sort(SoapResults.begin(), SoapResults.end(),less_soap_location());
		Fin.close();
		Fout.close();
	}//end try
	catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl;}
}

int binary_pos(vector<single_gap>& gap_regions, const int& regions_num, int& key_pos )
{
	int start =0, end = regions_num -1;
	while(start < end)
	{
		int k = (start + end)/2;
		if( key_pos >= gap_regions.at(k).pos1 && key_pos <= gap_regions.at(k).pos2)
			return k;
		if( key_pos <= gap_regions.at(k).pos1)
			end = k -1;
		else
			start = k + 1;
	}//end while
	return -1;
}

