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

#include <iostream>
#include <cstdlib>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include "/share/raid5/wenming/myProgram/myHeader/DNA.h"

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
      start=clock();

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
    int tmp_len, total_len = 0, N50_Total = 0, N50Size = 0;
    bool N50_set = false;


  /// sequence length statistics
     map<int,int> reads_len_counter;
     map<int,int>::const_iterator map_it;
     vector<int> reads_len;

    for(int i = 0; i < dna_num; i++)
	{
 	 tmp_len =  mydna.chain.at(i).length();
	 cout << mydna.chain.at(i).id << "\t" << tmp_len << endl;
	}

	}
  return EXIT_SUCCESS;
}
