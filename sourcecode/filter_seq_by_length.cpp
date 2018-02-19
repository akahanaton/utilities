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
#include "DNA.h"

using namespace std;

int main(int argc, char *argv[])
{
if(argc == 1)
{
  cout << "Filter the sequence by the specific length" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f fasta format: input file name, contain the sequence that the genbank do not cotained." << endl;
  cout << "\t\t-o fasta format: output file name" << endl;
  cout << "\t\t-l int format: longer than." << endl;
  cout << "\t\t-s int format: shorter than." << endl;
  cout << "\t\t-e int format: equal" << endl;
}
else
{
      clock_t start,finish;
      double totaltime;
      start=clock();

      int opt;
      int longer_len = 0, shorter_len = 0, equal_len = 0;
      char *in_genbank = 0, *out_fasta=0, *in_fasta= 0;

  while((opt = getopt(argc, argv, "f:o:l:s:e:")) != -1 )
  {
      switch(opt)
      {
      case 'f':
	  in_fasta = optarg;
	  break;
      case 'o':
	  out_fasta = optarg;
	  break;
      case 'l':
	  longer_len = boost::lexical_cast<int>(optarg);
	  break;
      case 's':
	  shorter_len = boost::lexical_cast<int>(optarg);
	  break;
      case 'e':
	  equal_len = boost::lexical_cast<int>(optarg);
	  break;
      case '?':
	  printf("unknown option: %c\n", optopt);
	  return EXIT_SUCCESS;
	  break;
      }
   }

  DNA mydna, my_out_dna;

  if (in_fasta)
   mydna.readFasta(in_fasta);

    int dna_num = mydna.getChainNum();
    int tmp_len, total_len = 0;
    cout << "Read " << dna_num << " sequence in."<< endl;
    cout << endl;

/// Filter the sequence whit specific length
    if (equal_len > 0)
    {
	  for (int i = 0; i < dna_num; ++i)
	    if (mydna.chain.at(i).length() == equal_len)
			my_out_dna.addChain(mydna.chain.at(i));
    }
    else if (equal_len == 0)
    {
      if (longer_len > 0 )
	  {
		  if (shorter_len == 0)
		  {
			for (int i = 0; i < dna_num; ++i)
				if (mydna.chain.at(i).length() >= longer_len)
					my_out_dna.addChain(mydna.chain.at(i));
		  }
		  else if (shorter_len > 0 )
		  {
			 for (int i = 0; i < dna_num; ++i)
			  {
				if ( mydna.chain.at(i).length() >= shorter_len && mydna.chain.at(i).length() <= longer_len )
				  my_out_dna.addChain(mydna.chain.at(i));
			  }
		 }
	  }
      else if (longer_len == 0 )
      {
        if (shorter_len > 0)
		  {
			for (int i = 0; i < dna_num; ++i)
			if (mydna.chain.at(i).length() <= shorter_len )
				my_out_dna.addChain(mydna.chain.at(i));
		  }
      }
    }
    if(out_fasta)
	{
	   cout<< my_out_dna.getChainNum() <<" sequence writing into file "<< out_fasta <<endl;
	   my_out_dna.writeFasta(out_fasta);
	}
    else
    {
      for (int i = 0; i < my_out_dna.getChainNum(); ++i)
		  cout << my_out_dna.chain.at(i).description << endl;
    }
}
  return EXIT_SUCCESS;
}
