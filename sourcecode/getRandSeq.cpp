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
#include "myHeader/DNA.h"

using namespace std;

int main(int argc, char *argv[])
{
if(argc == 1)
{
  cout << "getRandSeq, generat random sequence from a template sequence" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f *.fasta: fasta file name, containing the template sequence"<<endl;
  cout << "\t\t-o *.fasta: out file name"<<endl;
  cout << "\t\t-l int format: random sequence length"<<endl;
  cout << "\t\t-c int format:  coverage to cover the template sequence"<<endl;
}
else
{
      int opt, rand_len, coverage;
      char *in_file, *out_file;
      
  while((opt = getopt(argc, argv, "f:o:l:c:")) != -1 ) {
      switch(opt) {
      case 'f':
	  in_file = optarg;
	  break;
      case 'o':
	  out_file = optarg;
	  break;
      case 'l':
          rand_len = boost::lexical_cast<int>(optarg);
          break;
      case 'c':
          coverage = boost::lexical_cast<int>(optarg);
      case '?':
	  printf("unknown option: %c\n", optopt);
	  break;
      }
  }
  
  for(; optind < argc; optind++)
      printf("argument: %s\n", argv[optind]);  
      
  DNA mydna, RandDna;
  if( in_file != 0 && out_file != 0 )
  {
  mydna.readFasta( in_file );
  cout << mydna.chain.size()<<endl;
  cout << "........" << endl;
  cout << "generating sequence: sequence length: " << rand_len << endl;
  cout << "                     coverage: " << coverage << endl;
  if (coverage >= 1 &&  rand_len > 0)
    {
    getRandSeq( rand_len, coverage, mydna, &RandDna);
    RandDna.writeFasta( out_file );
    }
  }
}

  return EXIT_SUCCESS;
}
