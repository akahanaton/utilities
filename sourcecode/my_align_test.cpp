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
#include "alignment.h"
#include "DNA.h"
#include "blastparser.h"

using namespace std;


int main(int argc, char *argv[])
{
if(argc == 1)
{
  cout << "filter blast result" << endl;
  cout << "Use Threshold: Identity defalt 95 changable" << endl;
  cout << "               coverage defalt 90 unchangable" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f *.fasta:  contigs sequences"<<endl;
  cout << "\t\t-r *.fasta:  reference sequences"<<endl;
  cout << "\t\t-b tab string: blast result file, m 8 format"<<endl;
  cout << "\t\t-i int format: Identity, default 95.0"<<endl;
}
else
{
      int opt;
      float identity = 95.0;
      char *in_seq, *blast_file, *ref_seq, *out_seq;
  while((opt = getopt(argc, argv, "f:r:b:i:o:")) != -1 ) 
  {
      switch(opt) 
      {
      case 'f':
	  in_seq = optarg;
	  break;
      case 'r':
	  ref_seq = optarg;
	  break;
      case 'o':
	  out_seq = optarg;
	  break;
      case 'b':
	  blast_file = optarg;
	  break;
      case 'i':
         identity = boost::lexical_cast<float>(optarg);
          break;
      case '?':
	  printf("unknown option: %c\n", optopt);
	  return EXIT_SUCCESS;
	  break;
      }
   }
   
  if(!in_seq || !blast_file)
  {
    cout << "File read in error, please check the option. "<< endl;
    return EXIT_SUCCESS;
  }
  else{
  
  ///blast result filter
//   blastParser bp;
//   bp.SetReads(in_seq);
//   bp.readBlast8Result(blast_file,identity);
//   bp.RmShortAlignReads();
//   bp.getResult();
  
  /// map contig to ref
  blastParser bp;
  bp.SetRefSeq(ref_seq);
  bp.SetReads(in_seq);
  bp.readBlast8Result(blast_file,identity);
  if(out_seq)
    bp.mapContigsToRef(out_seq);
  }

//}
}

  return EXIT_SUCCESS;
}
