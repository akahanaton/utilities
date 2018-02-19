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
  cout << "Filter the CDS information whitin the genbank file" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-g genbank format: input file name" << endl;
  cout << "\t\t-f fasta format: input file name, contain the sequence that the genbank do not cotained." << endl;
  cout << "\t\t-o fasta format: output file name" << endl;
}
else
{
      clock_t start,finish;
      double totaltime;
      start=clock();
      
      int opt;
      int longer_len = 0, shorter_len = 0, equal_len = 0;
      char *in_genbank = 0, *out_fasta=0, *in_fasta= 0;
      
  while((opt = getopt(argc, argv, "g:f:o:")) != -1 ) 
  {
      switch(opt) 
      {
      case 'g':
	  in_genbank = optarg;
	  break;
      case 'f':
	  in_fasta = optarg;
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
  
  DNA mydna,mycds;
  if (in_genbank)
   mydna.readGenBank(in_genbank); 
  else 
    {
    cout << "Sequence file did not specificed, exit" << endl;
    return EXIT_SUCCESS;
    }
    
  if (in_fasta)
  { 
  	// erase the Sequence in the genbank file, and read it in the fasta files 
  	mydna.chain.clear();
  	mydna.readFasta(in_fasta);
  }
    int dna_num = mydna.getChainNum();
    int tmp_len, total_len = 0;
    cout << "Read " << dna_num << " sequence in."<< endl;
    cout << endl;
    
    fastaDNA tmp_seq;
    int cds_num = mydna.getCDSNum();
    
    for (int i = 0 ; i < cds_num; ++i)
      {
	    tmp_seq = mydna.getCDSSeq(0,i);
	    mycds.addChain(tmp_seq);
      }
    if(out_fasta)
    {
      mycds.writeFasta(out_fasta);
      cout << cds_num << " CDS are writern to " << out_fasta << endl;
   }
    else 
    {     
      for (int i = 0; i < mycds.getChainNum(); ++i)
	  cout << mycds.chain.at(i).description << endl;
    }
    
   finish=clock();
   totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
   
   cout<<"\nProgram takes "<<totaltime<<" seconds"<<endl;
   
}

  return EXIT_SUCCESS;
}
