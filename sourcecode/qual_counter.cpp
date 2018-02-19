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

#include <numeric>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include "../myHeader/DNA.h"

using namespace std;

int main(int argc, char *argv[])
{
if(argc == 1)
{
  cout << "Count the average qual of all sequences" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f fasta format: input file name, contain the sequence that the genbank do not cotained." << endl;
  cout << "\t\t-q fasta format: quality file name" << endl;
}
else
{
      clock_t start,finish;
      double totaltime;
      start=clock();
      
      int opt;
      char *in_qual = 0, *in_fasta = 0;
      
  while((opt = getopt(argc, argv, "f:q:")) != -1 ) 
  {
      switch(opt) 
      {
      case 'q':
	  	in_qual = optarg;
	  break;
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

   if(in_qual)
     mydna.readQual(in_qual);
	 long int totalqual = 0;
	 long int totalbase = 0; 
	 map< string, vector<int> >::iterator map_it(mydna.qual.begin());
	 	for( ; map_it != mydna.qual.end(); ++ map_it)
	 		{
	 			unsigned int tmq_qual = 0;
	 			vector<int>::iterator vec_it(map_it->second.begin());
	 		//			for ( ; vec_it != map_it->second.end(); ++ vec_it);
	 		//				totalqual += *vec_it;
	 			totalqual += accumulate( map_it->second.begin(), map_it->second.end(), 0);
	 			totalbase += map_it->second.size();
	 		}
	 		
	 cout << "average qual "<< totalqual / totalbase << endl;
	 
 finish=clock();
   totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
   
   cout<<"\nProgram takes "<<totaltime<<" seconds"<<endl;
   
}


   
  return EXIT_SUCCESS;
}
