/***************************************************************************
 *   Copyright (C) 2007 by bob   *
 *   bob@datarig   *
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
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "/genome/myProgram/BioC++.h"
using namespace std;

int main(int argc, char *argv[])
{
        DNA dna;
	if (argc != 3)
	{
    cout << "input 1: The Sequencs files, fasta format." << endl;
    cout << "input 2: The repeat sequencs files , like : " << endl;
    cout << "\tAATAATGATAATAAT" << endl;
    cout << "\tATATATATATATATA" << endl;
    cout << "\tGGCGGCGGCGGCG" << endl;
    cout << "\tAATAATGATAATAAT" << endl;
  	return EXIT_SUCCESS;
	}
        dna.readFasta(argv[1]);	
			
	
				
	ifstream Fin_rep(argv[2]);
	if(!Fin_rep) 
		cout << "Can't open repeats file: " << argv[2] << endl;
	string repeat;
	vector<string> repeats;
	for(;;)
	{
	getline(Fin_rep, repeat);
	
	if(Fin_rep.eof()) 
		break;
	else
	 {
		repeats.push_back(repeat);
	 }
	}// end for i
	int repeats_num = repeats.size();
	// regex repeat_tag(repeat);
 	 vector<fastaDNA>::iterator tmpIt;
	for( tmpIt = dna.chain.begin(); tmpIt != dna.chain.end(); tmpIt++ )
	  {
		for( int i = 0; i < repeats_num; i++)
		 {
		 if(dna.getOligo2Id(repeats.at(i),tmpIt))
			break;
		}
	  }
  return EXIT_SUCCESS;
}
