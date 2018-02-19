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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>
#include <iostream>
#include <cstdlib>
#include "/share/raid1/wenming/myHeader/DNA.h"
using namespace std;

int main(int argc, char *argv[])
{
	
        DNA dna;
   //     dna.readGenBank(argv[1], "CDS");
	if (argc != 3 )
	{
		cout << "Count the oligo repeats in the fasta sequence;"<< endl;
		cout << "Input 1: the fasta files contains the sequence."<< endl;
		cout << "Input 2: the repeat's length, int, like 13 or 15"<< endl;
		cout << "use > to write the result into."<< endl;
  		return EXIT_SUCCESS;
	}
        dna.readFasta(argv[1]);
	int dna_num = dna.getChainNum();
	int oligo_length = boost::lexical_cast<int>(argv[2]);
	for (int i = 0; i < dna_num; ++i)
	{
		dna.oligoCount(i,oligo_length);
	}
	dna.getOligoNum();
  return EXIT_SUCCESS;
}
