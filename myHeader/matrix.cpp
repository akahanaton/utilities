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
#include "nuc44.h"

NUC44::NUC44()
{
   Matrix = new int[15][15];
  int tmpMatrix[15][15] =
   { 5, -4, -4, -4,  1, -4, -4,  1, -4,  1, -4, -1, -1, -1, -2,
    -4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2,
    -4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2,
    -4, -4, -4,  5, -4,  1,  1, -4, -4,  1, -1, -1, -1, -4, -2,
     1, -4,  1, -4, -1, -4, -2, -2, -2, -2, -3, -1, -3, -1, -1,
    -4,  1, -4,  1, -4, -1, -2, -2, -2, -2, -1, -3, -1, -3, -1,
    -4, -4,  1,  1, -2, -2, -1, -4, -2, -2, -1, -1, -3, -3, -1,
     1,  1, -4, -4, -2, -2, -4, -1, -2, -2, -3, -3, -1, -1, -1,
    -4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1,
     1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1,
    -4, -1, -1, -1, -3, -1, -1, -3, -1, -3, -1, -2, -2, -2, -1,
    -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1,
    -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1,
    -1, -1, -1, -4, -1, -3, -3, -1, -1, -3, -2, -2, -2, -1, -1,
    -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };

	for (int i = 0; i < 15; i++)
	  for (int j = 0; j < 15; j++)
		Matrix[i][j] = tmpMatrix[i][j];
	Nuc2Index["A"] = 0;
	Nuc2Index["C"] = 1;
	Nuc2Index["G"] = 2;
	Nuc2Index["T"] = 3;
	Nuc2Index["R"] = 4;
	Nuc2Index["Y"] = 5;
	Nuc2Index["K"] = 6;
	Nuc2Index["M"] = 7;
	Nuc2Index["S"] = 8;
	Nuc2Index["W"] = 9;
	Nuc2Index["B"] = 10;
	Nuc2Index["D"] = 11;
	Nuc2Index["H"] = 12;
	Nuc2Index["V"] = 13;
	Nuc2Index["N"] = 14;

	GapPenalty = 8;
	GapExtPenalty = 2;
}

NUC44::~NUC44()
{
	delete [] Matrix;
}

void NUC44::setGapPenalty( int input_gappen)
{
	GapPenalty = input_gappen;
}

void NUC44::setGapExtPenalty( int input_gapextpen)
{
	GapExtPenalty = input_gapextpen;
}

int NUC44::getNuc44Score(int entry1,int entry2) const
{
 try{
	if (entry1 < 0 || entry1 > 14)
		throw entry1;
 	if( entry2 < 0 || entry2 > 14)
		throw entry2;
	return Matrix[entry1][entry2];
 }//end try
 catch (int Msg) { cerr << endl << "Entry for NCU44 must between 0-14,Exception: " << Msg<< endl; }
}


int NUC44::getNucIndex(string input_nuc)
{
	return Nuc2Index[input_nuc];
}


