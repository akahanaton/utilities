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
#ifndef MATRIX_H
#define MATRIX_H

#include "general.h"

using namespace std;


enum BACKPOS {ABOVE_LEFT, LEFT,ABOVE};

/**
	@author wenming <wenming@genomics.org.cn>
*/
class NUC44
{
private:
	//void Init();
	int (*Matrix)[15];
	map <string, int> Nuc2Index;
	int GapPenalty;
	int GapExtPenalty;
public:
	NUC44();
	~NUC44();
	int getNuc44Score(int entry1, int entry2) const;
	int getNucIndex(string input_nuc);
	void setGapPenalty( int input_gappen);
	void setGapExtPenalty( int input_gapextpen);
	int getGapPenalty()   { return GapPenalty; };
	int getGapExtPenalty() { return GapExtPenalty; };
};

#endif
