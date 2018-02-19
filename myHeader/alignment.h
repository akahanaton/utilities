#ifndef ALIGN_H
#define ALIGN_H

#include "general.h"
#include "matrix.h"

class Alignment
{
public:
	Alignment(string input_str1, string input_str2);
	Alignment() {};
	//void align();
	void BuiltMatrix( int input_row, int input_col);
	void Traceback( int input_row, int input_col);
	void GlobalAlign();
	void getAlignMatrix();
        void LocalAlign();
	
private:
	string RawSeq1,RawSeq2;
	int Row,Col;
	vector< vector<int> > AlignMatrix;	
	int BestAlignScore;
	NUC44 nuc44;
	struct BACKINFO
		{	
			int BackRow;
			int BackCol;
			BACKPOS BackPos;
			bool SeqEnd1, SeqEnd2;
		} BackInfo;
};

# endif