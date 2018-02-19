#include "alignment.h"


Alignment::Alignment(string input_str1,string input_str2)
:RawSeq1(input_str1),RawSeq2(input_str2)
{
	Row = input_str1.length() + 1;	
	Col = input_str2.length() + 1;
	BackInfo.SeqEnd1 = false;	
	BackInfo.SeqEnd2 = false;	
	BuiltMatrix(Row, Col);
}

void Alignment::BuiltMatrix(int input_row, int input_col)
{
	int d = nuc44.getGapPenalty();

	AlignMatrix.resize(input_row);
	for (int i = 0; i < input_row; i++)
		AlignMatrix[i].resize(input_col);

	// count the F(i,0) 
	for (int i = 0; i < input_row; i++)
		AlignMatrix[i][0] = -i * d;
	// count the F(0,j) 
	for (int j = 0; j < input_col; j++)
		AlignMatrix[0][j] = -j * d;

	// count the F(i,j)
	int Nuc44Score, row_index, col_index;
	for (int i = 1; i < input_row; i++)
	{
		row_index = nuc44.getNucIndex( boost::lexical_cast<string>(RawSeq1.at(i-1)) );
		for (int j = 1; j < input_col; j++)
		{
		col_index = nuc44.getNucIndex( boost::lexical_cast<string>(RawSeq2.at(j-1)) );
	 	Nuc44Score = nuc44.getNuc44Score(row_index, col_index);
		AlignMatrix[i][j] = MAX3(AlignMatrix[i-1][j-1] + Nuc44Score, AlignMatrix[i-1][j]-d, AlignMatrix[i][j-1]-d);
		}
	}
}

void Alignment::Traceback( int input_row, int input_col)
{	
//	cout << "input_row: " << input_row <<"\tinput_col: " << input_col << "\t" ;

  int row_index = nuc44.getNucIndex( boost::lexical_cast<string>(RawSeq1.at(input_row - 1)) );
  int col_index = nuc44.getNucIndex( boost::lexical_cast<string>(RawSeq2.at(input_col - 1)) );
	
	int Nuc44Score = nuc44.getNuc44Score(row_index, col_index);
	if (AlignMatrix[input_row][input_col] == AlignMatrix[input_row - 1][input_col - 1] + Nuc44Score )
		{
			if( !BackInfo.SeqEnd2 && !BackInfo.SeqEnd1 )
				{BackInfo.BackRow = input_row - 1; 
				BackInfo.BackCol = input_col - 1;
				BackInfo.BackPos = ABOVE_LEFT;}
	//		cout << AlignMatrix[input_row][input_col] << endl;
		}
	if (AlignMatrix[input_row][input_col] == AlignMatrix[input_row - 1][input_col] - nuc44.getGapPenalty())
		{
			if( !BackInfo.SeqEnd1 )
				{BackInfo.BackRow = input_row - 1;
				BackInfo.BackCol = input_col;
				BackInfo.BackPos = ABOVE;}
	//		cout << AlignMatrix[input_row][input_col] << endl;
		} 
	if (AlignMatrix[input_row][input_col] == AlignMatrix[input_row][input_col - 1] - nuc44.getGapPenalty() )
		{
			if( !BackInfo.SeqEnd2 )
				{BackInfo.BackCol = input_col - 1;
				 BackInfo.BackRow = input_row;
				 BackInfo.BackPos = LEFT;}
	//		cout << AlignMatrix[input_row][input_col] << endl;
		} 
	if (BackInfo.BackRow == 0)
		{BackInfo.SeqEnd1 = true;}
	if (BackInfo.BackCol == 0)
		{BackInfo.SeqEnd2 = true;}
}

void Alignment::GlobalAlign( )
{
	string AlignSeq1, AlignSeq2;
	BackInfo.BackRow = Row - 1; BackInfo.BackCol = Col - 1;
	AlignSeq1.append( boost::lexical_cast<string>(RawSeq1[ BackInfo.BackRow ]) );
	AlignSeq2.append( boost::lexical_cast<string>(RawSeq2[ BackInfo.BackCol ]) );
	while( 1 )
	{
		if (BackInfo.SeqEnd1 && BackInfo.SeqEnd2)
			break;
		else 
			Traceback( BackInfo.BackRow, BackInfo.BackCol);
	    switch (BackInfo.BackPos)
		{
	      case ABOVE_LEFT :
		AlignSeq1.append( boost::lexical_cast<string>(RawSeq1[ BackInfo.BackRow ]) );
		AlignSeq2.append( boost::lexical_cast<string>(RawSeq2[ BackInfo.BackCol ]) );
		    break;
	      case LEFT :
		AlignSeq1.append( boost::lexical_cast<string>(RawSeq1[ BackInfo.BackCol ]) );
		AlignSeq2.append("-");
		    break;
	      case ABOVE : 
		AlignSeq1.append("-");
		AlignSeq2.append( boost::lexical_cast<string>(RawSeq2[ BackInfo.BackRow ]) );
		    break;
		}
	//	cout << "in while  "<< BackInfo.SeqEnd1<<" " << BackInfo.SeqEnd2<<endl;
	}
	cout << "AligSeq1:\t";
	int len = AlignSeq1.length();
	for ( int i = len - 1; i >= 0; i--)
		cout <<  AlignSeq1.at(i);
	cout << endl;
	cout << "AligSeq2:\t";

	len = AlignSeq2.length();
	for ( int i = len - 1; i>= 0; i--)
		cout <<  AlignSeq2.at(i);
	cout << endl;

}

void Alignment::getAlignMatrix()
{
	cout << "  \t  \t";
	for(int i = 0; i < Col - 1 ; i++ )
	 	cout << RawSeq2.at(i) << "\t";
	cout<<endl;
	for (int i = 0; i < Row; i++)
	{
		
		if(i  ==  0 )
			cout << " \t";
		else
			cout << RawSeq1.at(i-1)<<"\t";
		for (int j = 0; j < Col; j++)
			{
			cout << AlignMatrix[i][j]<<"\t";	
		}
		cout<<endl;
	}
}



/*!
    \fn Alignment::LocalAlign
 */
void Alignment::LocalAlign()
{
    /// @todo implement me
}
