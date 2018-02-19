#include "blastparser.h"

////////////////// member of class BlastParser ///////////////////
BlastParser::BlastParser()
{
	reads_set = false;
	ref_set = false;
}

BlastParser::~BlastParser()
{
}

void BlastParser::AddQuery( const BLAST8RESULT& input_result)
{
  BlastResults.push_back(input_result);
}
void BlastParser::ReadBlastResult(char* resultfile, const int& results_num)
{
	try
	{
	ifstream Fin(resultfile);
	if(!Fin) throw strcat("Cannot open input result file ", resultfile);
	string buf;
	int tmp_pos = 0;
	int Query_len = 0;
	string str_num, str_QEnd, str_SEnd;
	bool QStart_set = false, SStart_set = false;
	BLAST8RESULT *tmp_result_ptr;
	BLAST8RESULT tmp_result;
	tmp_result_ptr = &tmp_result;
	if (results_num > 0)
		BlastResults.reserve(results_num);
		for(;;)
		{
			getline(Fin, buf);
			if(Fin.eof())
			{
				if(QStart_set)
				{
					tmp_result_ptr->QEnd = boost::lexical_cast<int>(str_QEnd);
					tmp_result_ptr->SEnd = boost::lexical_cast<int>(str_SEnd);
					if(tmp_result_ptr->QEnd < tmp_result_ptr->QStart)
						myswap( tmp_result_ptr->QEnd, tmp_result_ptr->QStart );
					if(tmp_result_ptr->SEnd < tmp_result_ptr->SStart)
						myswap( tmp_result_ptr->SEnd, tmp_result_ptr->SStart );
					QStart_set = false;
					SStart_set = false;
					AddQuery(tmp_result);
				}
				break;
			}
			else if(buf.find("Query=") == 0 )  //query name
			{
				if(QStart_set)
				{
					tmp_result_ptr->QEnd = boost::lexical_cast<int>(str_QEnd);
					tmp_result_ptr->SEnd = boost::lexical_cast<int>(str_SEnd);
					if(tmp_result_ptr->QEnd < tmp_result_ptr->QStart)
						myswap( tmp_result_ptr->QEnd, tmp_result_ptr->QStart );
					if(tmp_result_ptr->SEnd < tmp_result_ptr->SStart)
						myswap( tmp_result_ptr->SEnd, tmp_result_ptr->SStart );
					QStart_set = false;
					SStart_set = false;
					AddQuery(tmp_result);
				}
				tmp_result.Init();
				tmp_pos = buf.find(' ', 7);
				if (tmp_pos != string::npos)
					tmp_result_ptr->QueryID = buf.substr(7, tmp_pos - 7 );
				else
					tmp_result_ptr->QueryID = buf.substr(7);
				getline(Fin, buf);
				tmp_pos = buf.find("letters)");
				while( tmp_pos == string::npos  )
				{
					getline(Fin, buf);
					tmp_pos = buf.find("letters)");
				}
				tmp_pos = buf.find("(");
				str_num = buf.substr(tmp_pos + 1, buf.find(" ", tmp_pos) - tmp_pos);
                //--------------------------------------------------
                // printf("unknown option: %s\n", "here4");
                //--------------------------------------------------
                //--------------------------------------------------
                // cout<<"aaa"<<str_num<<"bbb"<<endl;
                //--------------------------------------------------
				tmp_result_ptr->QLen =  boost::lexical_cast<int>(str_num.substr(0, remove(str_num.begin(), str_num.end(),' ') - str_num.begin()) );
                //--------------------------------------------------
                // printf("unknown option: %s\n", "here5");
                //--------------------------------------------------
			}
			else if(buf[0] == '>') // database info
			{
				if(QStart_set)
				{
					tmp_result_ptr->QEnd = boost::lexical_cast<int>(str_QEnd);
					tmp_result_ptr->SEnd = boost::lexical_cast<int>(str_SEnd);
					if(tmp_result_ptr->QEnd < tmp_result_ptr->QStart)
						myswap( tmp_result_ptr->QEnd, tmp_result_ptr->QStart );
					if(tmp_result_ptr->SEnd < tmp_result_ptr->SStart)
						myswap( tmp_result_ptr->SEnd, tmp_result_ptr->SStart );
					QStart_set = false;
					SStart_set = false;
					AddQuery(tmp_result);
				}
				tmp_pos = buf.find(' ');
				if (tmp_pos != string::npos)
					tmp_result_ptr->SubjectID= buf.substr(1, tmp_pos - 1 );
				else
					tmp_result_ptr->SubjectID= buf.substr(1);
				getline(Fin, buf);
				tmp_pos = buf.find("Length =");
				while( tmp_pos == string::npos  )
				{
					getline(Fin, buf);
					tmp_pos = buf.find("Length =");
				}
				tmp_pos = buf.find("=");
				str_num = buf.substr(tmp_pos + 2);
				tmp_result_ptr->SLen = boost::lexical_cast<int>(str_num);
                //--------------------------------------------------
                // printf("unknown option: %s\n", "here6");
                //--------------------------------------------------
			}
			else if(buf.find("Score =") != string::npos)
			{
				if(QStart_set)
				{
					tmp_result_ptr->QEnd = boost::lexical_cast<int>(str_QEnd);
					tmp_result_ptr->SEnd = boost::lexical_cast<int>(str_SEnd);
					if(tmp_result_ptr->QEnd < tmp_result_ptr->QStart)
						myswap( tmp_result_ptr->QEnd, tmp_result_ptr->QStart );
					if(tmp_result_ptr->SEnd < tmp_result_ptr->SStart)
						myswap( tmp_result_ptr->SEnd, tmp_result_ptr->SStart );
					QStart_set = false;
					SStart_set = false;
					AddQuery(tmp_result);
				}
				//tmp_result.Init();
				tmp_pos = buf.find("="); // Score
				tmp_pos = buf.find_first_of("123456789", tmp_pos);
				str_num = buf.substr(tmp_pos, buf.find(" ", tmp_pos) - tmp_pos);
				tmp_result_ptr->BitScore = boost::lexical_cast<float>(str_num);
                //--------------------------------------------------
                // printf("unknown option: %s\n", "here7");
                //--------------------------------------------------

				getline(Fin, buf); //Identities and alignment length
				tmp_pos = buf.find("=");
				str_num = buf.substr(tmp_pos + 2, buf.find("/", tmp_pos) - tmp_pos - 2);
				Query_len = boost::lexical_cast<int>(str_num);
                //--------------------------------------------------
                // printf("unknown option: %s\n", "here8");
                //--------------------------------------------------
				tmp_pos = buf.find("/", tmp_pos);
				str_num = buf.substr(tmp_pos, buf.find(" ", tmp_pos) - tmp_pos);
				str_num = buf.substr(tmp_pos + 1, buf.find(" ", tmp_pos) - tmp_pos -1);
				tmp_result_ptr->AlignmentLength = boost::lexical_cast<int>(str_num);
                //--------------------------------------------------
                // printf("unknown option: %s\n", "here9");
                //--------------------------------------------------
				tmp_pos = buf.find("(", tmp_pos);
				str_num = buf.substr(tmp_pos + 1, buf.find("%", tmp_pos) - tmp_pos -1);
				tmp_result_ptr->Identity= boost::lexical_cast<int>(str_num);
                //--------------------------------------------------
                // printf("unknown option: %s\n", "here10");
                //--------------------------------------------------
				tmp_pos = buf.find("=", tmp_pos);
				if( tmp_pos != string::npos)// Gaps exist
				{
					str_num = buf.substr(tmp_pos + 2, buf.find("/", tmp_pos) - tmp_pos - 2);
					tmp_result_ptr->Gaps = boost::lexical_cast<int>(str_num);
					tmp_result_ptr->Mismatches = tmp_result_ptr->AlignmentLength - Query_len - tmp_result_ptr->Gaps;
				}
				else // Gaps doesn't exist
				{
					tmp_result_ptr->Gaps = 0;
					tmp_result_ptr->Mismatches = tmp_result_ptr->AlignmentLength - Query_len;
				}

				getline(Fin, buf); // RComplement
				if(buf.find("Minus") != string::npos)
					tmp_result_ptr->RComplement = true;
				else
					tmp_result_ptr->RComplement = false;
			}
			else if(buf.find("Query:") != string::npos)
			{
				if(!QStart_set)
				{
					tmp_pos = buf.find(" ");
					str_num = buf.substr(tmp_pos + 1, buf.find(" ", tmp_pos + 1) - tmp_pos -1);
					tmp_result_ptr->QStart = boost::lexical_cast<int>(str_num);
					QStart_set = true;
				}
                tmp_pos = buf.find_last_of(" ");
                str_QEnd = buf.substr(tmp_pos + 1);
			}
			else if(buf.find("Sbjct:") != string::npos)
			{
				if(!SStart_set)
				{
					tmp_pos = buf.find(" ");
					str_num = buf.substr(tmp_pos + 1, buf.find(" ", tmp_pos + 1) - tmp_pos -1);
					tmp_result_ptr->SStart = boost::lexical_cast<int>(str_num);
					SStart_set = true;
				}
                tmp_pos = buf.find_last_of(" ");
                str_SEnd= buf.substr(tmp_pos + 1);
			}
			else
				continue;
		}//end for(;;)
	}// end try
	catch(char* pMsg) {cerr << endl << "Exception:" << pMsg << endl;}
}

void BlastParser::MapReadsToRef( char *out_filename, float input_identity)
{
  /// pre:  blast results, Reads and Reference sequence have been read in
  /// post: according to the blast results, reads are mapped to the reference
  ///       and a map sequence writte into the file out_filename


  /// firstly:  check the reads and reference sequence are not empty
   if( RefSeq.empty() || Reads.empty() )
  {
    cout << "Map need enter the reads sequence, exit now!" << endl;
    return ;
  }

   /// Secondly: processing the blast result
   /// step 1: remove the short aligned reads , less than 90%;
   /// step 2: sort the blast result by SStart (Subject Start loc);
   /// step 3: remove the repeats,according to the SStart;

  unsigned int QueryNum = BlastResults.size();
  cout << "before UniqByQueryID " << QueryNum << endl;

  UniqBy_QueryID_SStart();
  QueryNum = BlastResults.size();
  cout << "after UniqByQueryID " << QueryNum << endl;


  cout << "\t\tbegin Remove Short Align Results"<< endl;
  QueryNum = BlastResults.size();
  //for (int i = 0; i < QueryNum; i++)
      //cout << BlastResults.at(i).QueryID<< endl;
      cout << QueryNum << endl;

  /// step 1
  // RmShortAlignReads();

  cout << "\t\tafter Remove Short Align Results"<< endl;
  QueryNum = BlastResults.size();
//   for (int i = 0; i < QueryNum; i++)
//       cout << BlastResults.at(i).QueryID<< endl;
//       cout << QueryNum << endl;

  /// step 2
   sort(BlastResults.begin(), BlastResults.end(),less_by_SbjctID_SStart());


   cout << "\t\tafter sort  ";
//     for (int i = 0; i < QueryNum; i++)
//       cout << BlastResults.at(i).QueryID << "\t\t" <<BlastResults.at(i).SStart << "\t\t"<<BlastResults.at(i).SEnd << endl;
    cout << BlastResults.size()<< endl;


  /// step 3
  RmRepeatReads();

  cout << "\t\tafter remove repeat "<< endl;

  QueryNum = BlastResults.size();
//     for (int i = 0; i < QueryNum; i++)
//       cout << BlastResults.at(i).QueryID << "\t\t" <<BlastResults.at(i).SStart << "\t\t"<<BlastResults.at(i).SEnd << endl;
    cout <<"useful blast result: " <<BlastResults.size()<< endl;


  /// thirdly: begin map
  /// need change for multiseq

  string myMapSeq, tmp_id,tmp_seq1, tmp_seq2;
  vector<BLAST8RESULT>:: iterator cur_it,blast_it;
  int i = 1;

  blast_it = BlastResults.begin(); /// the first reference
  if ( blast_it->SStart > 0 )
    myMapSeq.append( blast_it->SStart -1 ,'X');
  if( blast_it->RComplement )
  {
    tmp_seq1 = Reads[blast_it->QueryID].substr(blast_it->QStart - 1, blast_it->AlignmentLength );
    myMapSeq.append( tmp_seq1 );
  }
  if(blast_it->RComplement)
    {
      tmp_seq1 = Reads[blast_it->QueryID].substr(blast_it->QStart - 1 , blast_it->AlignmentLength );
      /// reverse compliment
      reverse_copy(tmp_seq1.begin(), tmp_seq1.end(), tmp_seq2.begin());
      transform(tmp_seq2.begin(), tmp_seq2.end(), tmp_seq1.begin(),compliment_nucleotide);
    } /// end if RComplement
      else
    {
        tmp_seq1 = Reads[blast_it->QueryID].substr(blast_it->QStart - 1,blast_it->AlignmentLength);
    }
    cout << "map  " << i++ << "\t" <<blast_it->QueryID << "\t" << tmp_seq1.length() << endl;
    myMapSeq.append(tmp_seq1);


  cur_it = blast_it + 1;
  for ( ; cur_it != BlastResults.end()-1; cur_it++)
     {
     if ( cur_it->SStart > blast_it->SEnd )
     {
        myMapSeq.append( cur_it->SStart - blast_it->SEnd ,'X');
     }
     else
     {
        cur_it->SStart = blast_it->SEnd + 1;
        cur_it->QStart = cur_it->QStart + blast_it->SEnd - cur_it->SStart + 1;
     }
     blast_it++;

     tmp_id =cur_it->QueryID;

      if( cur_it->RComplement )
      {
      tmp_seq1 = Reads[tmp_id].substr(cur_it->QStart - 1, cur_it->QEnd - cur_it->QStart + 1 );
      /// reverse compliment
      reverse_copy(tmp_seq1.begin(), tmp_seq1.end(), tmp_seq2.begin());
      transform(tmp_seq2.begin(), tmp_seq2.end(), tmp_seq1.begin(),compliment_nucleotide);
      }
      else
      {
        tmp_seq1 = Reads[tmp_id].substr(cur_it->QStart - 1, cur_it->QEnd - cur_it->QStart + 1 );
      }


      //cout << "map" << i++ << "\t" << cur_it->QueryID << "\t" << tmp_seq1.length() << endl;
      //cout << tmp_seq1 <<endl;
      myMapSeq.append( tmp_seq1 );
    } /// end for

  if(BlastResults.back().SEnd > blast_it->SStart )
  {
      myMapSeq.append( BlastResults.back().SStart - blast_it->SEnd ,'X');
  }
  else
  {
        BlastResults.back().SStart = blast_it->SEnd + 1;
        BlastResults.back().QStart = BlastResults.back().QStart +
				     blast_it->SEnd - BlastResults.back().SStart + 1;
  }
  tmp_id = BlastResults.back().QueryID;

  if( BlastResults.back().RComplement )
  {
  tmp_seq1 = Reads[tmp_id].substr(BlastResults.back().QStart - 1,
                                BlastResults.back().QEnd - BlastResults.back().QStart + 1 );
  /// reverse compliment
  reverse_copy(tmp_seq1.begin(), tmp_seq1.end(), tmp_seq2.begin());
  transform(tmp_seq2.begin(), tmp_seq2.end(), tmp_seq1.begin(),compliment_nucleotide);
  }
  else
  {
    tmp_seq1 = Reads[tmp_id].substr(BlastResults.back().QStart - 1,
                                  BlastResults.back().QEnd - BlastResults.back().QStart + 1 );
  }

  cout << cur_it->QueryID << "\t" << BlastResults.back().QueryID<<endl;
  cout << "last map " << i++ << "\t" << tmp_id << "\t" << tmp_seq1.length() << endl;
  myMapSeq.append( tmp_seq1 );

  if ( BlastResults.back().SEnd < RefSeq.begin()->second.length() ) ///float  back() return the last reference
   myMapSeq.append( RefSeq.begin()->second.length() - BlastResults.back().SEnd ,'X');

    cout << "Ref Seq length:" <<  RefSeq.begin()->second.length() << endl;
    cout << "my mapped seq len " <<myMapSeq.length() << endl;
    tmp_id = "my ref seq";
    DNA myRefDna;
    myRefDna.addChain( tmp_id, tmp_id, myMapSeq );
    myRefDna.writeFasta( out_filename );
}

void BlastParser::SetReads( char *infile )
{
    DNA tmpDna;
    tmpDna.readFasta( infile );
    unsigned int readnum = tmpDna.getChainNum();
    for (unsigned int i = 0; i < readnum; i++ )
      Reads.insert(map<string, string>::value_type(tmpDna.chain.at(i).id, tmpDna.chain.at(i).seq));
	reads_set = true;
}

void BlastParser::SetRefSeq( char *infile )
{
    DNA tmpDna;
    tmpDna.readFasta( infile );
    unsigned int readnum = tmpDna.getChainNum();
    for (unsigned int i = 0; i < readnum; i++ )
      RefSeq.insert(map<string, string>::value_type(tmpDna.chain.at(i).id, tmpDna.chain.at(i).seq));
	ref_set = true;
}

void BlastParser::UniqBy_QueryID( )
{
  if (BlastResults.size() > 0 )
  {
      vector<BLAST8RESULT>::iterator it;
       it = unique( BlastResults.begin(), BlastResults.end(), equal_by_QueryID());
       BlastResults.resize( it - BlastResults.begin());
  }
}

void BlastParser::UniqBy_SubjectID( )
{
  if (BlastResults.size() > 0 )
  {
      vector<BLAST8RESULT>::iterator it;
       it = unique( BlastResults.begin(), BlastResults.end(), equal_by_SubjectID());
       BlastResults.resize( it - BlastResults.begin());
  }
}

void BlastParser::UniqBy_QueryID_SStart( )
{
  if (BlastResults.size() > 0 )
  {
      vector<BLAST8RESULT>::iterator it;
       it = unique( BlastResults.begin(), BlastResults.end(), equal_by_QueryID_SStart());
       BlastResults.resize( it - BlastResults.begin());
  }
}

void BlastParser::RmRepeatReads(int inter )
{
    vector<BLAST8RESULT>:: iterator it;
    it = unique( BlastResults.begin(), BlastResults.end(), rm_repeat_by_SStart_read());
    BlastResults.resize( it - BlastResults.begin());
}

void BlastParser::GetResult(bool head_info)
{
	if(head_info)
		cout<<"QueryID\t"<<"SubjectID\t"<<"Identity\t"<<"AlignmentLength\t"
			<<"Mismatches\t"<<"Gaps\t"<<"QLen\t"<<"QStart\t"<<"QEnd\t"
			<<"SLen\t"<<"SStart\t"<<"SEnd"<<"BitScore\t"<<"RComplement\t"<<endl;
      cout<< BlastResults.at(0).QueryID << "\t"
	      << BlastResults.at(0).SubjectID << "\t"
	      << BlastResults.at(0).Identity << "\t"
	      << BlastResults.at(0).AlignmentLength << "\t"
		  << BlastResults.at(0).Mismatches << "\t"
	      << BlastResults.at(0).Gaps << "\t"
		  << BlastResults.at(0).QLen << "\t"
	      << BlastResults.at(0).QStart << "\t"
	      << BlastResults.at(0).QEnd << "\t"
		  << BlastResults.at(0).SLen << "\t"
	      << BlastResults.at(0).SStart << "\t"
	      << BlastResults.at(0).SEnd << "\t"
	     // << BlastResults.at(0).SStart - 1 << "\t"
		  //useless
	    //<< BlastResults.at(0).Expect << "\t"
	      << BlastResults.at(0).BitScore << "\t"
	      << BlastResults.at(0).RComplement << "\t"
	      << endl;
int QueryNum = BlastResults.size();
    for (int i = 1 ; i < QueryNum; i++)
    {
    cout<< BlastResults.at(i).QueryID << "\t"
	    << BlastResults.at(i).SubjectID << "\t"
	    << BlastResults.at(i).Identity << "\t"
	    << BlastResults.at(i).AlignmentLength << "\t"
	    << BlastResults.at(i).Mismatches << "\t"
	    << BlastResults.at(i).Gaps << "\t"
		<< BlastResults.at(i).QLen << "\t"
	    << BlastResults.at(i).QStart << "\t"
	    << BlastResults.at(i).QEnd << "\t"
		<< BlastResults.at(i).SLen << "\t"
	    << BlastResults.at(i).SStart << "\t"
	    << BlastResults.at(i).SEnd << "\t"
		//<< BlastResults.at(i).SStart - BlastResults.at(i-1).SEnd << "\t"
		//useless
	//	<< BlastResults.at(i).Expect << "\t"
	    << BlastResults.at(i).BitScore << "\t"
	    << BlastResults.at(i).RComplement << "\t"
	    << endl;
    }
}


void BlastParser::MapContigsToRef(char *out_seq, int align_len_threshold, bool gap_region, char *align_file, char *out_align_depth)
{
  /// pre:  blast results, Reads and Reference sequence have been read in
  /// post: according to the blast results, reads are mapped to the reference
  ///       and a map sequence writte into the file out_seq

	DNA tmp_dna;
	if (align_file)
			tmp_dna.Read454AlignDepth( align_file );

  /// firstly:  check the reads and reference sequence are not empty
   if( RefSeq.empty() || Reads.empty() )
  {
    cout << "Map need enter the reads sequence, exit now!" << endl;
    return ;
  }
   /// Secondly: processing the blast result
   /// step 1: sort the result by QueryID and QStart (Query Start loc);
   /// step 2: remove the short aligned reads,
   ///         according to blast result,for each contig,only longest part
   ///         aligned to the reference are keeped back,
   ///         and no overlap between each part;
   /// step 3: sort the blast result by SStart (Subject Start loc);
   /// step 4: rm the repeat Alignment in the reference according to the SStart;
   /// step 5: adjust the QStart and QEnd postion,
   ///				 and the SStart and SEnd postion is also adjusted according to that;


  /// step 1
  sort(BlastResults.begin(), BlastResults.end(),less_by_QueryID_QStart());
 // GetResult( true );

  /// step 2
  RmMussyAlignInContig(align_len_threshold);
 // GetResult( true );

  /// step 3
   sort(BlastResults.begin(), BlastResults.end(),less_by_SbjctID_SStart());
   //GetResult( true );

  /// step 4
  RmRepeatAlignBy_QStart();
  //GetResult(true);

  /// stem 5
  RmQueryInMiddle();
  GetResult(true);

  /// step 6
  vector<BLAST8RESULT>:: iterator cur_it, pre_it;

   /// thirdly: begin map
  /// need change for multi sbjct
  /// initialize the a raw sequence stand by 'X'
  /// and the depth at each position is 0;
  string myMapSeq(RefSeq.begin()->second.length(), 'X');

  vector<int> tmp_depth (RefSeq.begin()->second.length(),0);

  string tmp_id,tmp_des,tmp_seq1,tmp_seq2;
	int i = 1;

  for ( cur_it = BlastResults.begin(); cur_it != BlastResults.end(); ++cur_it)
  {
    if( cur_it->RComplement )
     {
		  tmp_seq1 = Reads[cur_it->SubjectID].substr(cur_it->SStart - 1 , cur_it->SEnd - cur_it->SStart + 1 );

		  /// reverse compliment
			reverse(tmp_seq1.begin(), tmp_seq1.end());
			tmp_seq2 = tmp_seq1;
			transform(tmp_seq2.begin(), tmp_seq2.end(), tmp_seq1.begin(),compliment_nucleotide);
     }
    else
     {
		  tmp_seq1 = Reads[cur_it->SubjectID].substr(cur_it->SStart - 1 , cur_it->SEnd - cur_it->SStart + 1);
	 }

      /// string& replace( size_type index, size_type num, const string& str );
     myMapSeq.replace(cur_it->QStart -1, cur_it->QEnd - cur_it->QStart + 1,tmp_seq1);

     if (align_file )
      copy( tmp_dna.depth[cur_it->QueryID].begin() + cur_it->QStart - 1,
			tmp_dna.depth[cur_it->QueryID].begin() + cur_it->QEnd - 1,
			tmp_depth.begin() + cur_it->SStart - 1 );

	 /// get the uncovered region
	 if(gap_region)
      RefSeq.begin()->second.replace(RefSeq.begin()->second.begin() + cur_it->QStart -1,
								     RefSeq.begin()->second.begin() + cur_it->QEnd - 1,
									 cur_it->QEnd - cur_it->QStart + 1,
									 'X');
   }

    tmp_id = "map seq";
    tmp_des = tmp_id.append("\tlength=");
    tmp_des.append(boost::lexical_cast<string>(myMapSeq.length()));
    tmp_dna.addChain( tmp_id, tmp_des, myMapSeq );
    if (align_file)
	{
		tmp_dna.depth.clear();
		tmp_dna.depth.insert(map<SEQID,DEPTH>::value_type(tmp_id, tmp_depth));
		tmp_dna.Write454AlignDepth(out_align_depth);
	}
    tmp_dna.writeFasta(out_seq);

	if(gap_region)
	{
		char gap_file[50];
		strcpy(gap_file, out_seq);
		strcat(gap_file, ".GapRegion");
		tmp_dna.clear();
		tmp_dna.addChain( tmp_id,tmp_des, RefSeq.begin()->second );
		tmp_dna.writeFasta( gap_file );
	}
    return ;
}

void BlastParser::GetUncoverSeq(char* out_seq)
{
  vector<BLAST8RESULT>:: iterator cur_it, pre_it;

  /// need change for multiseq
  /// initialize the a raw sequence stand by 'X'
  string myMapSeq(RefSeq.begin()->second.length(), 'X');
  string tmp_id,tmp_des,tmp_seq1,tmp_seq2;
  int i = 1;

  for ( cur_it = BlastResults.begin(); cur_it != BlastResults.end(); ++cur_it)
  {
    if( cur_it->RComplement )
     {
		  tmp_seq1 = Reads[cur_it->QueryID].substr(cur_it->QStart - 1 , cur_it->QEnd - cur_it->QStart + 1 );

		  /// reverse compliment
			reverse(tmp_seq1.begin(), tmp_seq1.end());
			tmp_seq2 = tmp_seq1;
			transform(tmp_seq2.begin(), tmp_seq2.end(), tmp_seq1.begin(),compliment_nucleotide);
     }
    else
     {
			tmp_seq1 = Reads[cur_it->QueryID].substr(cur_it->QStart - 1,cur_it->QEnd - cur_it->QStart + 1);
	 }

      /// string& replace( size_type index, size_type num, const string& str );
     myMapSeq.replace(cur_it->SStart -1, cur_it->SEnd - cur_it->SStart + 1,tmp_seq1);

      RefSeq.begin()->second.replace(RefSeq.begin()->second.begin() + cur_it->SStart -1,
								     RefSeq.begin()->second.begin() + cur_it->SEnd - 1,
									 cur_it->SEnd - cur_it->SStart + 1,
									 'X');
   }

	char gap_file[50];
	strcpy(gap_file, out_seq);
	strcat(gap_file, ".GapRegion");
	DNA tmp_dna;
	tmp_dna.clear();
	tmp_dna.addChain( tmp_id,tmp_des, RefSeq.begin()->second );
	tmp_dna.writeFasta( gap_file );
    return ;
}

void BlastParser::
GetEstCoverInfo(const int& query_len_threshold,const int& sbjct_len_threshold,const int& align_len_threshold,const int& identity_threshold,bool sort_tag)
{
  /// step 1
  if(query_len_threshold > 0)
		RmShortQuery(query_len_threshold);
  if(sbjct_len_threshold > 0)
		RmShortSbjct(sbjct_len_threshold);
  if(align_len_threshold > 0)
		RmShortAlign(align_len_threshold);
  if(identity_threshold > 0)
		RmLowIdentity(identity_threshold);
  if(sort_tag) {
	//--------------------------------------------------
	//   sort(BlastResults.begin(), BlastResults.end(),less_by_QueryID());
	//--------------------------------------------------
	  sort(BlastResults.begin(), BlastResults.end(),less_by_QueryID_SbjctID());
  }

  vector<BLAST8RESULT>:: iterator cur_it = BlastResults.begin();
  set<string> sbjct_id_set;
  for(; cur_it != BlastResults.end(); ++cur_it)
	  sbjct_id_set.insert( cur_it->SubjectID );
  cout << "Hit number " << sbjct_id_set.size() << endl;

   /// thirdly: begin map
  cur_it = BlastResults.begin();

  int cur_covered_base = 0, max_cover_base = 0;
  string max_cover_sbjct;
  map<string, int> ref_info;
  map<string, int> cov_info;
  map<string, string> qurey_sbjct;
  string cur_query_seq(cur_it->QLen, 'N');
  cur_query_seq.replace(cur_it->QStart - 1, cur_it->QEnd-cur_it->QStart + 1, cur_it->QEnd-cur_it->QStart + 1, 'X');
  cur_covered_base = count(cur_query_seq.begin(), cur_query_seq.end(), 'X');
  max_cover_base = cur_covered_base;

//--------------------------------------------------
//   string& replace ( size_t pos1, size_t n1, size_t n2, char c );
//   string& replace ( iterator i1, iterator i2, size_t n2, char c );
//   The section is replaced by a repetition of character c, n2 times.
//--------------------------------------------------

  vector<BLAST8RESULT>:: iterator pre_it = cur_it;
  ++cur_it;
  for ( ; cur_it != BlastResults.end(); ++cur_it)
  {
	  if(cur_it->QueryID == pre_it->QueryID) {
		  if(cur_it->SubjectID == pre_it->SubjectID)
		  {
			cur_query_seq.replace(cur_it->QStart - 1, cur_it->QEnd-cur_it->QStart + 1, cur_it->QEnd-cur_it->QStart + 1, 'X');
		  }
		  else //same query, but aligned to anthor subjcet seq
		  {
			  cur_covered_base = count(cur_query_seq.begin(), cur_query_seq.end(), 'X');
			  if(cur_covered_base > max_cover_base)
			  {
				  max_cover_base = cur_covered_base;
				  max_cover_sbjct = pre_it->SubjectID;
			  }
			  cur_query_seq.clear();
			  cur_query_seq = string(cur_it->QLen, 'N');
			  cur_query_seq.replace(cur_it->QStart - 1, cur_it->QEnd-cur_it->QStart + 1, cur_it->QEnd-cur_it->QStart + 1, 'X');
		  }
	  } else { // a new query
		  cur_covered_base = count(cur_query_seq.begin(), cur_query_seq.end(), 'X');
		  if(cur_covered_base > max_cover_base)
		  {
			  max_cover_base = cur_covered_base;
			  max_cover_sbjct = pre_it->SubjectID;
		  }
		  cov_info.insert(map<string, int>::value_type( pre_it->QueryID, max_cover_base));
		  ref_info.insert(map<string, int>::value_type( pre_it->QueryID, pre_it->QLen));
		  qurey_sbjct.insert(map<string, string>::value_type( pre_it->QueryID, max_cover_sbjct));
		  cur_query_seq.clear();
		  cur_query_seq = string(cur_it->QLen, 'N');
		  cur_query_seq.replace(cur_it->QStart -1, cur_it->QEnd-cur_it->QStart + 1, cur_it->QEnd-cur_it->QStart + 1, 'X');
		  cur_covered_base = count(cur_query_seq.begin(), cur_query_seq.end(), 'X');
		  max_cover_base = cur_covered_base;
		  max_cover_sbjct = cur_it->SubjectID;
	  }
	  ++pre_it;
  }// end for


  --cur_it;
  pre_it = cur_it -1;
  if(cur_it->QueryID == pre_it->QueryID)
  {
	  if(cur_it->SubjectID == pre_it->SubjectID)
	  {
		  cur_covered_base = count(cur_query_seq.begin(), cur_query_seq.end(), 'X');
		  if(cur_covered_base > max_cover_base)
		  {
			  max_cover_base = cur_covered_base;
			  max_cover_sbjct = pre_it->SubjectID;
		  }
	  }
	  cov_info.insert(map<string, int>::value_type( cur_it->QueryID, max_cover_base));
	  ref_info.insert(map<string, int>::value_type( cur_it->QueryID, cur_it->QLen));
	  qurey_sbjct.insert(map<string, string>::value_type( cur_it->QueryID, max_cover_sbjct));
  }
  else
  {
	  cov_info.insert(map<string, int>::value_type( cur_it->QueryID, max_cover_base));
	  ref_info.insert(map<string, int>::value_type( cur_it->QueryID, cur_it->QLen));
	  qurey_sbjct.insert(map<string, string>::value_type( cur_it->QueryID, max_cover_sbjct));
  }

		map<string, int>::iterator map_it(cov_info.begin() );
		int low_50 = 0, betw_50_60 = 0, betw_60_70 = 0, betw_70_80 = 0, betw_80_90 = 0, betw_90_100 = 0;
		int base_low_50 = 0, base_betw_50_60 = 0, base_betw_60_70 = 0, base_betw_70_80 = 0, base_betw_80_90 = 0, base_betw_90_100 = 0;
		int hit_seq_num =0;
		int total_base = 0, total_covered_base=0, query_len = 0;
		float max_covered = 0, min_covered = 100, cur_covered = 0;
		cout << "Sequence ID\tContig ID\tCoverage %"<<endl;

		for( ; map_it != cov_info.end(); ++map_it)
		{
			total_covered_base += map_it->second;
			query_len = ref_info[map_it->first];
			total_base += query_len;
			cur_covered = 100.00 * map_it->second / query_len;

			if(cur_covered <= 50)
			{
				++low_50;
				base_low_50 += map_it->second;
			}
			else if( cur_covered > 50 && cur_covered <= 60)
			{
				++betw_50_60;
				base_betw_50_60 += map_it->second;
			}
			else if( cur_covered > 60 && cur_covered <= 70)
			{
				++betw_60_70;
				base_betw_60_70 += map_it->second;
			}
			else if( cur_covered > 70 && cur_covered <= 80)
			{
				++betw_70_80;
				base_betw_70_80 += map_it->second;
			}
			else if( cur_covered > 80 && cur_covered <= 90)
			{
				++betw_80_90;
				base_betw_80_90 += map_it->second;
			}
			else
			{
				++betw_90_100;
				base_betw_90_100 += map_it->second;
			}

			if(cur_covered > max_covered)
				max_covered = cur_covered;
			if(cur_covered < min_covered)
				min_covered = cur_covered;
			cout << map_it->first << "\t" << qurey_sbjct[map_it->first] << "\t" << query_len << "\t" << cur_covered << endl;
		}

		cout << endl <<"Total info:" << endl;
		int total_seq = ref_info.size();
		cout << "Total hited sequence " << total_seq << "\tTotal base " << total_base << endl;
		cout << "Total covered bases " << total_covered_base << " , Coverage ratio "
			 << 100.00 * total_covered_base / total_base << endl;
		cout <<"\tCoverage between 90---100%\t"<< betw_90_100 << "\tcovered base\t"<< 100.00*base_betw_90_100 / total_base <<endl;
		cout <<"\tCoverage between 80---90%\t" << betw_80_90 << "\tcovered base\t"<< 100.00*base_betw_80_90 / total_base<<endl;
		cout <<"\tCoverage between 70---80%\t" << betw_70_80 << "\tcovered base\t"<< 100.00*base_betw_70_80/total_base<<endl;
		cout <<"\tCoverage between 60---70%\t" << betw_60_70 << "\tcovered base\t"<< 100.00*base_betw_60_70/total_base<<endl;
		cout <<"\tCoverage between 50---60%\t" << betw_50_60 << "\tcovered base\t"<< 100.00*base_betw_50_60/total_base<<endl;
		cout <<"\tCoverage lower than 50%\t\t" << low_50 << "\tcovered base\t"<<100.00* base_low_50/total_base <<endl;
		cout <<"\tCoverage higher than 50%\t" << total_seq - low_50 <<endl;
		cout << "\tMax coverage\t\t" << max_covered << endl;
		cout << "\tMin coverage\t\t" << min_covered << endl;
}

void BlastParser::
GetContigCoverInfo(const int& query_len_threshold,const int& sbjct_len_threshold,const int& align_len_threshold,const int& identity_threshold,bool sort_tag)
{
  /// step 1
  if(query_len_threshold > 0)
		RmShortQuery(query_len_threshold);
  if(sbjct_len_threshold > 0)
		RmShortSbjct(sbjct_len_threshold);
  if(align_len_threshold > 0)
		RmShortAlign(align_len_threshold);
  if(identity_threshold > 0)
		RmLowIdentity(identity_threshold);
  if(sort_tag)
	  sort(BlastResults.begin(), BlastResults.end(),less_by_QueryID());

  vector<BLAST8RESULT>:: iterator cur_it = BlastResults.begin();
  set<string> sbjct_id_set;
  for(; cur_it != BlastResults.end(); ++cur_it)
	  sbjct_id_set.insert( cur_it->SubjectID );
  cout << "Hit number " << sbjct_id_set.size() << endl;

   /// thirdly: begin map
  cur_it = BlastResults.begin();

  int cur_covered_base = 0;
  map<string, int> ref_info;
  map<string, int> cov_info;
  string cur_query_seq(cur_it->QLen, 'N');
  cur_query_seq.replace(cur_query_seq.begin() + cur_it->QStart -1, cur_query_seq.begin() + cur_it->QEnd -1,
						cur_it->QEnd - cur_it->QStart + 1, 'X');
  cur_covered_base = count(cur_query_seq.begin(), cur_query_seq.end(), 'X');

  vector<BLAST8RESULT>:: iterator pre_it;
  ++cur_it;
  for ( ; cur_it != BlastResults.end(); ++cur_it)
  {
	  pre_it = cur_it - 1;
	  if(cur_it->QueryID == pre_it->QueryID)
	  {
			cur_query_seq.replace(cur_it->QStart - 1, cur_it->QEnd-cur_it->QStart + 1, cur_it->QEnd-cur_it->QStart + 1, 'X');
	  }
	  else // a new query
	  {
		  cur_covered_base = count(cur_query_seq.begin(), cur_query_seq.end(), 'X');
		  cov_info.insert(map<string, int>::value_type( pre_it->QueryID, cur_covered_base));
		  ref_info.insert(map<string, int>::value_type( pre_it->QueryID, pre_it->QLen));
		  cur_query_seq = string(cur_it->QLen, 'N');
		  cur_query_seq.replace(cur_it->QStart -1, cur_it->QEnd-cur_it->QStart + 1, cur_it->QEnd-cur_it->QStart + 1, 'X');
	  }
  }// end for

  pre_it = cur_it - 1;
  cur_covered_base = count(cur_query_seq.begin(), cur_query_seq.end(), 'X');
  cov_info.insert(map<string, int>::value_type( pre_it->QueryID, cur_covered_base));
  ref_info.insert(map<string, int>::value_type( pre_it->QueryID, pre_it->QLen));

		map<string, int>::iterator map_it(cov_info.begin() );
		int low_50 = 0, betw_50_60 = 0, betw_60_70 = 0, betw_70_80 = 0, betw_80_90 = 0, betw_90_100 = 0;
		int hit_seq_num =0;
		int total_base = 0, total_covered_base=0, query_len = 0;
		float max_covered = 0, min_covered = 100, cur_covered = 0;
		cout << "Sequence ID\tContig ID\tCoverage %"<<endl;

		for( ; map_it != cov_info.end(); ++map_it)
		{
			total_covered_base += map_it->second;
			query_len = ref_info[map_it->first];
			total_base += query_len;
			cur_covered = 100.00 * map_it->second / query_len;

			if(cur_covered <= 50)
				++low_50;
			else if( cur_covered > 50 && cur_covered <= 60)
				++betw_50_60;
			else if( cur_covered > 60 && cur_covered <= 70)
				++betw_60_70;
			else if( cur_covered > 70 && cur_covered <= 80)
				++betw_70_80;
			else if( cur_covered > 80 && cur_covered <= 90)
				++betw_80_90;
			else
				++betw_90_100;

			if(cur_covered > max_covered)
				max_covered = cur_covered;
			if(cur_covered < min_covered)
				min_covered = cur_covered;
			cout << map_it->first << "\t" <<  query_len << "\t" << cur_covered << endl;
		}

		cout << endl <<"Total info:" << endl;
		int total_seq = ref_info.size();
		cout << "Total hited sequence " << total_seq << endl;
		cout << "Total covered bases " << total_covered_base << " , Coverage ratio "
			 << 100.00 * total_covered_base / total_base << endl;
		cout <<"\tCoverage between 90---100%\t"<< betw_90_100 <<endl;
		cout <<"\tCoverage between 80---90%\t" << betw_80_90 <<endl;
		cout <<"\tCoverage between 70---80%\t" << betw_70_80 <<endl;
		cout <<"\tCoverage between 60---70%\t" << betw_60_70 <<endl;
		cout <<"\tCoverage between 50---60%\t" << betw_50_60 <<endl;
		cout <<"\tCoverage lower than 50%\t\t" << low_50 <<endl;
		cout <<"\tCoverage higher than 50%\t" << total_seq - low_50 <<endl;
		cout << "\tMax coverage\t\t" << max_covered << endl;
		cout << "\tMin coverage\t\t" << min_covered << endl;
}


void BlastParser::RmLowIdentity(int identity_threshold)
{
  BlastResults.erase(
      remove_if( BlastResults.begin(), BlastResults.end(), bind2nd(check_align_identity(),identity_threshold) ),
      BlastResults.end()
      );
}

void BlastParser::RmShortAlign(int len_threshold)
{
  BlastResults.erase(
      remove_if( BlastResults.begin(), BlastResults.end(), bind2nd(check_align_len(),len_threshold) ),
      BlastResults.end()
      );
}

void BlastParser::RmShortQuery(int len_threshold)
{
  BlastResults.erase(
      remove_if( BlastResults.begin(), BlastResults.end(), bind2nd(check_qurey_len(),len_threshold) ),
      BlastResults.end()
      );
}

void BlastParser::RmShortSbjct(int len_threshold)
{
  BlastResults.erase(
      remove_if( BlastResults.begin(), BlastResults.end(), bind2nd(check_sbjct_len(),len_threshold) ),
      BlastResults.end()
      );
}

void BlastParser::RmLowQueryAlignRatio(int ratio_threshold)
{
	BlastResults.erase(
		remove_if( BlastResults.begin(), BlastResults.end(), bind2nd(check_query_align_ratio(), ratio_threshold) ),
		BlastResults.end()
		);
}

void BlastParser::RmLowSbjctAlignRatio(int ratio_threshold)
{
	BlastResults.erase(
		remove_if( BlastResults.begin(), BlastResults.end(), bind2nd(check_sbjct_align_ratio(), ratio_threshold) ),
		BlastResults.end()
		);
}

void BlastParser::RmMussyAlignInContig(int len_threshold)
{
  /// remove the mussy result
	BlastResults.erase(
		remove_if( BlastResults.begin(), BlastResults.end(), bind2nd(check_align_len_ratio(), len_threshold) ),
		BlastResults.end()
		);
	vector<BLAST8RESULT>::iterator it;
	    it = unique( BlastResults.begin(), BlastResults.end(), equal_by_QueryID_QStart());
	    BlastResults.resize( it - BlastResults.begin());
	}


void BlastParser::RmRepeatAlignBy_SStart()
{
  /// remove the mussy result
  vector<BLAST8RESULT>::iterator it;
       it = unique( BlastResults.begin(), BlastResults.end(), rm_repeat_by_SStart_contig());
       BlastResults.resize( it - BlastResults.begin());
}

void BlastParser::RmRepeatAlignBy_QStart()
{
  /// remove the mussy result
  vector<BLAST8RESULT>::iterator it;
       it = unique( BlastResults.begin(), BlastResults.end(), rm_repeat_by_QStart_contig());
       BlastResults.resize( it - BlastResults.begin());
}

void BlastParser::RmQueryInMiddle()
{
	vector<BLAST8RESULT>::iterator Begin = BlastResults.begin() + 1;
	vector<BLAST8RESULT>::iterator End = BlastResults.end() - 1;
	//--------------------------------------------------
	// BlastResults.erase(
	// 	remove_if( Begin, End, check_query_in_middle()),
	// 	BlastResults.end()
	// 	);
	//--------------------------------------------------
}

void BlastParser::AssemblyWrongCheck(int len_threshold, int identity_threshold)
{
	BlastParser pos_wrong_result;
	BlastParser dir_wrong_result;

	if(len_threshold > 0)
		RmShortAlign(len_threshold);
	if(identity_threshold > 0)
		RmLowIdentity(identity_threshold);

	RmRepeatAlignBy_QStart();
	sort(BlastResults.begin(), BlastResults.end(),less_by_SbjctID_SStart());
	vector<BLAST8RESULT>::iterator cur_it = BlastResults.begin() + 1;
	vector<BLAST8RESULT>::iterator pre_it = cur_it -1, post_it = cur_it +1;

	while( cur_it != BlastResults.end() -1)
	{
		if(cur_it->SubjectID == pre_it->SubjectID && cur_it->SubjectID == post_it->SubjectID)
		{
			if(cur_it->RComplement != pre_it->RComplement)
			{
				dir_wrong_result.AddQuery(*pre_it);
				dir_wrong_result.AddQuery(*cur_it);
			}

			if((cur_it->QEnd > pre_it->QEnd) && (cur_it->QEnd > post_it->QEnd ))
			{
				pos_wrong_result.AddQuery(*pre_it);
				pos_wrong_result.AddQuery(*cur_it);
				pos_wrong_result.AddQuery(*post_it);
			}
			else if((cur_it->QEnd < pre_it->QEnd) && (cur_it->QEnd >  post_it->QEnd ))
			{
				pos_wrong_result.AddQuery(*pre_it);
				pos_wrong_result.AddQuery(*cur_it);
				pos_wrong_result.AddQuery(*post_it);
			}
		}
		++cur_it;
		++pre_it;
		++post_it;
	}//end while
	sort(pos_wrong_result.BlastResults.begin(), pos_wrong_result.BlastResults.end(),less_by_SbjctID_SStart());
//	pos_wrong_result.RmRepeatAlignBy_QStart();
	pos_wrong_result.UniqBy_SubjectID();
	cout << "Number of postion wrong sequence:\t" <<pos_wrong_result.BlastResults.size()<<endl;

	sort(dir_wrong_result.BlastResults.begin(), dir_wrong_result.BlastResults.end(),less_by_SbjctID_SStart());
	//dir_wrong_result.RmRepeatAlignBy_QStart();
	dir_wrong_result.UniqBy_SubjectID();
	cout << "Number of direction wrong sequence:\t" <<dir_wrong_result.BlastResults.size()<<endl;

	cout << "Position wront sequence" << endl;
	for(cur_it = pos_wrong_result.BlastResults.begin(); cur_it != pos_wrong_result.BlastResults.end(); ++cur_it)
	 cout << "\t" << cur_it->QueryID << "\t" << cur_it->SubjectID << "\t" << cur_it->QStart << "\t" << cur_it->QEnd <<endl;
	cout << "Direction wront sequence" << endl;
	for(cur_it = dir_wrong_result.BlastResults.begin(); cur_it != dir_wrong_result.BlastResults.end(); ++cur_it)
	 cout << "\t" << cur_it->QueryID << "\t" << cur_it->SubjectID << "\t" << cur_it->QStart << "\t" << cur_it->QEnd <<endl;
}

void BlastParser::FilterBlastResult(bool contig_flag, bool sort_flag, int len_threshold, int ratio_threshold)
{
	// reads and contigs blast result are filtered by different tactics
	//--------------------------------------------------
	// if(ratio_threshold > 0)
	// 	RmLowAlignRatio(ratio_threshold);
	//--------------------------------------------------
	if(len_threshold > 0)
		RmShortAlign(len_threshold);
	// sort first, then call RmMussyAlignInContig()
	sort(BlastResults.begin(), BlastResults.end(),less_by_QueryID_QStart());
	if(contig_flag) // process the reads results
		  RmMussyAlignInContig(len_threshold);
	// sort by Sbjct start pos
	if(sort_flag)
		  sort(BlastResults.begin(), BlastResults.end(),less_by_SbjctID_SStart());
}

void BlastParser::GetGapRegion( vector<single_gap>& gr, int add_len )
{
	// pre: the blast resuslts must have been filtered by the function FilterBlast8Result(bool,int)
	// post: get the gap region between two contigs
	// reference:   --------------------------------------------------------------------------
	// contigs:     -----------------------                          -------------------------
	//                       |  add_len   |       real gap           |<-- add_len -->|
	//                       |                   gap region                       -->|
	single_gap tmp_sg;

	if( BlastResults.at(0).SStart > 1 )
	{
		//--------------------------------------------------
		// tmp_sg.real_gap = BlastResults.at(0).SStart - 1;
		// tmp_sg.id1 = BlastResults.at(0).QueryID;
		// tmp_sg.id2 = BlastResults.at(0).QueryID;
		// tmp_sg.pos1 = 1;
		// tmp_sg.pos2 = BlastResults.at(0).SStart -1 + add_len;
		//--------------------------------------------------
		gr.push_back(tmp_sg);
	}
	size_t results_num = BlastResults.size();
	for( size_t i = 1; i < results_num; ++i)
	{
		if( BlastResults.at(i).SubjectID == BlastResults.at(i-1).SubjectID ) // same sbjct(referecne) sequence
		{
			//--------------------------------------------------
			// tmp_sg.real_gap = BlastResults.at(i).SStart - BlastResults.at(i-1).SEnd;
			//--------------------------------------------------
			//--------------------------------------------------
			//  if( tmp_sg.real_gap < 0)
			//--------------------------------------------------
			 {
				//--------------------------------------------------
				// tmp_sg.id1 = BlastResults.at(i-1).QueryID;
				// tmp_sg.id2 = BlastResults.at(i).QueryID;
				// tmp_sg.pos1 = BlastResults.at(i).SStart - add_len;
				// tmp_sg.pos2 = BlastResults.at(i-1).SEnd + add_len;
				//--------------------------------------------------
				gr.push_back(tmp_sg);
			 }
			//--------------------------------------------------
			//  else
			//--------------------------------------------------
			 {
				//--------------------------------------------------
				// tmp_sg.id1 = BlastResults.at(i-1).QueryID;
				// tmp_sg.id2 = BlastResults.at(i).QueryID;
				// tmp_sg.pos1 = BlastResults.at(i-1).SEnd - add_len;
				// tmp_sg.pos2 = BlastResults.at(i).SStart + add_len;
				//--------------------------------------------------
				gr.push_back(tmp_sg);
			 }
		}
		else // new sbjct(reference) sequence
		{
			if( BlastResults.at(i).SStart > 1 )
			{
				//--------------------------------------------------
				// tmp_sg.real_gap = BlastResults.at(i).SStart - 1;
				// tmp_sg.id1 = BlastResults.at(i).QueryID;
				// tmp_sg.id2 = BlastResults.at(i).QueryID;
				// tmp_sg.pos1 = 1;
				// tmp_sg.pos2 = BlastResults.at(i).SStart -1 + add_len;
				//--------------------------------------------------
				gr.push_back(tmp_sg);
			}
		}
	}
}
