#include "DNA.h"
///////////////////member of class fastaDNA ////////////////////////////////

fastaDNA::fastaDNA() {
  seq = "";
  id = "";
  description = "";

}

fastaDNA::fastaDNA( const fastaDNA &other_fastadna):
	id(other_fastadna.id),
	description(other_fastadna.description),
	seq(other_fastadna.seq)
{  }

void fastaDNA::set_id(string input_id)
{
  id = input_id;
}

void fastaDNA::set_description(string input_description)
{
  description = input_description;
}

void fastaDNA::set_seq(string input_seq)
{
  seq = input_seq;
}

void fastaDNA::initial()
{
  seq = "";
  id = "";
  description = "";
}

bool not_nucleotide(char ch)
{
  string str = "ACGTatcg";
  string ts(1,ch);
  if( str.find(ts) == string::npos ) return true;
  else return false;
}

void fastaDNA::cleanSeq()
{
  string::iterator endp = remove_if(seq.begin(), seq.end(), not_nucleotide);
  seq.erase(endp, seq.end());
}

int fastaDNA::length()
{
  return seq.size();
}

void fastaDNA::reverse_compliment()
{
  string tmp_seq;
  tmp_seq.resize(seq.length());
  reverse_copy(seq.begin(), seq.end(), tmp_seq.begin());
  transform(tmp_seq.begin(), tmp_seq.end(), seq.begin(),
            compliment_nucleotide);
}

/*!
    \fn fastaDNA::longerThan(int input_len)
 */
bool fastaDNA::longerThan(int input_len)
{
    /// @todo implement me
  return seq.length() >= input_len;
}


/*!
    \fn fastaDNA::method_1()
 */
bool fastaDNA::shorterThan(int input_len)
{
    /// @todo implement me
  return seq.length() < input_len;
}


/*DNAchain DNAchain::assemble(DNAchain c2, vector<hsp> hsps) {

DNAchain DNAchain::assemble(DNAchain c1, DNAchain c2, int minmatch) {
  int i,j;
  // DNAchain query_chain = this->length() > c2.length() ? c2 : *this ;
  // DNAchain sbjct_chain = this->length() > c2.length() ? *this : c2 ;
  // string query = query_chain.seq;
  // string sbjct = sbjct_chain.seq;
  string query = c1.seq;
  string sbjct = c2.seq;

  // first, build a list of the query's substrings
  if( sbjct->find(sbjct, query, minmatch) )

  vector<int> subStrings;
  string::iterator p;
  for(i=0; i<query.size(); i++) {
      // if the substring doesn't start and end with 'N', put it in the list
      if( query[i] != 'N' && query[p+minmatch] != 'N' ) {
	  subStrings.push_back(i);
      }
  }

  // second, search every substring in the sbjct string
  for(i=0; i<subStrings.size(); i++) {
      string sbjct = sbjct_chain.seq;
      if( sbjct->find(query, subStrings[i], minmatch) )

  }
}
*/

///////////////////member of class CDS ////////////////////////////////

CDS::CDS()
{
  initial();
}

CDS::CDS(const CDS& other_cds):
	gene(other_cds.gene),
	locus_tag(other_cds.locus_tag),
	protein_id(other_cds.protein_id),
	db_xref(other_cds.db_xref),
	id(other_cds.id),
	cds_info(other_cds.cds_info),
	codon_start(other_cds.codon_start),
	recomplement(other_cds.recomplement),
	indices(other_cds.indices)
{  }

void CDS::initial()
{
  cds_info.clear();
  db_xref.clear();
  gene.clear();
  locus_tag.clear();
  protein_id.clear();
  indices.first = 0;
  indices.second =0;
  codon_start = 0;
  recomplement = false;
}

void CDS::set_indices(int input_pos1,int input_pos2)
{
  indices = make_pair ( input_pos1, input_pos2);
}

///////////////////member of class DNA ////////////////////////////////

DNA::DNA() { }

DNA::DNA(const DNA& other_dna):
	chain(other_dna.chain),
	cds(other_dna.cds),
	oligo_counter(other_dna.oligo_counter),
	qual(other_dna.qual),
	depth(other_dna.depth)
{  }

void DNA::readFasta( char *file ) // read Fasta format sequences
{
  try {
      ifstream Fin(file);
      if(!Fin) throw strcat("Cannot open input Fasta file ", file);

      bool fastaFormat = false;
      string buf;
      int newLine;
      fastaDNA tmpChain;


      for(;;) {
	  getline(Fin, buf);
	  /* The new line is marked with an integer, which corresponds to one of the following cases:
		1. New sequence line, which starts with "> ".
		2. Sequence line, which starts with any charactor except for ">"
		3. End of file, detected by function eof().
	  */

         if(buf[0] == '>' ) newLine=1;
	  else if(Fin.eof()) newLine=3;
	  else newLine=2;

          if( newLine == 1 )
		  { // this is a start of a new chain.
			  if( !fastaFormat )
			  {
				// if it is the first chain, just build a new chain object
				tmpChain.initial();
				fastaFormat = true;
			  }
			  else
			  {
				// otherwise, need to store the old chain first, then build a new chain
				//tmpChain.cleanSeq();
				addChain(tmpChain);
				tmpChain.initial();
			  }
			  tmpChain.set_description( buf.substr(1,buf.size()-1) );
			  int pos = buf.find_first_of(" ");
			  tmpChain.set_id( buf.substr(1,pos-1) );
		  }

	  if( newLine == 2 && fastaFormat ) { // sequence line

	      tmpChain.seq.append(buf);
	  }

	  if( newLine == 3 ) // end of file
	  {
	      if( !fastaFormat ) throw strcat(file, " is not in Fasta format.");
	      else
		  {
	          // otherwise, need to store the old chain
	          // tmpChain.cleanSeq();
	          addChain(tmpChain);
	          break;
		  }
	  }
      }//endl for
	  Fin.close();
  }// end try

  catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}

void readFastaMap(char *file, map_DNA& out_seq_map) // read Fasta format sequence
{
	ifstream Fin(file);
	if(!Fin)
	{
		cout << "Can not open file,"<<file <<" in function readFastaMap" << endl;
		return;
	}
	string buf, tmp_id, tmp_seq;
	int pos;
	bool fasta_format = false;
	for(;;)
	{
		getline(Fin, buf);
		if(Fin.eof())
		{
			if(!fasta_format)
				cout << file << "is not in Fasta format." << endl;
			else
			{
				out_seq_map.insert(make_pair( tmp_id, tmp_seq ) );
				break;
			}
		}
		else if(buf.empty())
			continue;
		else if(buf[0] == '>')
		{
			if(!fasta_format)
				fasta_format = true;
			else
			{
				// add last seq and reset the tmp_id and tmp_seq
				out_seq_map.insert(make_pair( tmp_id, tmp_seq ) );
				tmp_id.clear();
				tmp_seq.clear();
			}
			pos = buf.find_first_of(" \t");
			if(pos != string::npos)
				tmp_id = buf.substr(1, pos -1 );
			else
				tmp_id = buf.substr(1);
		}
		else if(fasta_format)
		{
			tmp_seq.append(buf);
		}
	} // end for
	Fin.close();
}

void readFastqMap(char *file, map_DNA& out_seq_map) // read Fastq format sequence
{
	ifstream Fin(file);
	if(!Fin)
	{
		cout << "Can not open file,"<<file <<" in function readFastaMap" << endl;
		return;
	}
	string buf, tmp_id, tmp_seq;
	int pos;
	for(;;)
	{
		getline(Fin, buf);
		if(Fin.eof())
			break;
		else if(buf[0] == '@')
		{
			tmp_id = buf.substr(1);
			getline(Fin, buf);
			tmp_seq = buf;
			out_seq_map.insert(make_pair( tmp_id, tmp_seq ) );
		}
	} // end for
	Fin.close();
}

void DNA::writeFasta ( char* file )
{
  try {
      ofstream Fout(file);
      if(!Fout) throw strcat("Cannot open output Fasta file ", file);

      int i,j;
      for(i=0; i<chain.size(); i++)
	  {
	  Fout << ">" << chain[i].description << endl;
		  for(j=0;j<chain[i].seq.size();j++)
		  {
			  Fout << chain[i].seq[j];
			  if((j+1)%60==0) Fout << endl;
		  }
	  Fout << endl;
      }
	  Fout.close();
  }// end try
  catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}

void writeFastaMap ( char* file, map_DNA& o_seq )
{
  try {
      ofstream Fout(file);
      if(!Fout) throw strcat("Cannot open output Fasta file ", file);

	  map<string, string>::iterator m_it;
	  for(m_it = o_seq.begin(); m_it != o_seq.end(); ++m_it)
	  {
	  Fout << ">" << m_it->first<< endl;
		  for(int j=0;j < m_it->second.size();j++)
		  {
			  Fout << m_it->second.at(j);
			  if((j+1)%60==0) Fout << endl;
		  }
	  Fout << endl;
      }
	  Fout.close();
  }// end try
  catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}

void DNA::readGenBank ( char* file )
{
try {
      ifstream Fin(file);
      if(!Fin) throw strcat("Cannot open input Fasta file ", file);

      bool newCDS = false,cds_begin = false,seq_begin = false;
      string buf,regstr,indices_pos1,indices_pos2,cds_info;
      int newLine;

      CDS tmpCDS;
      fastaDNA tmpfastaDNA;

      // read the first line
      getline(Fin, buf);
      // LOCUS       NT_113796             201709 bp    DNA     linear   CON 29-FEB-2008
      boost::smatch matchstr;
      regstr = "^LOCUS\\s+(\\w+)(_?)(\\d+)\\s+(\\d+)";
      boost::regex expression(regstr);
      if(boost::regex_search(buf, matchstr,expression))
        {
				string tempid(matchstr[0].first, matchstr[0].second);
				tmpfastaDNA.set_id(tempid);
				tmpfastaDNA.set_description( buf );

				string  seq_len(matchstr[4].first, matchstr[4].second);
				int seq_length = boost::lexical_cast<int>(seq_len);
				if (seq_length > 0)
				  	tmpfastaDNA.seq.reserve(seq_length);
        }

	// read the cds infomation
  for(;;) {

	// read a new line in every cycle
	getline(Fin, buf);

	//read all the cds already,break out the for loop
	boost::regex cds_end_tag("ORIGIN");
	if( regex_search(buf, cds_end_tag))
	{
		//cout<<"origin"<<endl;
		seq_begin = true;
		addCDS(tmpCDS);
	}

	 /// ^\s+(\w+)\s+(complement)?\(?<?\d+\.\.\d+\)?
	 string read_type = "^\\s+(\\w+)\\s+(complement)?\\(?<?\\d+\\.\\.\\d+\\)?";
	 boost::regex cds_tag(read_type);
	 boost::regex seq_tag("^\\s+\\d+\\s+(.+)[A-Za-z]$");
	      if (regex_search(buf, cds_tag)) newLine = 1;
   else if (Fin.eof()) newLine = 3;
	 else if (regex_search(buf,seq_tag)) newLine = 4;
	    else newLine = 2;


	if (newLine == 1)
	{
		//cout << " newLine 1" << endl;
		int pos1,pos2;
		boost::regex indices_tag("(\\d+)");
		vector<string> indices_pos;
		string::const_iterator start = buf.begin();
		string::const_iterator end = buf.end();
		while( boost::regex_search(start, end, matchstr, indices_tag) )
		{
		   string msg(matchstr[1].first, matchstr[1].second);
		   indices_pos.push_back(msg);
		   start = matchstr[1].second;
		}
		if ( indices_pos.size() == 2)
		{
		    pos1 = boost::lexical_cast<int>(indices_pos.at(0));
		    pos2 = boost::lexical_cast<int>(indices_pos.at(1));
		    if (pos1 > pos2)
		      myswap(pos1,pos2);
		}

	  if (!newCDS)
		{
		    newCDS = true;
		    cds_begin = true;
		}
	  else
	  {
       // firstly, need to check the current CDS is a new one
       // by compare the pos1, pos2 to the old one
			 // if is still the old one, just complete the reading process;
	      if ( pos1 == tmpCDS.indices.first && pos2 == tmpCDS.indices.second )
	            cds_begin = false;
		  	else
		  	{
	    	//if is a real new CDS, need to store the old chain first,then build a new chain
	     	  addCDS(tmpCDS);
		  		cds_begin = true;
		  	}
		}

		if ( cds_begin )
		{
		  tmpCDS.initial();
		  tmpCDS.set_indices(pos1,pos2);

		  boost::regex complement_tag("complement");
		  if( boost::regex_search(buf,complement_tag) )
		        tmpCDS.recomplement = true;

		  boost::smatch what;
		  boost::regex_match(buf, what, cds_tag);
		  string tmp_str( what[1].first, what[1].second );
		  tmpCDS.cds_info = tmp_str;
		  //cout << "cds_info " << tmp_str << endl;
		}
	} // end newLine == 1

	if (newLine == 2 )
	{
	     string read_tag;

	    /// /gene /locus_tag /protein_id /db_xref
	    if( tmpCDS.gene.empty() )
	    {
	      read_tag = "(/gene=\")(.+)(\")";
	      boost::regex gene_tag(read_tag);
	      boost::smatch matchstr;
	      if(regex_search(buf,matchstr,gene_tag))
	      {
					string tmp_1(matchstr[2].first, matchstr[2].second);
					tmpCDS.gene = tmp_1;
	      }
	    }

	    if ( tmpCDS.locus_tag.empty() )
	    {
	      boost::regex locus_tag_tag("(/locus_tag=\")(.+)(\")");
	      boost::smatch matchstr;
	      if(regex_search(buf,matchstr,locus_tag_tag))
	      {
					string tmp_2(matchstr[2].first, matchstr[2].second);
          tmpCDS.locus_tag = tmp_2;
	      }
	    }

	    if ( tmpCDS.protein_id.empty() )
	    {
	      boost::regex protein_id_tag("(/protein_id=\")(.+)(\")");
	      boost::smatch matchstr;
	      if(regex_search(buf,matchstr,protein_id_tag))
	      {
					 string tmp_3(matchstr[2].first, matchstr[2].second);
           tmpCDS.protein_id = tmp_3;
	      }
	    }

	    if ( tmpCDS.db_xref.empty() )
	    {
	      boost::regex db_xref_tag("(/db_xref=\")(.+)(\")");
	      boost::smatch matchstr;
	      if(regex_search(buf,matchstr,db_xref_tag))
	      {
					string tmp_4(matchstr[2].first, matchstr[2].second);
					tmpCDS.db_xref = tmp_4;
	      }
	    }
	}//end newLine ==2

	// read the sequence
	if( newLine == 4 && seq_begin)
	{
	 	//cout<<"seq:"<<buf<<endl;
		boost::regex s("([A-Za-z]+)");
		string::const_iterator start = buf.begin();
		string::const_iterator end = buf.end();
		while( regex_search(start, end, matchstr, s) )
		{
		   string msg(matchstr[1].first, matchstr[1].second);
		  // cout<<"seq_words:"<<msg<<endl;
		   tmpfastaDNA.seq.append(msg);
		   start = matchstr[1].second;
		}
	}

	if (newLine == 3)
	{
		// read all the infomation already
		//tmpfastaDNA.cleanSeq();
		addChain(tmpfastaDNA);
		break;
	}//end newLine ==3

  }// end for

	cds.erase(cds.begin());
	Fin.close();
  } // end try
  catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}

void DNA::addChain( const fastaDNA& dc)
{
  chain.push_back(dc);
}

void DNA::addCDS(CDS input_cds)
{
   cds.push_back(input_cds);
}

void DNA::addChain(const string& id, const string& description, const string& seq)
{
  fastaDNA tempChain;
  tempChain.set_id(id);
  tempChain.set_description(description);
  tempChain.set_seq(seq);
  chain.push_back(tempChain);
}

int DNA::getChainNum()
{
	return chain.size();
}


int DNA::getCDSNum()
{
	return cds.size();
}

fastaDNA DNA::getCDSSeq( const int& chain_index, const int& cds_index)
{
	fastaDNA tmpfastaDNA;
	string tmp_str= boost::lexical_cast<string>(cds_index);
	tmpfastaDNA.set_id(tmp_str);
	tmpfastaDNA.set_description(tmp_str);
	tmpfastaDNA.description.append("  cds_info=");
	tmpfastaDNA.description.append(cds.at(cds_index).cds_info);

	if(!cds.at(cds_index).gene.empty() )
	{
	   tmpfastaDNA.description.append("  gene=");
	   tmpfastaDNA.description.append(cds.at(cds_index).gene);
	}
	if(!cds.at(cds_index).protein_id.empty() )
	{
	   tmpfastaDNA.description.append("  protein_id=");
	   tmpfastaDNA.description.append(cds.at(cds_index).protein_id);
	}
	if(!cds.at(cds_index).locus_tag.empty() )
	{
	   tmpfastaDNA.description.append("  locus_tag=");
	   tmpfastaDNA.description.append(cds.at(cds_index).locus_tag);
	}
	if(!cds.at(cds_index).db_xref.empty() )
	{
	   tmpfastaDNA.description.append("  db_xref=");
	   tmpfastaDNA.description.append(cds.at(cds_index).db_xref);
  }

	tmp_str =  chain.at(chain_index).seq.substr
		( cds.at(cds_index).indices.first-1,
		 	cds.at(cds_index).indices.second - cds.at(cds_index).indices.first + 1 );
	tmpfastaDNA.set_seq(tmp_str);
	return tmpfastaDNA;
}

fastaDNA DNA::getChain(const int& index)
{
  try{
        int ChainNum = getChainNum();
	if( index > ChainNum) throw "Error: index out,index,array size.";
	return chain.at(index);
  }
 catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}

void DNA::oligoCount(const int& chain_index, const int& oligo_num)
{
  try {
    map<string,int> cur_oligo_counter;
	string tmpSeq = chain.at(chain_index).seq;
	boost::to_upper(tmpSeq);
	int seq_len = tmpSeq.length();

	for (int i = 0 ; i <= seq_len - oligo_num; ++i)
	{
	  string oligo = tmpSeq.substr(i,oligo_num);
	  if (cur_oligo_counter.count(oligo))
		++cur_oligo_counter[oligo];
	  else
		cur_oligo_counter.insert(map<string,int>::value_type(oligo,1));
	} //

	map<string,int>::iterator map_it = cur_oligo_counter.begin();
    for (; map_it != cur_oligo_counter.end(); ++map_it)
	  if (oligo_counter.count(map_it->first))
		   oligo_counter[map_it->first] += map_it->second;
	  else
		  oligo_counter.insert(map<string,int>::value_type(map_it->first, map_it->second));
  }
  catch (char* pMsg) { cerr << endl << "Exception:"<< pMsg << endl; }
}

void DNA::getOligoNum()
{
	cout << "Base information in all sequences " << endl;
	map<string,int>::const_iterator map_it = oligo_counter.begin();
	int all_base_num = 0;

	while ( map_it != oligo_counter.end())
	{
		all_base_num +=  map_it->second;
	    ++map_it;
	 }
	 cout  << " Total sequence: " << chain.size() <<endl;
	 cout << " Total base:" << all_base_num << endl;

	 map_it = oligo_counter.begin();
	while ( map_it != oligo_counter.end())
	 {
	  cout <<"\t"<< map_it->first << "\t"
	  << map_it->second << "\t"
	  << 100.0 * map_it->second / all_base_num
	  << " %" << endl;
	  ++map_it;
	 }

	//--------------------------------------------------
	//  if (oligo_len == 1)
	//  {
	// 	cout << "\tA + T" << "\t"
	// 		 << oligo_counter["A"] + oligo_counter["T"] <<"\t"
	// 		 << 100.0 * (oligo_counter["A"] + oligo_counter["T"] )/ all_base_num
	// 		 << " %"<< endl;
	// 	cout << "\tC + G" << "\t"
	// 		 << oligo_counter["C"] + oligo_counter["G"] <<"\t"
	// 		 << 100.0 * (oligo_counter["C"] + oligo_counter["G"] )/ all_base_num
	// 		 << " %"<< endl;
	//  }
	//--------------------------------------------------
}


bool DNA::getOligo2Id(string input_oligo, vector<fastaDNA>::iterator input_it)
{
 	//regex oligo_tag(input_oligo);
	string tmpSeq = input_it->seq;
	boost::to_upper(tmpSeq);
	int loc = tmpSeq.find(input_oligo,0);
  	if(loc != string::npos)
	{
		//input_it->id.erase(input_it.find_last_not_of(' ')+1, string::npos);
		cout << input_it->id << endl;
		return true;
        }
	else
		return false;
}



/*!
    \fn DNA::get_longerThan(int input_len)
 */

void DNA::readQual(char* file)
{
    try {
      ifstream Fin(file);
      if(!Fin) throw strcat("Cannot open input Fasta file ", file);

      bool fastaFormat = false;
      string buf;
      int newLine;
      SEQID tmp_id;
      QUAL tmp_qual;
      QUAL single_line;

      for(;;) {

	  // read a new line in every cycle
	  	getline(Fin, buf);

  	  /* The new line is marked with an integer, which corresponds to one of the following cases:

     		1. New sequence line, which starts with "> ".
     		2. Sequence line, which starts with any charactor except for ">"
     		3. End of file, detected by function eof().
	  */

         if(buf[0] == '>' ) newLine=1;
	  else if(Fin.eof()) newLine=3;
	  else newLine=2;

          if( newLine == 1 ) { // this is a start of a new chain.
	      if( !fastaFormat ) {
	          // if it is the first chain, just build a new chain object
		  				tmp_id.clear();
		  				tmp_qual.clear();
	          	fastaFormat = true;
	      }
	      else {
	         // otherwise, need to store the old chain first, then build a new chain
	     	  qual.insert(map<SEQID,QUAL>::value_type(tmp_id, tmp_qual));
	     	  tmp_id.clear();
          tmp_qual.clear();
	         }
	      int pos = buf.find_first_of(" ");
	      tmp_id =  buf.substr(1,pos-1);
	  }

	  if( newLine == 2 && fastaFormat )
	  { // qual line
	     /// transform string qual line to vector<int>
	      string tmp_str;
	      int tmp_q;
	     for(int i=0; i<buf.length(); i++)
	     {
	      if ( ! (buf.at(i) == ' '))
	         tmp_str.push_back(buf.at(i));
	       if ( buf.at(i) == ' ' || i == buf.length()-1 )
	        {
	           tmp_q = boost::lexical_cast<int>(tmp_str);
	           tmp_qual.push_back(tmp_q);
	           tmp_str.clear();
	         }
	      }
	  }

	  if( newLine == 3 ) { // end of file
	      if( !fastaFormat ) throw strcat(file, " is not in Fasta format.");
	      else {
	          // otherwise, need to store the old chain first, then to
	          // build a new chain
	          //tmpChain.cleanSeq();
	          qual.insert(map<SEQID,QUAL>::value_type(tmp_id, tmp_qual));
	          break;
	      		}
	  		}
      }
	  Fin.close();
  }

  catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}

void DNA::writeQual(char *file)
{
  /// @todo implement me
  try
  {
      ofstream Fout(file);
      if(!Fout) throw strcat("Cannot open output Fasta file ", file);

      map<SEQID,QUAL>::iterator it;
      int i,j;
      for(it=qual.begin(); it !=qual.end(); it++)
      {
	  Fout << ">" << it->first << "\tlength=" << it->second.size() << endl;
	  for( int j = 0; j < it->second.size(); j++)
	  {
	      Fout << it->second.at(j);
	      if((j+1)%60==0) Fout << endl;
	      else Fout << " ";
	  }
	  Fout << endl;
      }
	  Fout.close();
  }
  catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}

void DNA::Read454AlignDepth(char* file)
{
	/// get the Depth infomation in the contigs
	  try
  {
      ifstream Fin(file);
      if(!Fin) throw strcat("Cannot open output Fasta file ", file);

      bool new_contig = false;
      string buf;
      int tmp_pos_depth,newLine;
      SEQID tmp_id;
      DEPTH tmp_depth;

      map<SEQID,DEPTH>::iterator it;

      for(;;)
      {
      getline(Fin, buf);

       	 	 if(buf[0] == '>' ) newLine=1;
	  	else if(Fin.eof()) newLine=3;
	    else newLine=2;

      if( newLine == 1 )
          { // this is a start of a new chain.
	      		if( !new_contig )
	      			{
	          		// if it is the first chain, just build a new chain object
		  					tmp_depth.clear();
	         		 	new_contig = true;
	      			}
	      		else
	      			{
	         				// otherwise, need to store the old chain first, then build a new chain
		  					 //tmpChain.cleanSeq();
		  					//cout << tmp_id << "\t" << tmp_depth.size()<<endl;
	     	  			depth.insert(map<SEQID,QUAL>::value_type(tmp_id, tmp_depth));
	     	  			tmp_id.clear();
		  					tmp_depth.clear();
	             }
	      			int pos = buf.find_first_of(" ");
	      			tmp_id =  buf.substr(1,pos-1);
	  			}

	  if( newLine == 2 && new_contig )
	  {
     typedef vector< string > split_vector_type;
     split_vector_type SplitVec; // #2: Search for tokens
     boost::split( SplitVec, buf, boost::is_any_of(" \t") );
		 tmp_pos_depth = boost::lexical_cast<int>(SplitVec.at(3));
	   tmp_depth.push_back(tmp_pos_depth);
	   }

	  if( newLine == 3 )
	  	{ // end of file
	      if( !new_contig ) throw strcat(file, " is not in Fasta format.");
	      else
	      	{
	          // otherwise, need to store the old chain first, then to
	          // build a new chain
	      		//cout << tmp_id << "\t" << tmp_depth.size()<<endl;
	          depth.insert(map<SEQID,QUAL>::value_type(tmp_id, tmp_depth));
	          break;
	      	}
	  		}
      } // end for
	  Fin.close();
  }/// end try
  catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}

void DNA::Write454AlignDepth(char* file)
{
  try
  {
      ofstream Fout(file);
      if(!Fout) throw strcat("Cannot open output file ", file);

      map<SEQID,DEPTH>::iterator it;
      int i,j;
      for(it=depth.begin(); it !=depth.end(); it++)
      {

	  		Fout << ">" << it->first << "\tlength=" << it->second.size() << endl;
	  		//copy ( it->second.begin(), it->second.end(),ostream_iterator<int>(cout," "));
	  		for( int j = 0; j < it->second.size(); j++)
	  		{
	      	Fout << it->second.at(j)<< endl;
	  		}
	  		Fout << endl;
      }
	  Fout.close();
  }
  catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}

void DNA::trimBeginEndX(int short_len)
{
  /// @todo implement me
  unsigned int numReads = getChainNum();
  int raw_len = 0, trimed_begin_len = 0 ,trimed_end_len = 0;
  vector<int>::iterator it_begin,it_end;
  string XBegin = "(^X+)([ATCG]+)";
  string EndX = "([ATCG]+)(X+$)";
  boost::regex TrimBegin(XBegin);
  boost::regex TrimEnd(EndX);

  for(int i=0; i < numReads; i++ )
  {
   raw_len = chain.at(i).seq.length();

   /// trim Begin X
    chain.at(i).seq = boost::regex_replace( chain.at(i).seq, TrimBegin, "$2" );
    trimed_begin_len = raw_len - chain.at(i).seq.length();
    it_begin = qual[chain.at(i).id].begin();
    it_end = it_begin + trimed_begin_len;
    if (trimed_begin_len > 0)
      qual[chain.at(i).id].erase(it_begin, it_end);

    /// trim End X
    chain.at(i).seq = boost::regex_replace( chain.at(i).seq, TrimEnd, "$1" );
    trimed_end_len = raw_len - trimed_begin_len - chain.at(i).seq.length();
    it_begin = qual[chain.at(i).id].begin() + chain.at(i).seq.length();
    it_end = qual[chain.at(i).id].end();
    if (trimed_end_len > 0)
      qual[chain.at(i).id].erase(it_begin, it_end);

    cout << raw_len << "\t" << trimed_begin_len << "\t" << trimed_end_len << endl;
    cout << chain.at(i).seq.length() << "\t" << chain.at(i).seq  <<endl;

    cout << qual[chain.at(i).id].size() << "\t";
    for (int j = 0; j < qual[chain.at(i).id].size(); j++)
      cout << qual[chain.at(i).id].at(j) << " " ;
      cout << endl;


    if (chain.at(i).seq.length() < short_len)
        cout <<"After trim:\t"<< chain.at(i).id << "\tshorter than\t" << short_len << endl;
  }
}

void DNA::UniqByDescription()
{
	sort(chain.begin(), chain.end(), less_by_description());
	vector<fastaDNA>::iterator it;
	it = unique( chain.begin(), chain.end(), equal_by_description());
	chain.resize( it - chain.begin());
}

void DNA::clear()
{
	chain.clear();
	cds.clear();
	oligo_counter.clear();
	qual.clear();
	depth.clear();
}

void getRandSeq(int input_len,int input_coverage, DNA &input_dna,DNA *out_dna)
{

  fastaDNA tmpChain;
  for( int i = 0; i < input_dna.chain.size(); i++){
  int seq_len = input_dna.chain.at(i).length();
  int rand_seq_number = seq_len * input_coverage / input_len;
  cout << "rand_seq_number" << rand_seq_number <<endl;

  int * rand_num_ptr = getRandNum(seq_len, rand_seq_number);

    for( int j = 0; j <= rand_seq_number; j++ ){

    int pos = *(rand_num_ptr + j);

    srand(pos);
    int adjust_len = rand()%10 + 1; // to adjust the random sequences's length

    int plusORminus = rand()%2;
    if (plusORminus == 1)
        adjust_len = -1 * adjust_len;

    int rand_seq_len = input_len + adjust_len;
    //cout << "rand_seq_len"<< rand_seq_len << endl;
      if ( (seq_len - pos) >= input_len){

          tmpChain.seq = input_dna.chain.at(i).seq.substr(pos,rand_seq_len);

          /// set the seq reverse_compliment randomly
          if(plusORminus == 1)
            tmpChain.reverse_compliment();

	  /// set id
          tmpChain.id = input_dna.chain.at(i).id;
          tmpChain.id.append("_rand_");
          string tmpstr = boost::lexical_cast<string>(j);
          tmpChain.id.append(tmpstr);
          /// set description
          tmpChain.description = tmpChain.id.append("  length=");
          tmpstr = boost::lexical_cast<string>(rand_seq_len);
          tmpChain.description = tmpChain.id.append(tmpstr);
          tmpChain.description = tmpChain.id.append("  position=");
          tmpstr = boost::lexical_cast<string>(pos);
          tmpChain.description = tmpChain.id.append(tmpstr);

          /// save rand seq
          out_dna->chain.push_back(tmpChain);
      }//end if
    }//end for j
  }// end for i
}
