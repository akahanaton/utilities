#ifndef BLASTPARSERH
#define BLASTPARSERH
#include "DNA.h"
/**
	@author wenming <wenming@genomics.org.cn>
*/

typedef string SEQID;
typedef string SEQ;

class single_gap
{
public:
	int	real_gap_size;
	int	gap1_start_pos;
	int gap1_end_pos;
	int	gap2_start_pos;
	int gap2_end_pos;
	DNA	long_reads;
	DNA	short_reads;

	single_gap()
	{
		gap1_start_pos = gap1_end_pos = gap2_start_pos = gap2_end_pos = real_gap_size = 0;
	}
};

/////////////////////// class blastresult  ////////////////////
class BLAST8RESULT
{
//# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
public:
  string  QueryID;
  string  SubjectID;
  int	  Identity;
  int     AlignmentLength;
  int     Mismatches;
  int     Gaps;
  int     QLen;
  int     QStart;
  int     QEnd;
  int     SLen;
  int     SStart;
  int     SEnd;
  //double  Expect;
  float   BitScore;
  bool    RComplement;

  void Init()
  {
   QueryID = "";
   SubjectID = "";
   Identity = 0;
   AlignmentLength = 0;
   Mismatches = 0;
   Gaps = 0;
   QLen = 0;
   QStart = 0;
   QEnd = 0;
   SLen = 0;
   SStart = 0;
   SEnd = 0;
   //Expect = 0;
   BitScore = 0;
   RComplement = false;
  }

  BLAST8RESULT(){Init();}

};

struct less_by_SbjctID_SStart
{
  inline bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2) const
  {
	  if(r1.SubjectID != r2.SubjectID)
		  return r1.SubjectID < r2.SubjectID;
	  else
		  return r1.SStart < r2.SStart;
  }
};

struct less_by_QueryID_QStart
{
	// map contig
 bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2) const
    {
		if (r1.QueryID != r2.QueryID)
			return r1.QueryID < r2.QueryID;
		else
			return r1.QStart < r2.QStart;
    }
};

struct check_align_len:public std::binary_function<BLAST8RESULT,int,bool>
{
inline	bool operator()(const BLAST8RESULT &r, int align_len) const
	{
		return r.AlignmentLength < align_len;
	}
};

struct check_qurey_len:public std::binary_function<BLAST8RESULT,int,bool>
{
inline	bool operator()(const BLAST8RESULT &r, int query_len) const
	{
		return r.QLen < query_len;
	}
};

struct check_sbjct_len:public std::binary_function<BLAST8RESULT,int,bool>
{
inline	bool operator()(const BLAST8RESULT &r, int sbjct_len) const
	{
		return r.SLen < sbjct_len;
	}
};

struct check_query_in_middle:public std::unary_function< vector<BLAST8RESULT>::iterator, bool>
{
inline	bool operator()(const vector<BLAST8RESULT>::iterator r) const
	{
		vector<BLAST8RESULT>::iterator pre;
		vector<BLAST8RESULT>::iterator post;
		pre=r-1;
		post=r+1;
		return (r->QueryID == pre->QueryID && r->QueryID == post->QueryID);
	}
};

struct check_query_align_ratio:public std::binary_function<BLAST8RESULT, int, bool>
{
  /// pre: 1 blastResult are not empty,
  /// post: return ture if the align_ratio < 90.0
  bool operator() (const BLAST8RESULT &r, int align_ratio ) const
  {
        int tmp_ratio = 100 * r.AlignmentLength / r.QLen ;
          return tmp_ratio < align_ratio;
  }
};

struct check_sbjct_align_ratio:public std::binary_function<BLAST8RESULT, int, bool>
{
  /// pre: 1 blastResult are not empty,
  /// post: return ture if the align_ratio < 90.0
  bool operator() (const BLAST8RESULT &r, int align_ratio ) const
  {
        int tmp_ratio = 100 * r.AlignmentLength / r.SLen ;
          return tmp_ratio < align_ratio;
  }
};

struct check_align_len_ratio:public std::binary_function<BLAST8RESULT,int,bool>
{
inline	bool operator()(const BLAST8RESULT &r, int align_len) const
	{
		if (r.AlignmentLength < align_len)
		   {
			   int align_ratio = 100 * r.AlignmentLength / r.QLen;
			   if(align_ratio > 90)
				   return false;
			   else
				   return true;
		   }
		else
		  return false;
	}
};

struct check_align_identity:public std::binary_function<BLAST8RESULT,int,bool>
{
inline	bool operator()(const BLAST8RESULT &r, int align_identity) const
	{
		return r.Identity < align_identity;
	}
};

struct rm_repeat_by_SStart_contig
{
 /// pre:  1.  all reads must  have been sorted by SStart
 /// post: rm the reads maped to the same loc, or inside a first alignment
  bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2) const
  {
    if( r1.QueryID == r2.QueryID )
	{
		//if (r2.RComplement)
		return ( (r2.SStart >= r1.SStart) && (r2.SEnd <= r1.SEnd) );// r2 is inside the r1
        //else
			//return (r2.SEnd <= r1.SEnd) ;// r2 is inside the r1
	}
    else
        return false ;
  }
};

struct rm_repeat_by_QStart_contig
{
 /// pre:  all reads must  have been sorted by SStart
 /// post: rm the same region map to different loc in refence
  bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2) const
  {
    if( r1.SubjectID == r1.SubjectID  )
        return (r1.QStart == r2.QStart && r1.QEnd == r2.QEnd) ;// r2 is inside the r1, remove also
    else
        return false;
  }
};

struct rm_repeat_by_SStart_read
{
 /// pre:  1.  all reads must  have been sorted by SStart
 ///       2.  reads are uniqued by QueryID
 /// post: rm the reads maped to the same loc, or inside a first alignment
  bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2) const
  {
    if( (r2.SStart - r1.SStart ) > 2 )
        return (r2.SEnd < r1.SEnd) ;// r2 is inside the r1, remove also
    else
        return true ;
  }
};

struct equal_by_QueryID
{
  inline bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2 ) const
  {
    return ( r1.QueryID == r2.QueryID );
  }
};

struct less_by_QueryID
{
  inline bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2 ) const
  {
    return ( r1.QueryID <= r2.QueryID );
  }
};

struct less_by_QueryID_SbjctID
{
  inline bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2 ) const
  {
    if ( r1.QueryID != r2.QueryID )
		return ( r1.QueryID <= r2.QueryID );
	else
		return (r1.SubjectID <= r2.SubjectID );
  }
};

struct equal_by_SubjectID
{
  inline bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2 ) const
  {
    return ( r1.SubjectID == r2.SubjectID);
  }
};

struct equal_by_QueryID_SStart
{
  /// pre: BlastResults are not empty
  /// post: all reads mapped to different locs
  bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2 ) const
  {
    if ( r1.QueryID == r2.QueryID )
      {
        /// a read mapped to different locs
        return ( r2.SStart < r1.SEnd );
      }
    else
	return false;
  }
};

struct equal_by_QueryID_QStart
{
  /// pre:  BlastResults are not empty,
  /// post:
  bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2 ) const
  {
    if ( r1.QueryID == r2.QueryID )
      {
        if ( (r2.QStart >= r1.QStart) && (r2.QEnd < r1.QEnd)  ||  (r2.QStart > r1.QStart) && (r2.QEnd <= r1.QEnd)  )
				return true;
		else
			return false;
      }
    else
		return false;
  }
};

///////////////// class blastParser ////////////////////

class BlastParser
{
public:
  vector<BLAST8RESULT>      BlastResults;
  map<SEQID, SEQ> 			Reads, RefSeq;
  bool 						reads_set, ref_set;
public:
    BlastParser();
    ~BlastParser();
    void MapReadsToRef(char *out_filename, float input_identity);
    void AddQuery(const BLAST8RESULT& input_result);
    void ReadBlastResult(char *resultfile, const int& results_num);
    void UniqBy_QueryID( );
    void UniqBy_SubjectID( );
    void UniqBy_QueryID_SStart( );
    void SetReads ( char* infile);
    void SetRefSeq ( char* infile);
    void RmRepeatReads( int inter = 2 );
    void MapContigsToRef(char *out_seq, int align_len_threshold, bool gap_region, char *align_file = 0, char *out_align_depth = 0);
	void GetUncoverSeq(char* out_seq);
	void GetEstCoverInfo(const int& query_len_threshold,const int& sbjct_len_threshold,const int& align_len_threshold, const int& identity_threshold, bool sort_tag);
	void GetContigCoverInfo(const int& query_len_threshold,const int& sbjct_len_threshold,const int& align_len_threshold, const int& identity_threshold, bool sort_tag);
	void RmLowIdentity(int identity_threshold);
    void RmShortAlign(int len_threshold);
    void RmShortQuery(int len_threshold);
    void RmShortSbjct(int len_threshold);
	void RmLowQueryAlignRatio( int ratio_threshold );
	void RmLowSbjctAlignRatio( int ratio_threshold );
	void RmMussyAlignInContig(int len_threshold);
    void RmRepeatAlignBy_SStart();
    void RmRepeatAlignBy_QStart();
    void RmQueryInMiddle();
    void FilterBlastResult(bool contig_flag, bool sort_flag, int len_threshold, int ratio_threshold);
	void GetGapRegion( vector<single_gap>& gr, int add_len = 100);
	void AssemblyWrongCheck(int len_threshold, int identity_threshold);
    void GetResult(bool head_info);
};

#endif
