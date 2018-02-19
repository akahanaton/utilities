#include "/share/raid5/wenming/myProgram/myHeader/DNA.h"
#include "/share/raid5/wenming/myProgram/myHeader/blastparser.h"

struct equal_by_QueryID_QEnd
{
  /// pre: BlastResults are not empty
  /// post: all reads mapped to different locs
  bool operator() (const BLAST8RESULT &r1, const BLAST8RESULT &r2 ) const
  {
    if ( r1.QueryID == r2.QueryID )
      {
        /// a read mapped to different locs
        return ( r2.QEnd <= r1.QEnd );
      }
    else
		return false;
  }
};

using namespace std;
void usage();

int main(int argc, char *argv[])
{
if(argc == 1)
{
	usage();
}
else
 {
      int opt;
      int align_len_cutoff = 0, align_ratio_cutoff = 0, mussy_len_cutoff = 0, identity_cutoff = 0, results_number = 1000;
	  bool is_reads = false, head_info = false, uniq_info = false, sort_by_QueryID_QStart = false, repeat_remove = false;
	  bool sort_by_SbjctID_SStart = false;
      char *in_seq = 0,  *blast_file = 0 ;
  while((opt = getopt(argc, argv, "b:n:i:l:m:q")) != -1 )
  {
     switch(opt)
      {
      case 'b':
			blast_file = optarg;
			break;
      case 'n':
			results_number = boost::lexical_cast<int>(optarg);
			break;
      case 'i':
			identity_cutoff = boost::lexical_cast<int>(optarg);
			break;
      case 'l':
			align_len_cutoff = boost::lexical_cast<int>(optarg);
			break;
      case 'm':
			mussy_len_cutoff = boost::lexical_cast<int>(optarg);
			break;
      case 'q':
			sort_by_QueryID_QStart = true;
			break;
      case '?':
			printf("unknown option: %c\n", optopt);
			return EXIT_SUCCESS;
      }
   }

  if( !blast_file)
  {
    cout << "File read in error, please check the usage. "<< endl;
    return EXIT_SUCCESS;
  }
  else{
  ///blast result filter
   BlastParser bp;
   bp.ReadBlastResult(blast_file, results_number);

   //cout << bp.BlastResults.size() << endl;


   if(identity_cutoff > 0)
		bp.RmLowIdentity( identity_cutoff );

   if(align_len_cutoff > 0)
	   bp.RmShortAlign(align_len_cutoff);

   if(sort_by_QueryID_QStart)
	   sort(bp.BlastResults.begin(), bp.BlastResults.end(), less_by_QueryID_QStart());

   if(mussy_len_cutoff > 0)
	   if(sort_by_QueryID_QStart) // this operation muset be called after sort_by_QueryID_QStart, and this flag  imply the results have been sorted;
			bp.RmMussyAlignInContig(mussy_len_cutoff);
	   else   // otherwise, sort first
	   {
		sort(bp.BlastResults.begin(), bp.BlastResults.end(), less_by_QueryID_QStart());
		bp.RmMussyAlignInContig(mussy_len_cutoff);
	   }

   //cout << bp.BlastResults.size() << endl;


   vector<BLAST8RESULT>::iterator cur_it;
   cur_it =  unique(bp.BlastResults.begin(), bp.BlastResults.end(), equal_by_QueryID_QEnd());
   bp.BlastResults.resize(cur_it - bp.BlastResults.begin());

   cout << bp.BlastResults.size() << endl;

  bp.GetResult(false);

  bp.GetContigCoverInfo(0,0);

	// BAC    (BAC length)    (number of Contig)    (list of Contig)    (list of Contig length)    (length of all contigs)

  set<string> BAC_ID;
  map<string, int> BAC_len, contig_len;
  map<string, int>::iterator m_it;
  vector<BLAST8RESULT> same_BAC;
  vector<BLAST8RESULT>::iterator same_BAC_it;

  for(cur_it = bp.BlastResults.begin(); cur_it != bp.BlastResults.end(); ++cur_it)
  {
	  BAC_ID.insert(cur_it->QueryID);
	  m_it = BAC_len.find(cur_it->QueryID);
	  if(m_it == BAC_len.end())
		  BAC_len.insert(pair<string, int>(cur_it->QueryID, cur_it->QLen));
  }

  cout << BAC_ID.size() << "\t" << BAC_len.size() << endl;

  set<string>::iterator s_it;
	  for(s_it = BAC_ID.begin(); s_it != BAC_ID.end(); ++s_it)
	  {
		  for(cur_it = bp.BlastResults.begin(); cur_it != bp.BlastResults.end(); ++cur_it)
		  {
			  if(cur_it->QueryID == *s_it)
			  {
				  same_BAC.push_back(*cur_it);
			  }
		  }
		  for(same_BAC_it= same_BAC.begin(); same_BAC_it != same_BAC.end(); ++same_BAC_it)
		  {
			  m_it = contig_len.find(same_BAC_it->SubjectID);
			  if(m_it == contig_len.end())
				  contig_len.insert(pair<string, int>(same_BAC_it->SubjectID, same_BAC_it->SLen));
		  }
		  int total_contig_len = 0;
		 for(m_it = contig_len.begin(); m_it != contig_len.end(); ++m_it)
			 total_contig_len += m_it->second;
		cout << *s_it << "\t" << BAC_len[*s_it] << "\t" << contig_len.size() << "\t"<<total_contig_len << endl;

		 for(m_it = contig_len.begin(); m_it != contig_len.end(); ++m_it)
				cout << "\t" << m_it->first << "\t" << m_it ->second << endl;
		contig_len.clear();
		  same_BAC.clear();
	  }
   }
}
  return EXIT_SUCCESS;
}

void usage()
{
  cout << "Convert blast result into tab format" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-b STR: Blast result file"<<endl;
  cout << "\t\t-n INT: Estimated results number, default 1000, usefull for speed up"<<endl;
  cout << "\t\t-i INT: Identity, default 0"<<endl;
  cout << "\t\t-l INT: Alingment length cutoff, default 0"<<endl;
  cout << "\t\t-m INT: Rm mussy alignment in blast result, useful for long sequence alignment cleaning"<<endl;
  cout << "\t\t-q    : Sort by Query start position"<<endl;
}
