#include "/share/raid5/wenming/myProgram/myHeader/DNA.h"
#include "/share/raid5/wenming/myProgram/myHeader/blastparser.h"

using namespace std;
void usage();

int main(int argc, char *argv[])

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
  while((opt = getopt(argc, argv, "b:n:i:l:c:m:huspq")) != -1 )
  {
     switch(opt)
      {
      case 'b':
			blast_file = optarg;
			break;
      case 'n':
			results_number = boost::lexical_cast<int>(optarg);
			break;
      case 'c':
			align_ratio_cutoff = boost::lexical_cast<int>(optarg);
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
      case 'h':
			head_info = true;
			break;
      case 'u':
			uniq_info = true;
			break;
      case 's':
			sort_by_SbjctID_SStart = true;
			break;
      case 'p':
			repeat_remove = true;
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

   if(uniq_info)
	   bp.UniqBy_QueryID();

   if(identity_cutoff > 0)
		bp.RmLowIdentity( identity_cutoff );

   if(align_len_cutoff > 0)
	   bp.RmShortAlign(align_len_cutoff);

   if ( (align_ratio_cutoff > 0) && ( !sort_by_SbjctID_SStart && !sort_by_QueryID_QStart ))
   {
	   sort(bp.BlastResults.begin(), bp.BlastResults.end(), less_by_QueryID_QStart());
	   bp.RmLowQueryAlignRatio(align_ratio_cutoff);
   }
   if(sort_by_QueryID_QStart)
   {
	   sort(bp.BlastResults.begin(), bp.BlastResults.end(), less_by_QueryID_QStart());
	   if( align_ratio_cutoff > 0)
		   bp.RmLowQueryAlignRatio(align_ratio_cutoff);
   }
   if(sort_by_SbjctID_SStart)
   {
	   sort(bp.BlastResults.begin(), bp.BlastResults.end(), less_by_SbjctID_SStart());
	   if( align_ratio_cutoff > 0)
		   bp.RmLowSbjctAlignRatio(align_ratio_cutoff);

   }
   if(mussy_len_cutoff > 0)
	   if(sort_by_QueryID_QStart) // this operation muset be called after sort_by_QueryID_QStart, and this flag  imply the results have been sorted;
			bp.RmMussyAlignInContig(mussy_len_cutoff);
	   else   // otherwise, sort first
	   {
		sort(bp.BlastResults.begin(), bp.BlastResults.end(), less_by_QueryID_QStart());
		bp.RmMussyAlignInContig(mussy_len_cutoff);
	   }


   if(repeat_remove)
	   bp.RmRepeatAlignBy_QStart();

   bp.GetResult(head_info);
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
  cout << "\t\t-c INT: Integrity Coverage ratio threshold, work with -q (query) -s (subjcet). default:-q , 0 "<<endl;
  cout << "\t\t-l INT: Alingment length cutoff, default 0"<<endl;
  cout << "\t\t-m INT: Rm mussy alignment in blast result, useful for long sequence alignment cleaning"<<endl;
  cout << "\t\t-h    : Output the comment info at the head "<<endl;
  cout << "\t\t-u    : Uniq by query id: one result for one query"<<endl;
  cout << "\t\t-s    : Sort by Sbjct start position"<<endl;
  cout << "\t\t-p    : Revmove rePeat by QStart and QEnd"<<endl;
  cout << "\t\t-q    : Sort by Query start position"<<endl;
}
