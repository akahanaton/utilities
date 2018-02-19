#include "/home/gmswenm/myProgram/myHeader/general.h"
#include "/home/gmswenm/myProgram/myHeader/DNA.h"
#include "/home/gmswenm/myProgram/myHeader/blastparser.h"

using namespace std;

void usage();

      int opt;
	  int query_len_cutoff =0, sbjct_len_cutoff =0, align_len_cutoff = 0, identity_cutoff = 0, results_number = 1000;
	  bool get_gap_region = false, get_contig_cov_info = false, get_est_cov_info = false, show_blast_result = false, sort_flag = false;
      char *in_seq = 0, *blast_file = 0, *ref_seq = 0, *out_seq = 0, *in_align = 0, *out_depth = 0;


int main(int argc, char *argv[])
{
	if(argc == 1)
	{
		usage();
		return EXIT_SUCCESS;
	}
	else
	{
		while((opt = getopt(argc, argv, "f:r:b:n:i:o:a:l:q:s:vwceS")) != -1 )
		{
			  switch(opt)
			  {
				case 'f':
					in_seq = optarg;
					break;
				case 'r':
					ref_seq = optarg;
					break;
				case 'o':
					out_seq = optarg;
					break;
				case 'a':
					in_align = optarg;
					break;
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
				case 'q':
				  query_len_cutoff = boost::lexical_cast<int>(optarg);
				  break;
				case 's':
				  sbjct_len_cutoff = boost::lexical_cast<int>(optarg);
				  break;
				case 'w':
				  show_blast_result = true;
				  break;
				case 'S':
				  sort_flag = true;
				  break;
				case 'v':
				  get_gap_region = true;
				  break;
				case 'c':
				  get_contig_cov_info = true;
				  break;
				case 'e':
				  get_est_cov_info = true;
				  break;
				case '?':
					printf("unknown option: %c\n", optopt);
					return EXIT_SUCCESS;
				break;
			}
		} // end while
	}// end if-else

  /// map contig to ref
  BlastParser bp;
  //--------------------------------------------------
  // printf("unknown option: %s\n", "here2");
  //--------------------------------------------------
  bp.ReadBlastResult(blast_file, results_number);
  //--------------------------------------------------
  // printf("unknown option: %s\n", "here3");
  //--------------------------------------------------
  if(get_contig_cov_info)
  {
	bp.GetContigCoverInfo(query_len_cutoff, sbjct_len_cutoff,  align_len_cutoff, identity_cutoff, sort_flag);
	  if(show_blast_result)
		  bp.GetResult( false);
	return EXIT_SUCCESS;
  }
  if(get_est_cov_info){
      //--------------------------------------------------
      // printf("unknown option: %s\n", "here1");
      //--------------------------------------------------
	bp.GetEstCoverInfo(query_len_cutoff, sbjct_len_cutoff,  align_len_cutoff, identity_cutoff, sort_flag);
	  if(show_blast_result)
		  bp.GetResult( false);
		return EXIT_SUCCESS;
  }
  if(get_gap_region){
		;
  }else{
		bp.SetRefSeq(ref_seq);
		bp.SetReads(in_seq);
		bp.MapContigsToRef(out_seq, align_len_cutoff, get_gap_region, in_align, out_depth);
  }
  return EXIT_SUCCESS;
}


void usage()
{
  cout << "Map assembly contigs to reference genome." << endl;
  cout << "According to the blast resualt,map the contigs to reference genome,uncovered region stand by 'X'" << endl;
  cout << "\tOptions used for get coverage info:" << endl;
  cout << "\t\t-b STR:  blast result file"<<endl;
  cout << "\t\t-i INT:  identity of alingment, default 0"<<endl;
  cout << "\t\t-l INT:  alingment length threshold use for map contig to reference, default 0"<<endl;
  cout << "\t\t-q INT:  query length threshold use for filter blast result, default 0"<<endl;
  cout << "\t\t-s INT:  sbjct length threshold use for filter blast result, default 0"<<endl;
  cout << "\t\t-S 	 :  sort flag: sort by query id and sbjct id, only for parted database"<<endl;
  cout << "\t\t-w    :  show the actual blast result used for map"<<endl;
  cout << "\t\t-c    :  get the coverage infomation"<<endl;
  cout << "\t\t-e    :  get the est coverage infomation, one query in on  contig"<<endl;
  cout << "\t\t-n INT:  Estimated results number, default 1000, usefull for speed up"<<endl;
  cout << "\tOptions additional, used for get mapped sequence:"<<endl;
  cout << "\t\t-f STR:  contigs sequences,fasta format"<<endl;
  cout << "\t\t-r STR:  reference sequences, fasta format"<<endl;
  cout << "\t\t-o STR:  out sequences, fasta format, optional"<<endl;
  cout << "\t\t-v    :  reversely get the uncovered region "<<endl;
}
