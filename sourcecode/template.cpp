#include "/share/raid5/wenming/myProgram/myHeader/DNA.h"
#include "/share/raid5/wenming/myProgram/myHeader/blastparser.h"

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
      int align_len_cutoff = 200;
      int align_ratio_cutoff = 0;
	  bool is_reads = false;
      float identity = 90.0;
      char *in_seq = 0,  *blast_file = 0 ;
  while((opt = getopt(argc, argv, "f:b:i:l:a:r")) != -1 )
  {
     switch(opt)
      {
      case 'f':
			in_seq = optarg;
			break;
      case 'b':
			blast_file = optarg;
			break;
      case 'a':
			align_ratio_cutoff = boost::lexical_cast<int>(optarg);
			break;
      case 'i':
			identity = boost::lexical_cast<float>(optarg);
			break;
      case 'l':
			align_len_cutoff = boost::lexical_cast<int>(optarg);
			break;
      case 'r':
			is_reads = true;
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
   if(in_seq)
		bp.SetReads(in_seq);
   bp.ReadBlast8Result(blast_file, identity);
   bp.FilterBlast8Result(is_reads, align_len_cutoff, align_ratio_cutoff);
   if(in_seq)
		bp.GetResult(true);
   else
		bp.GetResult(false);

   }
}

  return EXIT_SUCCESS;
}

void usage()
{
  cout << "filter blast result" << endl;
  cout << "Use Threshold: Identity defalt 95 changable" << endl;
  cout << "               coverage defalt 90 unchangable" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f STR: sequences file"<<endl;
  cout << "\t\t-r    : r stand for reads, otherwise is contig "<<endl;
  cout << "\t\t-b STR: blast result file, m 8 format"<<endl;
  cout << "\t\t-a FLOAT: alignment ratio threshold, default 0"<<endl;
  cout << "\t\t-i INT: Identity, default 90.0"<<endl;
  cout << "\t\t-l INT: alingment length cutoff, default 200"<<endl;
}
