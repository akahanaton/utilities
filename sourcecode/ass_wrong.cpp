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
      int align_len_cutoff = 0, identity_cutoff = 0, results_number = 1000;
      char *blast_file = 0 ;
  while((opt = getopt(argc, argv, "b:n:i:l:")) != -1 )
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
   bp.AssemblyWrongCheck(align_len_cutoff, identity_cutoff);
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
}
