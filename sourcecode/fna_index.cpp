#include "/share/raid6/wenming/myProgram/myHeader/DNA.h"

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
      char *in_seq = 0,  *out_seq = 0, *prefix = 0;
  while((opt = getopt(argc, argv, "f:o:p:")) != -1 )
  {
     switch(opt)
      {
      case 'f':
			in_seq = optarg;
			break;
      case 'o':
			out_seq = optarg;
			break;
      case 'p':
			prefix = optarg;
			break;
      case '?':
			printf("unknown option: %c\n", optopt);
			return EXIT_SUCCESS;
      }
   }// endl while

DNA mydna;
mydna.readFasta(in_seq);
int dna_num = mydna.chain.size();
string index_name;
if(prefix)
	for(int i = 0; i < dna_num; ++i)
	{
		index_name.append(prefix);
		index_name.append(boost::lexical_cast<string>(i));
		mydna.chain.at(i).description = index_name;
		index_name.clear();
	}
else
	for(int i = 0; i < dna_num; ++i)
	{
		index_name = boost::lexical_cast<string>(i);
		mydna.chain.at(i).description.insert(0,index_name);
	}
	mydna.writeFasta(out_seq);
 }// endl if
  return EXIT_SUCCESS;
}

void usage()
{
  cout << "Index the sequences" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f STR: input sequences file"<<endl;
  cout << "\t\t-o STR: output sequences file"<<endl;
  cout << "\t\t-p STR: the prefix of the index name of the sequences"<<endl;
}
