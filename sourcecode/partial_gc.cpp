#include "/share/raid5/wenming/myProgram/myHeader/DNA.h"

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
      int  n_times = 1;
      char *in_seq = 0;
  while((opt = getopt(argc, argv, "f:n:")) != -1 )
  {
     switch(opt)
      {
      case 'f':
			in_seq = optarg;
			break;
      case 'n':
			n_times = boost::lexical_cast<int>(optarg);
			break;
      case '?':
			printf("unknown option: %c\n", optopt);
			return EXIT_SUCCESS;
      }
   }

  if( !in_seq)
  {
    cout << "File read in error, please check the usage. "<< endl;
    return EXIT_SUCCESS;
  }
  else
  {
	DNA mydna;
	mydna.readFasta(in_seq);
	int seq_num = mydna.chain.size();
	for(int i = 0; i < seq_num; ++i)
		mydna.oligoCount(i, 1);
	mydna.getOligoNum();

	int all_base = 0;
	map<string, int>::iterator map_it = mydna.oligo_counter.begin();
	while(map_it != mydna.oligo_counter.end())
	{
		all_base += map_it->second;
		++map_it;
	}

	float all_gc_content = 100.00 * (mydna.oligo_counter["C"] + mydna.oligo_counter["G"]) / all_base;
	mydna.oligo_counter.clear();
	cout << "Tolal GC_content  " << all_gc_content << endl;
	cout << "Random select sequences ......"<<endl;
	 cout << "Sequence_number\tGC_content"<< endl;
	 for(int k = 0; k < n_times; ++k)
	{
	for(int i = 5; i < seq_num; ++i)
	{
		int* rand_num_ptr = getRandNum(seq_num -1 , i);

		for(int j = 0; j < i; ++j)
		{
			int chain_index = *(rand_num_ptr+j);
			if( chain_index == 0)
				cout << chain_index << "  ";
			mydna.oligoCount(chain_index, 1);
		}
		int cur_all_base = 0;
		map_it = mydna.oligo_counter.begin();
		while(map_it != mydna.oligo_counter.end())
		{
			cur_all_base += map_it->second;
			++map_it;
		}
		float cur_gc_content = 100.00 * (mydna.oligo_counter["C"] + mydna.oligo_counter["G"]) / cur_all_base;
		 if( abs(cur_gc_content - all_gc_content) < 0.1)
		 {
			 cout << i <<"\t"<< cur_gc_content << endl;
		 }

		delete [] rand_num_ptr;
		mydna.oligo_counter.clear();
	}
	cout << endl;
  }
   }//end if else
}

  return EXIT_SUCCESS;
}

void usage()
{
  cout << "Options:" << endl;
  cout << "\t\t-f STR: sequences file"<<endl;
  cout << "\t\t-n INT: repeat the calculation n times"<<endl;
}
