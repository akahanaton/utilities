#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <set>
#include <unistd.h>
#include <sys/stat.h>

#include "../myHeader/DNA.h"

struct solar_result
{
	string q_id;
	int q_len;
	int q_start;
	int q_end;
	string strand;
	string s_id;
	int s_len;
	int s_start;
	int s_end;
	int total_score;
	string q_pos_info;
	string s_pos_info;
	string score_info;
};

struct less_by_QID_QStart
{
  inline bool operator() (const solar_result& r1, const solar_result& r2) const
  {
	  if(r1.q_id != r2.q_id)
		  return r1.q_id < r2.q_id;
	  else
		  return r1.q_start< r2.q_start;
  }
};

using namespace std;

void usage();
void DealSolarResualt( char* in_file, const int& a_len, set<string>& q_id, multimap<string, solar_result>& r);

int main(int argc, char *argv[])
{
	if(argc == 1)
		usage();
	else
	{
      clock_t start,finish;
      double totaltime;
      start=clock();

      int opt, align_len = 0;
      char *scaftig_file = 0, *solar_file = 0;
	  while((opt = getopt(argc, argv, "s:l:")) != -1 )
	  {
		 switch(opt)
		  {
			  case 's':
				solar_file = optarg;
				break;
			  case 'l':
				align_len = boost::lexical_cast<int>(optarg);
				break;
			  case '?':
				 printf("unknown option: %c\n", optopt);
				 return EXIT_SUCCESS;
				 break;
		  }
	   }


		finish=clock();
		totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
		cout<<"\nProgram takes "<<totaltime<<" seconds"<<endl;
		set<string> query_id;
		multimap<string, solar_result> solar_results;
		DealSolarResualt( solar_file, align_len, query_id, solar_results);

		cout << query_id.size()<<endl;
		cout << solar_results.size() << endl;

		set<string>::iterator s_it;
		multimap<string, solar_result>::iterator mm_it;
		pair<multimap<string, solar_result>::iterator, multimap<string, solar_result>::iterator> ret;

		vector<solar_result> vec_results;

		for(s_it = query_id.begin(); s_it != query_id.end(); ++s_it)
		{
			ret = solar_results.equal_range(*s_it);
			for(mm_it = ret.first; mm_it != ret.second; ++mm_it)
				vec_results.push_back(mm_it->second);
			sort(vec_results.begin(), vec_results.end(), less_by_QID_QStart());
			for(int i = 0; i < vec_results.size(); ++i)
				cout << vec_results.at(i).q_id << "\t"<< vec_results.at(i).q_start << endl;
		}
		finish=clock();
		totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
		cout<<"\nProgram takes "<<totaltime<<" seconds"<<endl;
	}// end else

  return EXIT_SUCCESS;
}

void usage()
{
  cout << "Get the contig coverage info according to solar results" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-s STR: soap result file "<<endl;
  cout << "\t\t-l INT: minimum alignment length"<<endl;
}

void DealSolarResualt( char* in_file, const int& a_len, set<string>& q_id, multimap<string, solar_result>& r)
{
	try
	{
		ifstream Fin(in_file);
		if(!Fin) throw strcat( "Cannot open input result file", in_file);

		string buf;
		solar_result tmp_s_result;
		vector<string> split_vec;

		for(;;)
		{
			getline(Fin, buf);
			if(Fin.eof())
				break;
			else
			{
				boost::split( split_vec, buf, boost::is_any_of("\t") );
					tmp_s_result.q_id     = split_vec.at(0);
					tmp_s_result.q_len    = boost::lexical_cast<int>(split_vec.at(1));
					tmp_s_result.q_start  = boost::lexical_cast<int>(split_vec.at(2));
					tmp_s_result.q_end    = boost::lexical_cast<int>(split_vec.at(3));
					tmp_s_result.strand   = split_vec.at(4);
					tmp_s_result.s_id     = split_vec.at(5);
					tmp_s_result.s_len    = boost::lexical_cast<int>(split_vec.at(7));
					tmp_s_result.s_start  = boost::lexical_cast<int>(split_vec.at(8));
					tmp_s_result.s_end    = boost::lexical_cast<int>(split_vec.at(9));
					tmp_s_result.total_score    = boost::lexical_cast<int>(split_vec.at(10));
					tmp_s_result.q_pos_info = split_vec.at(11);
					tmp_s_result.s_pos_info = split_vec.at(12);
					tmp_s_result.score_info = split_vec.at(13);
					r.insert(pair<string, solar_result>(tmp_s_result.q_id, tmp_s_result));
					q_id.insert(tmp_s_result.q_id);
			}// end else
		}//end for
		Fin.close();
	}//end try
	catch(char* pMsg) { cerr << endl << "Exception:" << pMsg << endl;}
}
