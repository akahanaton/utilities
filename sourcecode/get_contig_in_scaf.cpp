#include "../myHeader/general.h"

using namespace std;
struct contig_info
{
	string contig_name;
	int pos;
	string scaf_rcom;
	string blat_rcom;
	int len;
	int align_len;
	string q_start;
	string q_end;
};

struct less_pos
{
	inline  bool operator() ( const contig_info& c1, const contig_info& c2) const
	{ return ( c1.pos < c2.pos );}
};


typedef map<string, contig_info> contigs;
void ReadScafList(char* filename, map<string, string>& c_s_n, map<string, contigs >& s_info);
void ReadBlatResult(char* filename, map<string, string>& c_s_n, map<string, contigs >& s_info);

int main(int argc, char * argv[])
{
	map<string, string> contig_scaf_map;
	map<string, contigs> scaf_info;

	ReadScafList(argv[1], contig_scaf_map, scaf_info);
	cout << contig_scaf_map.size() << "\t" << scaf_info.size() << endl;
	ReadBlatResult(argv[2], contig_scaf_map, scaf_info);

	return 0;
}

void ReadBlatResult(char* filename, map<string, string>& c_s_n, map<string, contigs >& s_info)
{
	ifstream Fin(filename);
	string buf, cur_contig_name, cur_bac_name;
	string contig_in_scaf_name;

	string regstr = "(\\d+)\\s+(.)\\s+(.+)\\s+(.+)\\s+(\\d+)\\s+(\\d+)";
	boost::regex exp(regstr);
	boost::smatch what;

	set<string> bac_name;
	multimap<string, string> bac_scaf;

	for(;;)
	{
		getline(Fin, buf);
		if(Fin.eof())
			break;
		else
		{
			if(boost::regex_match(buf, what, exp))
			{
				string msg_1(what[1].first, what[1].second);
				string msg_2(what[2].first, what[2].second);
				string msg_3(what[3].first, what[3].second);
				string msg_4(what[4].first, what[4].second);
				string msg_5(what[5].first, what[5].second);
				string msg_6(what[6].first, what[6].second);

				//--------------------------------------------------
				// cout << msg_1 << " " << msg_2 << " " << msg_3 << " " << msg_4 << " " << msg_5 << " "<< msg_6 << endl;
				//--------------------------------------------------

				int align_len = boost::lexical_cast<int>(msg_1);
				if( align_len > 0)
				{
					contig_in_scaf_name =  c_s_n[msg_3];
					s_info[contig_in_scaf_name][msg_3].align_len = align_len;
					s_info[contig_in_scaf_name][msg_3].blat_rcom = msg_2;
					s_info[contig_in_scaf_name][msg_3].q_start = msg_5;
					s_info[contig_in_scaf_name][msg_3].q_end= msg_6;
					bac_name.insert(msg_4);
					bac_scaf.insert(pair<string, string>( msg_4, contig_in_scaf_name ));
				}
			}

		}
	} //end for

	pair<multimap<string, string>::iterator, multimap<string, string>::iterator> ret;
	multimap<string, string>::iterator mul_it;
	set<string>::iterator set_it;
	map<string, contig_info>::iterator contig_info_it;
	vector<contig_info> sort_vector;
	vector<contig_info>::iterator vec_it;
	for(set_it = bac_name.begin(); set_it != bac_name.end(); ++set_it)
	{
		cout << ">" << *set_it << endl;
		ret = bac_scaf.equal_range(*set_it);
		for (mul_it = ret.first; mul_it != ret.second; ++mul_it)
		{
			if (! mul_it->second.empty())
			{
			cout << "#\t" <<mul_it->second << endl;
			sort_vector.clear();
			sort_vector.reserve(s_info[mul_it->second].size());   // allocate space for 7 elements

			for (contig_info_it = s_info[mul_it->second].begin(); contig_info_it != s_info[mul_it->second].end(); ++ contig_info_it)
				sort_vector.push_back(contig_info_it->second);
			sort(sort_vector.begin(), sort_vector.end(), less_pos());
			for(vec_it = sort_vector.begin(); vec_it != sort_vector.end(); ++vec_it)
				cout << "\t*\t"<< vec_it->contig_name << "\t"  << vec_it->pos  << "\t"  << vec_it->len << "\t"  << vec_it->scaf_rcom
					 << "\t|\t" << vec_it->blat_rcom << "\t" << vec_it->align_len << "\t"  << vec_it->q_start << "\t"  << vec_it->q_end <<endl;
			}
		}
	}
}

void ReadScafList(char* filename, map<string, string>& c_s_n, map<string, contigs >& s_info)
{
    try
    {
		ifstream Fin(filename);
		if(!Fin) throw strcat("Cannot open input Fasta file ", filename);

		int newLine = 0;
		string		buf, cur_scaf_name, cur_contig_name;
		contigs 	tmp_contigs;
		contig_info tmp_contig_info;

		//5880690_2962_18.7       0       +       2962    2962
		string regstr = "(\\d+)\\s+([\\+\\-])\\s+(\\d+)";
		boost::regex exp(regstr);
		boost::smatch what;

		for(;;)
		{

			getline(Fin, buf);

			if(buf[0] == '>' ) 		newLine = 1;
			else if(Fin.eof()) 		newLine = 3;
			else if(buf.empty()) 	newLine = 4;
			else 					newLine = 2;

			if (newLine == 1)
			{
				int tmp_pos =  buf.find_first_of(" \t");
				cur_scaf_name = buf.substr(1, tmp_pos -1);
				tmp_contigs.clear();
				s_info.insert(pair<string, contigs>(cur_scaf_name,tmp_contigs));
			}
			if (newLine == 2)
			{
				int tmp_pos =  buf.find_first_of(" \t");
				cur_contig_name = buf.substr(0, tmp_pos);
				string tmp_str = buf.substr(tmp_pos + 1);
				int tmp_pos1 = tmp_str.find_first_of(" \t");
				int tmp_pos2 = tmp_str.find_last_of(" \t");
				tmp_contig_info.contig_name = cur_contig_name;
				tmp_contig_info.pos = boost::lexical_cast<int>(tmp_str.substr(0, tmp_pos1));
				tmp_contig_info.scaf_rcom = tmp_str.at(tmp_pos1 + 1);
				tmp_contig_info.len = boost::lexical_cast<int>(tmp_str.substr(tmp_pos2+1));
				c_s_n.insert(pair<string, string>(cur_contig_name,cur_scaf_name));
				s_info[cur_scaf_name].insert(pair<string, contig_info>( cur_contig_name, tmp_contig_info));
			}
			if(newLine == 3)
				break;
			if(newLine == 4)
				continue;

		} // end for
	}// end try
	catch (char* pMsg) { cerr << endl << "Exception:" << pMsg << endl; }
}



