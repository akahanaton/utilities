#include "../myHeader/DNA.h"
using namespace std;

int main( int argc, char* argv[] )
{
	DNA mydna,gene;
	fastaDNA tmpDNA;
	string buf;
	mydna.readFasta(argv[1]);
	ifstream Fin( argv[2]);
    map<string, string> ref;
//	for(int i = 0; i< mydna.chain.size();++i)
//	{
//		ref.insert( map<string,string>::value_type(mydna.chain.at(i).id, mydna.chain.at(i).seq) );
//	}
	int s,e;
	int gene_index  = 1 ;
	string gene_index_str;
	for(;;)
	{
		getline(Fin, buf);
		if(buf[0] == '>')
			continue;
		else if (Fin.eof())
			break;
		else
		{
		//vector<string> vs;
		//boost::split(vs, buf, boost::is_any_of(" \t"));
		boost::regex des("(\\w+\\d+)\\s+(\\d+)\\s+(\\d+)");
		boost::smatch what;
		string::const_iterator start = buf.begin();
		string::const_iterator end = buf.end();
		if(boost::regex_search(start, end, what, des))
		{
			string tmp_s(what[2].first, what[2].second);
			cout << tmp_s<<endl;
			string tmp_e(what[3].first, what[3].second);
			cout << tmp_e<<endl;
		 s = boost::lexical_cast<int>(tmp_s) -1;
		 e = boost::lexical_cast<int>(tmp_e) -1;
		}
		//copy(vs.begin(), vs.end(), ostream_iterator<string >(cout, " "));
		if (s > e )
//		 gene.addChain( buf, buf, ref[vs.at(0)].substr(e,s-e+1 ) );
		 gene.addChain( buf, buf, mydna.chain.at(0).seq.substr(e,s-e+1) );
		else
//		 gene.addChain( buf, buf, ref[vs.at(0)].substr(s,e-s+1 ) );
		 gene.addChain( buf, buf, mydna.chain.at(0).seq.substr(s,e-s+1) );
		cout << endl;
		}
	}
	gene.writeFasta(argv[3]);
	return 0;
}
