#include "../myHeader/DNA.h"
using namespace std;

int main( int argc, char* argv[] )
{
	DNA mydna,gene;
	string buf;
	mydna.readFasta(argv[1]);
	ifstream Fin( argv[2]);
    map<string, string> ref;
	for(int i = 0; i< mydna.chain.size();++i)
	{
		ref.insert( map<string,string>::value_type(mydna.chain.at(i).id, mydna.chain.at(i).seq) );
	}

	int gene_index  = 1 ;
	string last_gene_id, cur_gene_id;
	for(;;)
	{
		getline(Fin, buf);
		if (Fin.eof())
			break;
		else
		{
		string des;
		vector<string> vs;
		boost::split(vs,buf, boost::is_any_of(" \t"));
		copy(vs.begin(), vs.end(), ostream_iterator<string >(cout, " "));
		cur_gene_id = *vs.begin();
		int s = boost::lexical_cast<int>(vs.at(7) ) -1;
		int e = boost::lexical_cast<int>(vs.at(8) ) -1;
		des.append(cur_gene_id);
		des.append(" ");
		des.append(vs.at(5));
		des.append(" ");
		des.append( vs.at(7)  );
		des.append(" ");
		des.append(  vs.at(8) );
		if (s > e )
		 gene.addChain( vs.at(0), des, ref[vs.at(5)].substr(e,s-e+1 ) );
		else
		 gene.addChain( vs.at(0), des, ref[vs.at(5)].substr(s,e-s+1 ) );
		cout << endl;
		}
	}
	gene.writeFasta(argv[3]);
	return 0;
}
