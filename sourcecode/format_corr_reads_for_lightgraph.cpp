#include "../myHeader/DNA.h"
using namespace std;

int main( int argc, char* argv[] )
{
	DNA  dna1, dna2, out_seq;
	dna1.readFasta(argv[1]);
	dna2.readFasta(argv[2]);
	int dna1_num = dna1.chain.size();
	int dna2_num = dna2.chain.size();

    map<string, string> ref1, ref2;
	map<string, string>::iterator m_it;
	string tmp_id, dna_id, dna_des;
	int pos;

	cout << "test 1" << endl;

	for(int i = 0; i<dna1_num; ++i)
	{
		pos = dna1.chain.at(i).id.find_last_of("/");
		tmp_id = dna1.chain.at(i).id.substr(0,pos);
		ref1.insert( map<string,string>::value_type(tmp_id, dna1.chain.at(i).seq) );
	}
	dna1.chain.clear();

	for(int i = 0; i<dna2_num; ++i)
	{
		pos = dna2.chain.at(i).id.find_last_of("/");
		tmp_id = dna2.chain.at(i).id.substr(0,pos);
		m_it = ref1.find(tmp_id);
		if(m_it != ref1.end())
		{
			dna_id = tmp_id;
			dna_id.append("/1");
			dna_des = dna_id.append(" ");
			dna_des.append(argv[4]);
			out_seq.addChain(dna_id, dna_des, ref1[tmp_id]);
			ref1.erase(m_it);

			dna_id = tmp_id;
			dna_id.append("/2");
			dna_des = dna_id.append(" ");
			dna_des.append(argv[4]);
			out_seq.addChain(dna_id, dna_des, dna2.chain.at(i).seq);
		}
		else
			ref2.insert( map<string,string>::value_type( tmp_id, dna2.chain.at(i).seq) );
	}
	dna2.chain.clear();
	cout << "test 2" << endl;

	if(!ref1.empty())
		for(m_it = ref1.begin(); m_it != ref1.end(); ++m_it )
		{
			dna_id = m_it->first;
			dna_id.append("/1");
			out_seq.addChain(dna_id, dna_id, m_it->second);
		}
	if(!ref2.empty())
		for(m_it = ref2.begin(); m_it != ref2.end(); ++m_it )
		{
			dna_id = m_it->first;
			dna_id.append("/2");
			out_seq.addChain(dna_id, dna_id, m_it->second);
		}

	cout << "test 3" << endl;
	out_seq.writeFasta(argv[3]);
	return 0;
}
