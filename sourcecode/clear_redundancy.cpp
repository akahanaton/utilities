#include "../myHeader/general.h"
#include "../myHeader/DNA.h"
#include "../myHeader/blastparser.h"

using namespace std;

void makeset(const vector<int>& v1, vector<int>& v2);
int findset(const int& query, const vector<int>& p);
void mergetrees(const int& q1, const int& q2, vector<int>& p);
void unionset(const int& q1, const int& q2, vector<int>& p);

struct myclass {
	  inline bool operator() (const int& i,const int& j) const { return (i<j);}
} myobject;

struct myclass_2 {
	  inline bool operator() (const BLAST8RESULT& r) const
	  {
		    return (100 * r.AlignmentLength / r.QLen < 95);
	  }
} self_or_low;

int main(int argc, char* argv[])
{
	BlastParser bp;
	DNA mydna;
	int tmp_1, tmp_2;
	map<int, int>query_len;
	multimap<int, int> redundant_query;
	vector<BLAST8RESULT>::iterator blast_it;
	mydna.readFasta("/share/raid5/wenming/cucumber/est.no.redundancy/combine.index.est");
	bp.ReadBlastResult("/share/raid5/wenming/cucumber/est.no.redundancy/combine.est.blast",10000);
	//--------------------------------------------------
	// mydna.readFasta(argv[1]);
	// bp.ReadBlastResult(argv[2],10000);
	//--------------------------------------------------
	cout << "begin remove "<<bp.BlastResults.size()<<endl;
    bp.BlastResults.erase(
		remove_if( bp.BlastResults.begin(), bp.BlastResults.end(), self_or_low ),
		bp.BlastResults.end() );
	cout << "after remove "<<bp.BlastResults.size()<<endl;
	for(blast_it = bp.BlastResults.begin(); blast_it != bp.BlastResults.end(); ++blast_it)
	{
		tmp_1 = boost::lexical_cast<int>(blast_it->QueryID);
		tmp_2 = boost::lexical_cast<int>(blast_it->SubjectID);
		query_len.insert(pair<int,int>(tmp_1, blast_it->QLen ) );
		if( tmp_1 > tmp_2)
			myswap(tmp_1, tmp_2);
		redundant_query.insert(pair<int, int>(tmp_1, tmp_2) );
	}
	bp.BlastResults.clear();

	map<int, int>::iterator sm_it;
	multimap<int,int>::iterator mm_it;
	vector<int>::iterator v_it;

	int query_index = 0;
	for(sm_it = query_len.begin(); sm_it != query_len.end(); ++sm_it)
		if(sm_it->first > query_index)
			query_index = sm_it->first;

	vector<int> parent;
	parent.reserve(query_index + 1);
	for(int i = 0; i <= query_index; ++i)
		parent.push_back(i);

	for(mm_it = redundant_query.begin(); mm_it != redundant_query.end(); ++mm_it)
	{
		unionset(mm_it->first, mm_it->second, parent);
		//cout << mm_it->second << "\t" <<mm_it->first<< endl;
	}

	redundant_query.clear();

	map<int, int>clean_query;
	for(int i = 0; i<=query_index; ++i)
	{
		if( i > parent[i])
		{
			// first is query, second is parent
			cerr<<i<<"\t"<<parent[i]<<endl;
			clean_query.insert(pair<int, int>( i,parent[i]));
		}
	}
	cout << "total repeat sequence\t" <<clean_query.size() << endl;

	vector<int> removeable_query;
	vector<int> same_cluster;
	for(int i = 0; i<=query_index; ++i)
	{
		same_cluster.clear();
		same_cluster.push_back(i);
		for(sm_it = clean_query.begin(); sm_it != clean_query.end(); ++sm_it)
		{
			if(sm_it->second == i) // have same parent, belong to the same cluster
				same_cluster.push_back(sm_it->first);
		}
		if(same_cluster.size() > 1)
		{
			vector<int>::iterator saved_query_it;
			int saved_query_len = 0;
			for(v_it = same_cluster.begin(); v_it != same_cluster.end(); ++v_it)
				if(query_len[*v_it] > saved_query_len)
					{
						saved_query_it = v_it;
						saved_query_len = query_len[*v_it];
					}
			same_cluster.erase(saved_query_it);
			removeable_query.insert(removeable_query.begin(),same_cluster.begin(),same_cluster.begin()+same_cluster.size());
		}
	}

	sort(removeable_query.begin(), removeable_query.end(), myobject);
	v_it = unique(removeable_query.begin(), removeable_query.end());
	removeable_query.resize( v_it - removeable_query.begin() );

	vector<int>::reverse_iterator rev_it;
	for(rev_it = removeable_query.rbegin(); rev_it != removeable_query.rend(); ++rev_it)
	{
		mydna.chain.erase(mydna.chain.begin() + *rev_it);
	}

	cout << "remove sequence\t" << removeable_query.size() << endl;
	cout <<"last no repeat sequence\t" << mydna.chain.size()<<endl;
	mydna.writeFasta("test.fa");
	return 0;
}

void makeset(const vector<int>& v1, vector<int>& v2)
{
	int vec_num = v1.size();
	v2.clear();
	v2.reserve(vec_num);
	for(int i=0; i < vec_num; ++i)
		v2.push_back(v1[i]);
}

int findset(const int& query, const vector<int>& p)
{
	int i = query;
	while( i!= p[i])
		i = p[i];
	return i;
}

void mergetrees(const int& q1, const int& q2, vector<int>& p)
{
	//set the smaller the root
	if(q1 > q2)
		p[q1] = q2;
	else
		p[q2] = q1;
}

void unionset(const int& q1, const int& q2, vector<int>& p)
{
	mergetrees(findset(q1, p), findset(q2, p), p);
}
