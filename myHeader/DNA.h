#ifndef DNA_H_
#define DNA_H_

#include "general.h"
//#include <algrithom>

using namespace std;


// handling DNA molecules data
class DNA;

// similar to class Chain, but handling a single DNA chain
class fastaDNA;

// qual infomation
typedef string SEQID;
typedef vector<int> QUAL;
typedef vector<int> DEPTH;

/// class fastaDNA
class fastaDNA
{
public:
  string id;
  string description;
  string seq;

  fastaDNA(); // constructor
  fastaDNA(const fastaDNA& other_fastadna); // copy constructor
  void initial(); // clear the content of object
  void cleanSeq(); // clean the sequence, remove non "AGCT" charactors
  void set_id(string input_id);
  void set_description(string input_description);
  void set_seq(string input_seq);
  void reverse_compliment(); // return the reverse compliment
  fastaDNA assemble(fastaDNA c1, fastaDNA c2, int minmatch);
  int length();
  bool longerThan(int input_len);
  bool shorterThan(int input_len);

};

/// class CDS
class CDS
{
public:
/// /gene /locus_tag /protein_id /db_xref
  string gene;
  string locus_tag;
  string protein_id;
  string db_xref;

  string id;
  string cds_info;
  int  codon_start;
  pair<int, int> indices;
  bool recomplement;

  CDS();  // constructor
  CDS(const CDS& other_cds); // copy constructor
  void initial();
  void set_indices(int input_pos1,int input_pos2);

};


class DNA
{
public:
    vector<fastaDNA> 	chain;
    vector<CDS>	   		cds;
    map<string,int> 	oligo_counter;
    map<SEQID,QUAL> 	qual;
    map<SEQID,DEPTH> 	depth;

    DNA(); // constructor
    DNA(const DNA& other_dna); // copy constructor
    void readFasta(char *file); // read Fasta format sequence
    void writeFasta(char *file); // write Fasta format sequence
    void addChain( const fastaDNA& dc); // add a new chain member
    // another way to add a new chain member
    void addChain(const string& id,const string& description,const string& seq);
    void addCDS(CDS input_cds);
    void readGenBank(char *file);// read GenBank format sequence

    void readQual(char* file);
    void writeQual(char *file);

    void trimBeginEndX(int short_len);

    void getOligoNum();
    void oligoCount(const int& chain_index, const int& oligo_num);

    int  getChainNum();
    int  getCDSNum();

    fastaDNA getChain(const int& chain_index);
    fastaDNA getCDSSeq(const int& chain_index, const int& cds_index);

    // get the sequence have the input oligo
    bool getOligo2Id(string input_oligo, vector<fastaDNA>::iterator input_it);

	void Read454AlignDepth(char* file);
	void Write454AlignDepth(char* file);

	void UniqByDescription();

	void clear();
};

void getRandSeq(int input_len,int input_coverage, DNA &input_dna,DNA *out_dna);
typedef map<string, string> map_DNA;
void readFastaMap(char *file, map_DNA& out_seq_map); // read Fasta format sequence
void readFastqMap(char *file, map_DNA& out_seq_map); // read Fastq format sequence
void writeFastaMap(char *file, map_DNA& out_seq_map); // write Fasta format sequence

struct less_by_description
{
	inline bool operator() (const fastaDNA& f1, const fastaDNA& f2) const
	{
		return f1.description < f2.description;
	}
};

struct equal_by_description
{
	inline bool operator() (const fastaDNA& f1, const fastaDNA& f2) const
	{
		return f1.description ==  f2.description;
	}
};
#endif /*DNA_H_*/
