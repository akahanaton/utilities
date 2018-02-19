#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "/genome/cpptest/BioC++.h"
using namespace std;

int main(int argc, char *argv[])
{
        DNA dna;
	if (argc != 4)
	{
	cout << "input 1: The repeat sequencs files , like : " << endl;
	cout << "\tAATAATGATAATAAT" << endl;
	cout << "\tATATATATATATATA" << endl;
	cout << "\tGGCGGCGGCGGCG" << endl;
	cout << "\tAATAATGATAATAAT" << endl;
	cout << "input 2: The sequence ID." << endl;
	cout << "input 3: The out Sequencs files, fasta format." << endl;
 	 return EXIT_SUCCESS;
	}

	ifstream Fin_rep(argv[1]);
	if(!Fin_rep) 
		cout << "Can't open repeats file: " << argv[2] << endl;
	string repeat,tmpid;
        
	int i = 0;	
	for(;;)
	{
	getline(Fin_rep, repeat);
	++i;
	if(Fin_rep.eof()) 
		break;
	else
	 {
	// regex repeat_tag(repeat);
	string index = lexical_cast<string>(i);
	tmpid = lexical_cast<string>(argv[2]);
	tmpid.append("-");
	tmpid.append(index);
	dna.addChain(tmpid,tmpid,repeat);
       	 }
	}// end for
	dna.writeFasta(argv[3]);
  return EXIT_SUCCESS;
}
