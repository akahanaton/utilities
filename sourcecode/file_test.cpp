#include <iostream>
#include <fstream>
#include <list>
#include <iterator>
#include "../myHeader/general.h"
using namespace std;

//ifstream::pos_type size;
char * memblock;
int size, line_gegin, line_end;

int main ()
{
	clock_t start, finish;
	double totaltime;
	start = clock();

	ifstream file ("ysp.soap.p.50200", ios::in|ios::binary|ios::ate);
	if (file.is_open())
	{
		size = file.tellg();
		memblock = new char [size];
		file.seekg (0, ios::beg);
		file.read (memblock, size);
		finish = clock();
		totaltime = (double)(finish - start)/CLOCKS_PER_SEC;
		cout << "program takes  " << totaltime << "  seconds"<<endl;
		//string str(memblock, 10000);
		//cout << str.length() << endl;
		cout << memblock[size -2] << "memblock[szie -2 ]" << endl;
		cout << memblock[size -1] << "memblock[szie -1 ]" << endl;
		delete[] memblock;
		finish = clock();
		totaltime = (double)(finish - start)/CLOCKS_PER_SEC;
		cout << "program takes  " << totaltime << "  seconds"<<endl;
		file.close();
	}
	else cout << "Unable to open file"<< endl;

	finish = clock();
	totaltime = (double)(finish - start)/CLOCKS_PER_SEC;
	cout << "program takes  " << totaltime << "  seconds"<<endl;

	return 0;
}

