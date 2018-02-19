#include<stdio.h>
#include<iostream>
#include<string>
#include<unistd.h>
#include<sys/stat.h>
using namespace std;
int main (int argc, char *argv[])
{
	string name;
	name = strcat( const_cast<char*>(name.c_str()), argv[1]);
	name.append("/");
	cout <<"argv[1] " << name << endl;
    if(mkdir(argv[1], S_IRWXU | S_IRWXG | S_IRWXO) != 0)
			cout << "maked folder " << *argv[1]<< endl;
	return EXIT_SUCCESS;
}
