#include <iostream>
using namespace std;

int Convert(int n, int Base);

int main ()
{
	int i;
	printf ("Enter a number: ");
	scanf ("%d",&i);
	Convert(i, 2);
	cout << endl;
	return 0;
}
int Convert(int n, int Base)
{
	if(n!=0)
	{
		Convert( n / Base, Base);
		cout << n%Base;
	}
}
