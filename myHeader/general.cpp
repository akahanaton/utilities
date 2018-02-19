#include "general.h"
#include <string>


int my_lowercase(char *s) {
  while(*s) {
	*s=tolower(*s);
	s++;
  }
  return 0;
}

int my_chopSpace(char *s) {
  for(;*s;s++){
	if(*s==' ') {
		char *s2=s;
		for(;*s2;s2++) {
			*(s2)=*(s2+1);
		}
	}
  }
  return 0;
}

char compliment_nucleotide(char ch) {
  std::string str = "ACGT";
  std::string ts(1,ch);
  if( ch == 'A' ) return 'T';
  else if( ch == 'T' ) return 'A';
  else if( ch == 'C' ) return 'G';
  else return 'C';
}

int* getRandNum(int upper_bound,int numbers)
{
    /// const size_t N=25223,M=170000;
    int* ia=new int[upper_bound];
    for(size_t i=0;i!=upper_bound;++i)
         *(ia+i)=i+1;
    srand((unsigned)time(NULL));
    for (size_t i=0;i!=numbers;++i)
    {
      myswap(*(ia+i),*(ia+i+rand()*rand()*rand()%( upper_bound - i ) ) );
       //swap(*(ia+i),*(ia+rand()*rand()*rand()%M));
       //cout << *(ia+i) << endl;
    }
    return ia;
}
