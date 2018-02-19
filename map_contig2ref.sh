rm map_contig2ref
g++ -c myHeader/DNA.cpp -Wno-deprecated
g++ -c myHeader/general.cpp -Wno-deprecated
g++ -c myHeader/blastparser.cpp -Wno-deprecated
g++ -v *.o -I ~/lib/include -L ~/lib/lib sourcecode/map_contig2ref.cpp -lboost_regex -o map_contig2ref
./map_contig2ref -b /share/raid5/wenming/cucumber/eva.1113/melon_est_unigene.1000.VS.scaffold_step3.seq.blat -i 90  -e | grep MU3451

