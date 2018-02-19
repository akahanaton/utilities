#include "/share/raid5/wenming/myProgram/myHeader/DNA.h"
#include "/share/raid5/wenming/myProgram/myHeader/blastparser.h"

using namespace std;
void usage();

int main(int argc, char *argv[])
{
if(argc == 1)
{
	usage();
}
else
 {
      int opt;
      int align_len_cutoff = 0;
      int identity_cutoff = 0;
      char *in_seq = 0,  *blast_file = 0 , *out_seq = 0;
  while((opt = getopt(argc, argv, "f:o:b:i:l:a:r")) != -1 )
  {
     switch(opt)
      {
      case 'f':
			in_seq = optarg;
			break;
      case 'b':
			blast_file = optarg;
			break;
      case 'o':
			out_seq = optarg;
			break;
      case 'i':
			identity_cutoff = boost::lexical_cast<int>(optarg);
			break;
      case 'l':
			align_len_cutoff = boost::lexical_cast<int>(optarg);
			break;
      case '?':
			printf("unknown option: %c\n", optopt);
			return EXIT_SUCCESS;
      }
   } // end while

	set<string> erase_id;

	map_DNA contigs;
	readFastaMap(in_seq, contigs);

	BlastParser bp;
	bp.ReadBlastResult(blast_file, 100);
   if(identity_cutoff > 0)
		bp.RmLowIdentity( identity_cutoff );
   if(align_len_cutoff > 0)
	   bp.RmShortAlign(align_len_cutoff);

	sort(bp.BlastResults.begin(), bp.BlastResults.end(), less_by_QueryID_QStart());

	//--------------------------------------------------
	// bp.GetResult(false);
	//--------------------------------------------------

   map<string, int> appended_tag;
   map<string, int> removale_tag;
   map<string, int> to_be_appended_tag;
   map<string, string>::iterator m_it = contigs.begin();
   for(; m_it != contigs.end(); ++m_it)
	   removale_tag.insert( map<string, int>::value_type(m_it->first, 0) );

   vector<BLAST8RESULT>::iterator cur_it = bp.BlastResults.begin(), next_it = bp.BlastResults.begin() + 1;
   for ( ;next_it!= bp.BlastResults.end(); ++cur_it, ++next_it)
   {
		if( cur_it->SubjectID != next_it->SubjectID)
		{
			if(cur_it->QEnd > next_it->QStart && cur_it->SEnd == cur_it->SLen && next_it->SStart == 1 )
			{
				to_be_appended_tag[next_it->SubjectID] += 1;
			}
		}
   }
   cur_it = bp.BlastResults.begin();
   next_it = bp.BlastResults.begin() + 1;
   string append_str;
   bool append_begin = false;
   for ( ;next_it!= bp.BlastResults.end(); ++cur_it, ++next_it)
   {
		if( cur_it->SubjectID != next_it->SubjectID)
		{
			int overlap_len = cur_it->QEnd - next_it->QStart + 1;
			if( overlap_len > 1 && cur_it->SEnd == cur_it->SLen && next_it->SStart == 1 && removale_tag[next_it->SubjectID] == 0)
			{
				if(cur_it->RComplement != next_it->RComplement )
				{
				  string tmp_seq;
				  tmp_seq.resize(contigs[next_it->SubjectID].length());
				  reverse_copy(contigs[next_it->SubjectID].begin(), contigs[next_it->SubjectID].end(), tmp_seq.begin());
				  transform(tmp_seq.begin(), tmp_seq.end(), contigs[next_it->SubjectID].begin(), compliment_nucleotide);
				}
				string app_str = contigs[next_it->SubjectID].substr(overlap_len);
				cout << cur_it->SubjectID<< "\t" << contigs[cur_it->SubjectID].length() <<"\t"<< overlap_len << "\t" << endl;
				if(append_begin)
					append_str.append(app_str);
				else
				{
					append_str = contigs[cur_it->SubjectID];
					append_str.append(app_str);
					append_begin = true;
				}
				next_it->RComplement = cur_it->RComplement;
				if (removale_tag[cur_it->SubjectID] != 2)
					removale_tag[cur_it->SubjectID] = 1;
			}
			else
			{
				if(append_begin)
				{
				contigs[cur_it->SubjectID] = append_str;
				removale_tag[cur_it->SubjectID] = 2;
				append_begin = false;
				append_str.clear();
				}
			}
		}
   }//end for

   for(m_it = contigs.begin(); m_it != contigs.end(); ++m_it)
		if(removale_tag[m_it->first] == 1)
		{
			cout << "\terase id\t" << m_it->first << endl;
			contigs.erase(m_it->first);
		}

   writeFastaMap(out_seq, contigs);
}

  return EXIT_SUCCESS;
}

void usage()
{
  cout << "linker two contigs which have overlap according to blast result" << endl;
  cout << "Options:" << endl;
  cout << "\t\t-f STR: read in sequences file"<<endl;
  cout << "\t\t-o STR: write out sequences file"<<endl;
  cout << "\t\t-b STR: blast result file"<<endl;
  cout << "\t\t-i INT: Identity, 0"<<endl;
  cout << "\t\t-l INT: alingment length cutoff, default 0"<<endl;
}
