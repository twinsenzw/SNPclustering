#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>

using namespace std;


int main()
{
	ifstream prodigal_faa("block_fasta_protein.faa"); //the prodigal output
	ifstream original_fa("block_fasta.fa"); //original block sequences
	ofstream formatted_faa("block_fasta_protein_formatted.faa");

	map< string, vector<int> > pos_map; // fasta_header -> prodigal_pos -1 , original_pos !!!note that -1
	map< string, vector<int> >::iterator pmi;
	map< string, string> prodigal; // fasta_header -> sequence

	//read in original, build pos_map
	string line;
	for(;getline(original_fa,line);)
	{
		if(line=="") continue;

		string fasta_header;
		string sequence;
		if(line[0]=='>') 
		{
			fasta_header=line;
			getline(original_fa,sequence);
			
			int ori_i; //iterator of original sequence
			for(ori_i=0;ori_i<sequence.length();ori_i++)
			{
				if(sequence[ori_i]!='-')
				{
					pos_map[fasta_header].push_back(ori_i);
				}
			}
		}
	}

	//read in prodigal, format header
	
	for(;getline(prodigal_faa,line);)
	{	

		if(line=="") continue;

		string fasta_header;
		string sequence;
		string start;
		string end;
		if(line[0]=='>') 
		{
			istringstream linestream(line);
			
			//deal with header
			getline(linestream,fasta_header,'#');
			istringstream headerstream(fasta_header);
			string part; vector<string> parts;
			for(;getline(headerstream,part,'_');) parts.push_back(part);
			string reassembly="";
			reassembly+=parts[0];
			for(int i=1;i<parts.size()-1;++i)
			{
				reassembly+="_";
				reassembly+=parts[i];
			}
			
			
			//deal with positions
			getline(linestream,start,'#');
			getline(linestream,end,'#');
			int start_base=stoi(start);
			int end_base=stoi(end);

			
			pmi=pos_map.find(reassembly);
			int ori_start_base=-1;
			if(pmi!=pos_map.end())
			{
				ori_start_base=pmi->second[start_base-1];
			}
			int ori_end_base=-1;
			if(pmi!=pos_map.end())
			{
				ori_end_base=pmi->second[end_base-1];
			}
			formatted_faa<<">"<<reassembly<<"|"<<ori_start_base<<"|"<<ori_end_base<<endl;
		}
		else
		{
			formatted_faa<<line<<endl;
		}
	}
}			
	
				
		
	
