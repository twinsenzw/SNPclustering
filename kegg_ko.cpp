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

class feature
{
public:
	int start;
	int end;
	string feature;
};

class feature_count
{
public:
	int count;
	string feature;
};
int main(int argc, char *argv[])
{
	string argv1=argv[1];
	string argv2=argv[2];
	int threshold=stod(argv1);
	//parse ko_genes mapping file
	cout<<"*parsing ko_genes mapping file..."<<endl;
	ifstream ko_genes(argv2.c_str());
	map < string, string> gene_ko;
	map < string, string>::iterator gki;
	string line;

	for(;getline(ko_genes,line);)
	{
		
		string ko;
		string temp;
		string gene;
		istringstream linestream(line);
		linestream >> temp >> gene;
		istringstream tempstream(temp);
		getline(tempstream, ko, ':');
		getline(tempstream, ko, ':');
		
		gene_ko[gene]=ko;
	}

	//parse blast output
	ifstream blast_out("kegg_blast_output.txt");
	cout<<"*parsing blast output..."<<endl;
	ofstream output("ko_mapping_output.txt");

	for(;getline(blast_out,line);)
	{
		istringstream linestream(line);
		string header;
		string gene;
		string temp;
		string id;
		
		getline(linestream,header,'\t');
		getline(linestream,gene,'\t');
		getline(linestream,temp,'\t');
		getline(linestream,id,'\t');

		double did=stod(id);
		if(did>=threshold)
		{
			auto mi=gene_ko.find(gene);
			if(mi!=gene_ko.end()) output<<header<<"\t"<<gene<<"\t"<<mi->second<<endl;
			else output<<header<<"\t"<<gene<<"\t"<<endl;
		}
		
	}



	
}			
	
				
		
	
