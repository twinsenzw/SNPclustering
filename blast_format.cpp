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
	int cluster_number=stoi(argv1);
	//parse ko_genes mapping file
	cout<<"*parsing ko_genes mapping file..."<<endl;
	ifstream ko_genes("/home/zhouw/Kegg_db/ko_genes.list");
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
	map <int, vector<feature>> block_feature; //int block_id -> feature vector
	map <int, vector<feature>>::iterator bfi;
	for(;getline(blast_out,line);)
	{
		string start_str;
		string end_str;
		string blockid_str;
		string gene_str;
		string temp;

		istringstream linestream(line);
		linestream >> temp >> gene_str;

		istringstream tempstream(temp);
		getline(tempstream,blockid_str,'|');
		getline(tempstream,blockid_str,'|');

		getline(tempstream,start_str,'|');
		getline(tempstream,start_str,'|');
		getline(tempstream,start_str,'|');

		getline(tempstream,end_str,'|');

		feature newfeature;
		newfeature.start=stoi(start_str);
		newfeature.end=stoi(end_str);
		newfeature.feature=gene_str;
		
		bfi=block_feature.find(stoi(blockid_str));
		if(bfi==block_feature.end())
		{
			vector<feature> features;
			features.push_back(newfeature);
			block_feature[stoi(blockid_str)]=features;
		}
		else
		{
			bfi->second.push_back(newfeature);
		}
	}


	//parse snp_cluster
	cout<<"*mapping feature to snp..."<<endl;
	ifstream snp_cluster("snp_clusters.txt");
	ofstream output("snp_feature.txt");


	for(;getline(snp_cluster,line);)
	{
		string cluster_str;
		string blockid_str;
		string position_str;
		istringstream linestream(line);
		linestream >> cluster_str >> blockid_str >> position_str;
		int position=stoi(position_str);
		int cluster=stoi(cluster_str);
		vector<string> ko;
		//map by block
		bfi=block_feature.find(stoi(blockid_str));
		if(bfi==block_feature.end())
		{
			ko.push_back("block_has_no_feature");
		}
		else
		{
			for(int i=0; i<bfi->second.size(); ++i)
			{
				if(position>=(bfi->second[i]).start && position<=(bfi->second[i]).end) //snp mapped to this feature position
				{
					gki=gene_ko.find((bfi->second[i]).feature); //if the feature is mapped to a ko
					if(gki!=gene_ko.end())
					{
						ko.push_back(gki->second);

					}
				}
			}
		}
		if(ko.size()!=0)
		{
			
			output<<cluster_str<<"\t"<<blockid_str<<"\t"<<position;
			for(int j=0;j<ko.size();j++)
			{
				output<<"\t"<<ko[j];
			}
			output<<endl;
		}
	}

	output.close();
	//ko count for cluster
	cout<<"*counting kos..."<<endl;
	ofstream output1("cluster_ko_counts.txt");
	ifstream snp_feature("snp_feature.txt");
	
	int present_cluster=-1;

	vector<feature_count> count;
	for(;getline(snp_feature,line);)
	{
		istringstream linestream(line);
		string field;

		string cluster_str;
		string block_str;
		string position_str;
		
		getline(linestream,cluster_str,'\t');
		getline(linestream,block_str,'\t');
		getline(linestream,position_str,'\t');

		int cluster=stoi(cluster_str);

		if(present_cluster!=cluster) //new cluster started
		{
			
			for(int i=0;i<count.size();i++)
			{
				output1<<present_cluster<<"\t"<<count[i].feature<<"\t"<<count[i].count<<endl;
			}
			present_cluster=cluster;
			count.clear();
		}
			
		for(;getline(linestream,field,'\t');)
		{
			int newfeature_bool=1;
			for(int i=0;i<count.size();i++)
			{
				if(count[i].feature==field)
				{	
					count[i].count++;
					
					newfeature_bool=0;
					break;
				}
			}
			if(newfeature_bool==1)
			{
				feature_count newfeaturecount;
				newfeaturecount.feature=field;
				newfeaturecount.count=1;
				count.push_back(newfeaturecount);
			}
		}
		
	}

	
}			
	
				
		
	
