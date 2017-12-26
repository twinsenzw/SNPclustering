#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <ctime>

using namespace std;

class clusters
{
public:
	map<string, int> A_counts;
	map<string, int> B_counts;
	map<string, int> C_counts;

};
	
int main()
{

	ifstream ko_counts("cluster_ko_counts.txt");
	ifstream ko1("/home/zhouw/Kegg_db/ko/ko00001.keg");
	

	//generate maps of K number to A,B,C level features
	cout<<"*generating hierarchy..."<<endl;
	map<string, set<string>> B_map; //mapping K number to B level brite
	map<string, set<string>> A_map;
	map<string, set<string>> C_map;

	string line;

	string present_A;
	string present_B;
	string present_C;

	for(;getline(ko1,line);)
	{
		if(line=="") continue;
		if(line[0]=='#' || line[0]=='!' || line[0]=='+') continue;

		istringstream linestream(line);
		string temp;

		if(line[0]=='A') //new A level feature
		{
			getline(linestream,temp,'>');
			getline(linestream,temp,'<');
			present_A=temp;
		}
		else if(line[0]=='B') //new B level feature
		{
			getline(linestream,temp,'>');
			getline(linestream,temp,'<');
			present_B=temp;
		}
		else if(line[0]=='C') //new C level feature
		{
			present_C=line.substr(5,line.length()-5);
		}		
		else if(line[0]=='D') //K number found
		{
			string K;
			linestream >> temp >> K;
			
			if(A_map.find(K)!=A_map.end()) A_map[K].insert(present_A);
			else
			{
				set<string> new_set;
				new_set.insert(present_A);
				A_map[K]=new_set;
			}

			if(B_map.find(K)!=B_map.end()) B_map[K].insert(present_B);
			else
			{
				set<string> new_set;
				new_set.insert(present_B);
				B_map[K]=new_set;
			}
			
			if(C_map.find(K)!=C_map.end()) C_map[K].insert(present_C);
			else
			{
				set<string> new_set;
				new_set.insert(present_C);
				C_map[K]=new_set;
			}			
			
		}
	}

	//reading and counting cluster groups
	cout<<"*mapping Kegg orthologs..."<<endl;
	int present_cluster=-1;
	vector<clusters> all_clusters;
	for(;getline(ko_counts,line);)
	{
		istringstream linestream(line);
		int cluster;
		string K;
		int count;

		linestream >> cluster >> K >> count;
		/*if(cluster==0)
		{
			for(map<string, int>::iterator mi=new_cluster.A_counts.begin();mi!=new_cluster.A_counts.end();mi++)
			{
				cout<<mi->first<<"."<<mi->second<<"\t";
			}
			cout<<endl;
		}*/

		if(cluster==present_cluster) //adding counts
		{
			if(A_map.find(K)!=A_map.end()) //if the K number can be mapped to some A feature
			{
				set<string>::iterator si;
				for(si=A_map[K].begin();si!=A_map[K].end();si++) //*si is the A feature string
				{
					if(all_clusters.back().A_counts.find(*si)!=all_clusters.back().A_counts.end()) //the A feature already recorded
					{
						all_clusters.back().A_counts[*si]=all_clusters.back().A_counts[*si]+count;
					}
					else //a new A feature
					{
						all_clusters.back().A_counts[*si]=count;
					}
				}
			}
	
			if(B_map.find(K)!=B_map.end()) //if the K number can be mapped to some B feature
			{
				set<string>::iterator si;
				for(si=B_map[K].begin();si!=B_map[K].end();si++) //*si is the B feature string
				{
					if(all_clusters.back().B_counts.find(*si)!=all_clusters.back().B_counts.end()) //the B feature already recorded
					{
						all_clusters.back().B_counts[*si]=all_clusters.back().B_counts[*si]+count;
					}
					else //a new B feature
					{
						all_clusters.back().B_counts[*si]=count;
					}
				}
			}	

			if(C_map.find(K)!=C_map.end()) //if the K number can be mapped to some C feature
			{
				set<string>::iterator si;
				for(si=C_map[K].begin();si!=C_map[K].end();si++) //*si is the C feature string
				{
					if(all_clusters.back().C_counts.find(*si)!=all_clusters.back().C_counts.end()) //the C feature already recorded
					{
						all_clusters.back().C_counts[*si]=all_clusters.back().C_counts[*si]+count;
					}
					else //a new c feature
					{
						all_clusters.back().C_counts[*si]=count;
					}
				}
			}
		}
		else //a new cluster is met
		{

			clusters new_cluster;

			if(A_map.find(K)!=A_map.end()) //if the K number can be mapped to some A feature
			{
				set<string>::iterator si;
				for(si=A_map[K].begin();si!=A_map[K].end();si++) //*si is the A feature string
				{
					if(new_cluster.A_counts.find(*si)!=new_cluster.A_counts.end()) //the A feature already recorded
					{
						new_cluster.A_counts[*si]=new_cluster.A_counts[*si]+count;
					}
					else //a new A feature
					{
						new_cluster.A_counts[*si]=count;
					}
				}
			}
	
			if(B_map.find(K)!=B_map.end()) //if the K number can be mapped to some B feature
			{
				set<string>::iterator si;
				for(si=B_map[K].begin();si!=B_map[K].end();si++) //*si is the B feature string
				{
					if(new_cluster.B_counts.find(*si)!=new_cluster.B_counts.end()) //the B feature already recorded
					{
						new_cluster.B_counts[*si]=new_cluster.B_counts[*si]+count;
					}
					else //a new B feature
					{
						new_cluster.B_counts[*si]=count;
					}
				}
			}	

			if(C_map.find(K)!=C_map.end()) //if the K number can be mapped to some C feature
			{
				set<string>::iterator si;
				for(si=C_map[K].begin();si!=C_map[K].end();si++) //*si is the C feature string
				{
					if(new_cluster.C_counts.find(*si)!=new_cluster.C_counts.end()) //the C feature already recorded
					{
						new_cluster.C_counts[*si]=new_cluster.C_counts[*si]+count;
					}
					else //a new c feature
					{
						new_cluster.C_counts[*si]=count;
					}
				}
			}			
			present_cluster=cluster;
			all_clusters.push_back(new_cluster);
		}
	}
	
	//OUTPUT
	ofstream Aout("Aout.txt");
	ofstream Bout("Bout.txt");
	ofstream Cout("Cout.txt");

	for(int i=0;i<all_clusters.size();i++)
	{
		map<string,int>::iterator mi;
		int total_counts=0;
		for(mi=all_clusters[i].A_counts.begin();mi!=all_clusters[i].A_counts.end();mi++)
		{
			total_counts+=mi->second;
		}
		for(mi=all_clusters[i].A_counts.begin();mi!=all_clusters[i].A_counts.end();mi++)
		{
			Aout<<i<<"\t"<<mi->first<<"\t"<<(double) mi->second/total_counts<<endl;
		}
		
		total_counts=0;
		for(mi=all_clusters[i].B_counts.begin();mi!=all_clusters[i].B_counts.end();mi++)
		{
			total_counts+=mi->second;
		}
		for(mi=all_clusters[i].B_counts.begin();mi!=all_clusters[i].B_counts.end();mi++)
		{
			Bout<<i<<"\t"<<mi->first<<"\t"<<(double) mi->second/total_counts<<endl;
		}
		
		total_counts=0;
		for(mi=all_clusters[i].C_counts.begin();mi!=all_clusters[i].C_counts.end();mi++)
		{
			total_counts+=mi->second;
		}
		for(mi=all_clusters[i].C_counts.begin();mi!=all_clusters[i].C_counts.end();mi++)
		{
			Cout<<i<<"\t"<<mi->first<<"\t"<<(double) mi->second/total_counts<<endl;
		}
	}
}						
