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

typedef vector <vector<int>> d_vec;
typedef vector <d_vec> t_vec;
typedef map <t_vec, double> t_map;

class BLOCK
{
public:
	vector<string> strain;
	vector<string> sequence;
};	

class SNP
{
public:
	int block;
	int pos; //0-based

	//**The sequence in the vector is based on a predefined map of strain names to index.
	
	map<string, int> src; //strain -> nttype
};

class CLUSTER
{
public:
	vector<SNP> snp;
	int block;
	SNP center;
};

int nt_type(char nt)
{
	int type;
	if(nt=='A'||nt=='a') type=0;
	else if(nt=='T'||nt=='t') type=1;
	else if(nt=='G'||nt=='g') type=2;
	else if(nt=='C'||nt=='c') type=3;
	else type=4;

return type;
}

double calc_dis(vector<int> &v1, vector<int> &v2)
{
	double dis=0.0;
	for(int i=0;i<v1.size();i++)
	{
		dis=dis+(double)(v1[i]-v2[i])*(v1[i]-v2[i]);
	}
return dis;
}

double min_dis(d_vec &v1, d_vec &v2) //brute force
{
	double dis=-1; //return value
	int mem[5][5];
	for(int k=0;k<5;k++)
	{
		for(int l=0;l<5;l++)
		{
			mem[k][l]=-1;
		}
	}


	int p[]={0,1,2,3,4}; sort(p,p+5);
	do
	{
		double present_dis=0.0;
		for(int i=0;i<5;i++)
		{
			
				if(mem[i][p[i]]==-1) mem[i][p[i]]=calc_dis(v1[i],v2[p[i]]);
				present_dis=present_dis+mem[i][p[i]];

		}
		if(dis==-1) dis=present_dis;
		else if(dis>present_dis) dis=present_dis;
  	} while(next_permutation(p,p+5));

return dis;
}


/*

double min_dis_DP(t_vec &v, map<t_vec,double> mem) ; //DP, f declare. t_vec contains two snp_vector results of two snp sites.

double min_dis_DP(t_vec &v, map<t_vec,double> mem) 
{	
	double dis=-1;
	
	//cout<<"Calculating "<<v[0].size()<<","<<v[1].size()<<" nts..."<<endl;
	
	if(mem.find(v)!=mem.end()) dis=(mem.find(v)->second);
	else //not in memory
	{
		

		if(v[0].size()==1)  //initial condition
		{
			dis=calc_dis(v[0][0],v[1][0]); 
			t_vec vsub;
			vsub.push_back(v[0]);
			vsub.push_back(v[1]);
			mem.insert(make_pair(vsub, dis)); 
			goto end;
		}
			
		for(int i=0;i<v[0].size();i++)
		{
			d_vec v1sub=v[0];
			v1sub.erase(v1sub.begin()+i);
			
			for(int j=0;j<v[1].size();j++)
			{
				d_vec v2sub=v[1];
				v2sub.erase(v2sub.begin()+j);
				t_vec vsub;
				vsub.push_back(v1sub);
				vsub.push_back(v2sub);

				double pdis=calc_dis(v[0][i],v[1][j])+min_dis_DP(vsub,mem);
				if(dis==-1||dis>pdis) dis=pdis;
				mem.insert(make_pair(vsub, dis));
				//if(v[0].size()==3) cout<<pdis<<"\t"<<calc_dis(v[0][i],v[1][j])<<endl;
			}

		}

		
	}

end:
return dis;

}
*/

void maf_parse(ifstream &maf, vector<BLOCK> &block, int total_strain) //only look at blocks with total_strain number of sequences
{
	string line;
	for(;getline(maf,line);)
	{
		if(line[0]=='#') continue;
		if(line[0]=='a')
		{
			istringstream linestream(line);
			string field;
			for(;getline(linestream,field,'=');) ;
			if(stoi(field)==total_strain)
			{
				BLOCK newblock;
				for(int i=0;i<total_strain;i++)
				{
					getline(maf,line);

					istringstream blockstream(line);
					string strain_contig; string sequence; string temp;

					blockstream >> temp >> strain_contig >> temp >> temp >> temp >> temp >> sequence;
					istringstream strainstream(strain_contig);
					string strain;
					getline(strainstream,strain,'.');

					newblock.strain.push_back(strain);
					newblock.sequence.push_back(sequence);
				}
				block.push_back(newblock);
			}
		}
	}
}
					
void block_snp(vector<BLOCK> block, vector<SNP> &snp, int allow_indel_bool)
{
	for(int i=0; i<block.size(); i++)
	{
		for(int bl=0; bl<block[i].sequence[0].size(); bl++) // size of the block sequence
		{
			SNP newsnp; //only add if poly_bool==1
			
			int init_type=nt_type(block[i].sequence[0][bl]);

			if(allow_indel_bool==0) // do not look at indels
			{
				if(init_type==4) continue; 
			}

			int poly_bool=0;

			for(int si=0; si<block[i].strain.size(); si++) // traverse all strains
			{
				if(allow_indel_bool==0)// do not look at indels
				{
					if(nt_type(block[i].sequence[si][bl])==4) 
					{
						poly_bool=0;
						break;
					}
				}
				if(init_type!=nt_type(block[i].sequence[si][bl])) poly_bool=1;
				newsnp.pos=bl;
				newsnp.block=i;
				newsnp.src[block[i].strain[si]]=nt_type(block[i].sequence[si][bl]);
			}
			
			if(poly_bool==1) snp.push_back(newsnp);
		}
	}
}

void snp_vector(SNP s1, d_vec &v1) // input snp, output double vector for min_dis
{
	map<string,int>::iterator mi; 
	
	for(int type=0;type<5;type++) //each nt type gives one vector<int> for each snp
	{
		vector<int> vnt1;
		for(mi=s1.src.begin();mi!=s1.src.end();mi++) vnt1.push_back((mi->second)==type);
		v1.push_back(vnt1);
	}
}


void snp_cluster(vector<SNP> snp, vector<CLUSTER> &cluster, double threshold)
{
	vector<SNP> center;

	int si; int ci; //si-> present snp index, ci-> center snp index

	//initialize
	
	for(si=0; si<snp.size(); si++)
	{
		if(si%10000==1) 
		{
			cout<<"**"<<si<<"snps in "<<cluster.size()<<" clusters..."<<endl;
		}

		d_vec vsi;
		snp_vector(snp[si],vsi); //prepare dvec for min_dis computing

		double best_dis=-1.0; //the center that gives smallest distance
		int best_cluster=-1;

		for(ci=0;ci<cluster.size();ci++)
		{
			d_vec vci;
			snp_vector(cluster[ci].center,vci);

			double dis=min_dis(vsi, vci);

			if(dis<=threshold)
			{
				if(dis<best_dis||best_dis<0) 
				{
					best_dis=dis;
					best_cluster=ci;
				}
			}
		}
	
		if(best_dis<0) //no clustering assigned, generate new cluster
		{
			CLUSTER newcluster;
			newcluster.center.block=snp[si].block;
			newcluster.center.pos=snp[si].pos;
			newcluster.center.src=snp[si].src;


			cluster.push_back(newcluster);
		}
		else
		{
			cluster[best_cluster].snp.push_back(snp[si]);
		}
	}
}
				
void plot_prepare(CLUSTER cluster) //output data.frame for plotting

int main()
{

	ifstream maf("alignment.maf");
	vector<BLOCK> block;
	cout<<"*Parsing maf..."<<endl;
	maf_parse(maf,block,32);
	ofstream block_size("block_size.txt");
	
	for(int i=0;i<block.size();i++)
	{
		block_size<<i<<"\t"<<block[i].sequence[0].size()<<endl;
	}

	vector<SNP> snp;
	cout<<"*Finding snp..."<<endl;
	block_snp(block, snp, 0);
	vector< vector<int>> v1;
	
	cout<<"**"<<snp.size()<<" snps found."<<endl;

	vector<CLUSTER> cluster;
	

	double d=20.0;

	cout<<"*Clustering snp at d2="<<d<<"..."<<endl;
	snp_cluster(snp, cluster, d);

	cout<<"*Outputing results..."<<endl;
	
	ofstream snp_cluster("snp_clusters20.txt");
	ofstream block_fasta("block_fasta20.fa");
	ofstream cluster_size("cluster_size20.txt");

	cluster_size<<cluster.size()<<" clusters generated at d2="<<d<<endl;

	for(int i=0; i<cluster.size();i++)
	{
		for(int j=0; j<cluster[i].snp.size();j++)
		{
			snp_cluster<<i<<"\t"<<cluster[i].snp[j].block<<"\t"<<cluster[i].snp[j].pos<<endl;
		}
	}

	for(int i=0; i<block.size();i++)
	{
		for(int j=0; j<block[i].strain.size();j++)
		{
			block_fasta<<">Block|"<<i<<"|strain|"<<block[i].strain[j]<<endl;
			block_fasta<<block[i].sequence[j]<<endl;
		}
	}




	
	
}
