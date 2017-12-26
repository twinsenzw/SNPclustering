#PBS -l walltime=72:00:00,nodes=1:ppn=1,mem=48gb

module load khmer/2.0
cd /home/zhouw/sepi/SNP_mapping

./SNP_cluster
