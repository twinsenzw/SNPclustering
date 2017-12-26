#PBS -l walltime=64:00:00,mem=96gb,nodes=1:ppn=1

module load khmer/2.0

cd /home/zhouw/sepi/SNP_mapping
./blast_format 143
