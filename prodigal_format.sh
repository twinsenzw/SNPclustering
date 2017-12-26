#PBS -l walltime=48:00:00,mem=48gb,nodes=1:ppn=1

module load khmer/2.0
cd /home/zhouw/sepi/SNP_mapping

#--1. run prodigal--
/home/zhouw/prodigal/prodigal.linux -i block_fasta.fa -a block_fasta_protein.faa

#--2. format header from sequence position to block sequence position--
# prodigal output: 1-based, snp_cluster output: 0-based, considering indel positions

./prodigal_format
cat block_fasta_protein_formatted.faa | tr -d '*' > temp
mv temp block_fasta_protein_formatted.faa 

#--3. blast to kegg--
module load usearch
usearch -ublast block_fasta_protein_formatted.faa  -db /home/zhouw/Kegg_db/species_prokaryotes.udb -evalue 1e-9 -accel 0.5 -maxhits 1 -userout kegg_blast_output.txt -userfields query+target+evalue+id
