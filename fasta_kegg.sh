#PBS -l walltime=12:00:00,mem=16gb,nodes=1:ppn=1

module load khmer/2.0

module load usearch
usearch -ublast block_fasta_protein_formatted.faa  -db /home/zhouw/Kegg_db/species_prokaryotes.udb -evalue 1e-9 -accel 0.5 -maxhits 1 -userout kegg_blast_output.txt -userfields query+target+evalue+id

./kegg_ko $1 "/home/zhouw/Kegg_db/ko_genes.list"
