
Example usage

input_file_path = "1148_MHCI/1148_MHCI/1148_MHCI.xlsx"
MHC_class = '1'
gibbs_cluster = [0,1,2,3,4]
file_name = '1148_MHCI'
output_file_path = "D:/Period_6/out/" + file_name  + "/"    # full file path, has to have "/" at the end

immuno_search(work_dir, input_file_path, MHC_class, gibbs_cluster, output_file_path, file_name)

work_dir
	working directory where run-time files will be created
	should have a subdirectory "/db" with all the required databases

input_file_path
	Full path to excel file
	Peptide table from PEAKS sheet same name as sample
	Gibbs clustering results in excel sheet name "gibbs_clustering"
MHC_class
	"1", "2" or "E"
gibbs_cluster	
	array with numbers - [0,1,2,3,..]
	if empty, all clusters are selected 
output_file_name
	full output directory path with a "/" at the end
	

Dependencies
	Python - pandas, numpy, csv, SeqIO from Biopython
	NCBI blast or blast from biopython
	Seqkit

FASTA files for blast database and parsing

	Known HLAs - https://peptideatlas.org/builds/human/hla/202311/APD_Hs_all.fasta

	Canonical proteins - uniprot canonical with isoforms or Refseq assembly proteins

	Single amino acid polymorphisms - http://119.3.41.228/dbSAP/download.html

	Six-frame translated protein - translated GRCh38p14_genomic using seqkit 2.8.2 
	
	seqkit translate db/GRCh38.p14_genomic.fna --append-frame -x -f 6 -M -m 8 --transl-table 	1 -s > db/human_6FT_m.fasta


Making blast databases

Command line queries for making required blast databases

#make databse from known HLA peptides fasta file
makeblastdb -in db/APD_Hs_all.fasta -dbtype prot -out db/APD_Hs_all

#make human canonical protein blast database from fasta
makeblastdb -in db/human_canonical.fasta  -dbtype prot -parse_seqids -out db/human_canonical

#make SAAV blast database from fasta
makeblastdb -in db/sap_db.fa  -dbtype prot -parse_seqids -out db/sap_db

#sixframe 

