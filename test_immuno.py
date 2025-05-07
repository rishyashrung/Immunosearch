import pandas as pd
import numpy as np
import csv
from Bio import SeqIO
import os
import subprocess
import getopt
import sys


#functions 

#fn for creating files from a list and file_name
def create_files(pep, file_name):

    ofile = open(file_name + ".fasta", "w")

    for i in range(len(pep)):
        ofile.write('>' + '\n' + pep[i] + '\n')
    ofile.close()

    # peptide list for parsing
    ofile = open(file_name + '_2.fasta', "w")

    for i in range(len(pep)):
        ofile.write(pep[i] + '\n')
    ofile.close()
    
    return print("\t" + "files", file_name, ",", file_name + '_2', "created" )

#fn for parsing and categorizing blast output
def parse_categorize(database_fasta, blast_out, fasta2): #requires inputs with file extension

    input1= SeqIO.parse(database_fasta,"fasta") # blastp reference database
    seqdb={}
    for record in input1:
        seq=str(record.seq)
        if record.id not in seqdb:
            seqdb[record.id]=seq

    input2= open(blast_out,"r") # blastp output
    input3= open(fasta2,"r") # novel peptide tab txt
    output= open('categorized_' + blast_out,"w")


    blastout={}
    hits_dic={}
    for line in input2:
        row=line.strip().split("\t")
        qid=row[-2]
        sid=row[1]
        sseq=seqdb[sid]
        ident=row[2]
        peplen=int(row[3])
        mismatch=row[4]
        alignlen=int(row[6])-int(row[5])+1
        sstart=int(row[7])
        send=int(row[8])
        gap=row[11]
        evalue=float(row[-5])
        alignseq=row[-1]
        category="NA"
        single_sub_pos="NA"
        if sstart>3:
            Nterm_seq=sseq[sstart-4:sstart+2] #check up 3 amino acid before N-term of this peptide
        else:
            Nterm_seq=sseq[:sstart]

        if len(sseq)-send<3:
            Cterm_seq=sseq[send-1:]
        else:
            Cterm_seq=sseq[send-3:send+3]

        if alignlen==peplen:
            if float(ident)==100:
                category="match to known protein"
            
            elif int(gap)==0 and int(mismatch)==1:
                category="map to known protein with 1 aa mismatch"
                for i in range(peplen):
                    if qid[i]!=alignseq[i]:
                        single_sub_pos=str(i+1)

            elif int(gap)==1 and int(mismatch)==0:
                category="map to known protein with 1 aa insertion"
            else:
                category="novelpep (map to known protein with more than 2 mismatched aa)"
        elif peplen-alignlen==1 and float(ident)==100:
            category="map to known protein with 1 aa deletion"

        else:
            category="novelpep (map to known protein with more than 2 mismatched aa)"
        
        if qid not in hits_dic:
            hits_dic[qid]=evalue
            blastout[qid]=[category,sid,ident,peplen,single_sub_pos,Nterm_seq,alignseq,Cterm_seq,alignlen,mismatch,gap]
        else:
            if evalue<hits_dic[qid]:
                hits_dic[qid]=evalue
                blastout[qid]=[category,sid,ident,peplen,single_sub_pos,Nterm_seq,alignseq,Cterm_seq,alignlen,mismatch,gap]

    #header=input3.readline().strip().split("\t")

    header=["Query","blastp_category","blastp_match","identity","peplen","sub_pos","Nterm-seq(3aa)","aligned_seq","Cterm-seq(3aa)","alignlen","mismatch","gap"]
    output.write("\t".join(header)+"\n")

    for line in input3:
        row=line.strip().split("\t")
        peptide=row[0]
        if peptide in blastout:
            results=blastout[row[0]]
            newrow=row+results
            output.write("\t".join(map(str,newrow))+"\n")
        else:
            newrow=row+["novelpep (no match to known protein found)","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"]
            output.write("\t".join(map(str,newrow))+"\n")


    input2.close()
    input3.close()
    output.close()
   
    #adding headers to the novel blast output
    header_names = ["qseqid", "sseqid", "pident", "qlen", "mismatch", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "gaps", "qseq", "sseq"]
    novel_file = pd.read_table('categorized_' + blast_out, sep ='\t', names = header_names)
    header_names = ["qseqid", "sseqid", "pident", "qlen", "mismatch", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "gaps", "qseq", "sseq"]
    blast_file = pd.read_table(blast_out, sep ='\t', names = header_names)
    blast_file
    #to check if any keys are missing in the dictionary seqdb[]
    missing = blast_file[~blast_file["sseqid"].isin(seqdb.keys())]
    if (len(missing) == 0):
        return print("\t"+ 'all is gud for',blast_out)
    else:
        return print('dictionary error') #make sure the dictionary key is the same as the accession in
    

def get_pos(substring, string):
# Find the starting position
    start_position = string.find(substring)

    if start_position != -1:
        # Calculate the ending position
        end_position = start_position + len(substring) - 1
        return(start_position, end_position)


#fn that gets matches and position of matches from peptide_list to fasta
def find_matches(peptide, db, query):
    search_6ft = {}
    pos_6ft = {}
    for string in peptide:
        search_6ft[string] = []
        for key, value in db.items():
            if string in value:
                if (query == '6FT'):
                    combi = (key, get_pos(string,value))
                    search_6ft[string].append(combi)
                elif (query == 'PCPS'):
                    search_6ft[string].append(key)
                else:
                    return print("specify find matches query: 6FT or PCPS")
    unmatched = []
    #fn for filtering the 6FT matched dictionary
    def my_filtering_function(pair):
        key, value = pair
        if value == []:
            return False  # filter pair out of the dictionary
        else:
            return True  # keep pair in the filtered dictionary
        
    matched = dict(filter(my_filtering_function, search_6ft.items()))
    positions = dict(filter(my_filtering_function, pos_6ft.items()))

    for key, value in search_6ft.items():
        if (value == []):
            unmatched.append(key)
    if (len(matched) + len(unmatched) == len(search_6ft)):  #to make sure the dictionary filter works fine
        return (matched, unmatched)
    else:
        return print("\t" + "error in separating matches")
    
def split_string(s):
    splits_dict = {}
    # Ensure the string is long enough to be split into two parts each with at least `min_length` characters
    if len(s) < 6:
        return splits_dict
    
    for i in range(3, len(s) - 2):
        part1 = s[:i]
        part2 = s[i:]
        if len(part2) >= 3:
            splits_dict[part1] = part2
    return splits_dict

def PCPS(input_file):
        
    input1= SeqIO.parse('db/human_canonical.fasta',"fasta") 
    seqdb={}
    for record in input1:
            seq=str(record.seq)
            if record.id not in seqdb:
                    seqdb[record.id]=seq
    
    input = open(input_file, "r")
    output = open('cis_PCPS', "w")
    output_2 = open('trans_PCPS', "w")
    for line in input:
            row=line.strip().split("\t")
            pep = row[0]
            split_combinations = split_string(pep)
            part1_list = split_combinations.keys()
            part2_list = split_combinations.values()
            match_1, unmatch_1 = find_matches(part1_list, seqdb, 'PCPS')
            match_2, unmatch_2 = find_matches(part2_list, seqdb, 'PCPS')
            row = []
            for key,value in match_1.items():
                    for pro_1 in match_1[key]:
                            if (split_combinations[key] in match_2.keys()):
                                    for pro_2 in match_2[split_combinations[key]]:
                                            if pro_1 == pro_2:
                                                result = (f"{key}|{split_combinations[key]}" + "\t" + f"{pro_2}")
                                                output.write(str(pep) + "\t" + result + "\n")
                                            else:
                                                result_2 = (f"{key}|{split_combinations[key]}" + "\t" + f"{pro_1}:{pro_2}")
                                                output_2.write(str(pep) + "\t" + result_2 + "\n")
    input.close()
    output.close()
    output_2.close()
    return print("\t" + "separate files cis_PCPS and trans_PCPS created for spliced peptides")



if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "w:i:m:g:o:f:", 
                                   ["work_dir=", "input_file_path=", "MHC_class=", 
                                    "gibbs_cluster=", "output_file_path=", "file_name="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    work_dir = input_file_path = MHC_class = gibbs_cluster = output_file_path = file_name = None

    for opt, arg in opts:
        if opt in ("-w", "--work_dir"):
            work_dir = arg
        elif opt in ("-i", "--input_file_path"):
            input_file_path = arg
        elif opt in ("-m", "--MHC_class"):
            MHC_class = str(arg)
        elif opt in ("-g", "--gibbs_cluster"):
            gibbs_cluster = [int(x) for x in arg.split(',')]  # Convert comma-separated string to a list of integers
        elif opt in ("-o", "--output_file_path"):
            output_file_path = arg
        elif opt in ("-f", "--file_name"):
            file_name = arg

    # Check if any required arguments are missing
    if None in (work_dir, input_file_path, MHC_class, gibbs_cluster, output_file_path, file_name):
        print("Missing required arguments")
        sys.exit(2)

    # Run the main function with the provided arguments

    os.chdir(work_dir)
    print("loading into work directory at", work_dir)
    
    # peptide length 8-11, removing gibbs junk and DB matched
    all_data = pd.ExcelFile(input_file_path)
    data = pd.read_excel(all_data, file_name)
    gibbs = pd.read_excel(all_data, 'gibbs_clustering')
    #data = pd.read_csv(path_peptides) #PSM table from PEAKS
    #gibbs = pd.read_csv(path_gibbs) #Gibbs clustering CSV
    print("filtering peptides found by DB search")
    
    data = data[data["Found By"] != 'DB Search']
    if (len(gibbs_cluster) != 0):
        print("Selected gibbs clusters based on input, gibbs clusters", gibbs_cluster)
        gibbs = gibbs[gibbs["Gn"].isin(gibbs_cluster)]
    else:
        print("selecting all gibbs clusters")
        gibbs = gibbs
    
    if (MHC_class == '1'):
        print("selecting peptides with length between 8 and 11 AA")
        gibbs = gibbs[gibbs["Sequence"].str.len().between(8,11)]
    elif (MHC_class == '2'):
        print("selecting peptides with length between 12 and 17 AA")
        gibbs = gibbs[gibbs["Sequence"].str.len().between(12,17)]
    elif (MHC_class == "E"):
        print("selecting peptides with length between 8 and 15 AA")
        gibbs = gibbs[gibbs["Sequence"].str.len().between(8,15)]
    else:
        print("!!! not filtered by length")

    data = data[data["Peptide"].isin(gibbs["Sequence"])]
    list = pd.unique(data["Peptide"])
    pep = list.tolist()
    blast_p = pep
    print("generating lists and files for further analysis")
    create_files(pep,'peptides')


    #blast against known HLA peptides
    if(len(pep) != 0):

        if(len(pep) == 1): #coz I'm a grammar Nazi :)
            print(str(len(pep)) + " peptide being searched for known HLAs")
        else:
            print( str(len(pep)) + " peptides being searched for known HLAs")
        
        #conda run -n Bio-stats blastp -task blastp-short -query peptides.fasta -db db/APD_Hs_all -out HLA_blast_out  -evalue 10.0 -outfmt \"6 qseqid sseqid pident qlen mismatch qstart qend sstart send evalue bitscore gaps qseq sseq\"
        subprocess.run("blastp -task blastp-short -query peptides.fasta -db db/APD_Hs_all -out HLA_blast_out -evalue 10.0 -outfmt \"6 qseqid sseqid pident qlen mismatch qstart qend sstart send evalue bitscore gaps qseq sseq\"", shell=True)

        #parsing and catergorizing HLA_blast output
        print("\t"+ "reading output")
        parse_categorize('db/APD_Hs_all.fasta', 'HLA_blast_out', 'peptides_2.fasta')

        #known HLA
        print("\t"+ "generating lists and files for further analysis")
        output_HLA = pd.read_table('categorized_HLA_blast_out')
        known = output_HLA[output_HLA["blastp_category"] == 'match to known protein']
        list = known["Query"] 
        pep = list.to_list()

        if (len(pep) != 0):
            print("\t"+ "known HLA found")
            ofile = open("known_HLA.fasta", "w")

            for i in range(len(pep)):
                ofile.write('>' + '\n' + pep[i] + '\n')
            ofile.close()
        else:
            print("\t"+ "no known HLA found")

        known = output_HLA[output_HLA["blastp_category"] != 'match to known protein']
        list = known["Query"] 
        pep = list.to_list()
    else:
        print("No peptides to search for known HLAs")

    
    if (len(pep) != 0): #to prevent blast with empty query
        
        create_files(pep, 'to_blastp')

        if (len(pep) == 1):
            print(str(len(pep)) + " peptide being searched for human canonical proteins")
        else:
            print(str(len(pep)) + " peptides being searched for human canonical proteins")
        
        #blast all proteins against human canonical proteins
        subprocess.run("blastp -task blastp-short -query to_blastp.fasta -db db/human_canonical -out blastp_out_human_canonical -evalue 10.0 -outfmt \"6 qseqid sseqid pident qlen mismatch qstart qend sstart send evalue bitscore gaps qseq sseq\"", shell=True)
        
        print("\t"+ "reading output")
        parse_categorize('db/human_canonical.fasta', 'blastp_out_human_canonical', 'to_blastp_2.fasta')
        output  = pd.read_table('categorized_blastp_out_human_canonical')
        
        #preparing fasta files of proteins to search 6FT database
        
        to_6ft = output[(output["blastp_category"] != 'map to known protein with 1 aa mismatch') & (output["blastp_category"] != 'match to known protein')]
        list = to_6ft["Query"] 
        pep = list.to_list()

        create_files(pep, 'to_6ft')
       
        #preparing fasta files of proteins with 1 AA mismatch
        mismatched = output[output["blastp_category"] == 'map to known protein with 1 aa mismatch']
        list = mismatched["Query"] 
        pep = list.to_list()
        SAAV = pep
    else:
        print("blastp input empty")

    print("\t"+ "generating lists and files for further analysis")
    
    
    if (len(pep) != 0): #to prevent blast with empty query

        create_files(pep,'SAAV')
        
        if (len(pep) == 1):  #more grammar Nazi
            print(str(len(pep)) + " peptide being searched for Single Amino Acid Variants")
        else:
            print(str(len(pep)) + " peptides being searched for Single Amino Acid Variants")
            
        #blastp against db_SAP (single amino acids polymorphisms)
        #conda run -n Bio-stats blastp -task blastp-short -query SAAV.fasta -db db/sap_db -out blast_out_SAAV -evalue 10.0 -outfmt \"6 qseqid saccver pident qlen mismatch qstart qend sstart send evalue bitscore gaps qseq sseq\"
        subprocess.run("blastp -task blastp-short -query SAAV.fasta -db db/sap_db -out blast_out_SAAV -evalue 10.0 -outfmt \"6 qseqid sseqid pident qlen mismatch qstart qend sstart send evalue bitscore gaps qseq sseq\"", shell=True)

        print("\t"+ "reading output")
        parse_categorize('db/sap_db.fa', 'blast_out_SAAV', 'SAAV_2.fasta')


        output = pd.read_table("categorized_blast_out_SAAV")
        not_SNP = output[output["blastp_category"] != 'match to known protein']
        pep = not_SNP["Query"]
        pep = pep.to_list()

        ofile = open("not_SAAV.fasta", "w")

        for i in range(len(pep)):
            ofile.write('>' + '\n' + pep[i] + '\n')
        ofile.close()

    else:
        print("no potential SAAVs, proceeding to 6FT search")

    #searches peptides in 6FT db
    
    pep = to_6ft["Query"]

    if(len(pep) != 0):

        if(len(pep) == 1):
            
            print(str(len(pep)) + " peptide being searched in the six frame translated human genome")
            
            
            subprocess.run("seqkit grep --by-seq --ignore-case --threads 12 --seq-type protein --pattern-file to_6ft_2.fasta db/human_6FT_m.fasta > 6ft_out", shell=True)


            input1= SeqIO.parse('6ft_out',"fasta") # 6FT results to dict
            seqdb={}
            for record in input1:
                seq=str(record.seq)
                if record.description not in seqdb:
                    seqdb[record.description]=seq
            
            if (len(seqdb) == 0):
                print("no match to 6FT")
                unmatched = pep
                matched = []
            else:
                matched = pep
                unmatched = []

        else:
            print(str(len(pep)) + " peptides being searched in the six frame translated human genome")

            
            subprocess.run("seqkit grep --by-seq --ignore-case --threads 12 --seq-type protein --pattern-file to_6ft_2.fasta db/human_6FT_m.fasta > 6ft_out", shell=True)


            print("\t"+ "writing 6FT results to dictionary")

            input1= SeqIO.parse('6ft_out',"fasta") # 6FT results to dict
            seqdb={}
            for record in input1:
                seq=str(record.seq)
                if record.description not in seqdb:
                    seqdb[record.description]=seq
                    
            print("\t"+ "number of matches from 6FT:", len(seqdb))

            input  = open('to_6ft_2.fasta', 'r') #list of peptides, has to be in a list
            pep = []
            for line in input:
                    pep = pep + line.strip().split("\t")

            print("\t"+ "locating and matching peptides to 6FT")

            matched, unmatched = find_matches(pep, seqdb, '6FT')

            if len(matched) == 0:
                print("\t"+ "no peptides matched to 6ft")
            else:
                print("\t"+ "writing " + str(len(matched)) + " peptide matches to 6FT")
                with open('matches_to_6ft.csv','w') as f:
                    w = csv.writer(f)
                    w.writerow(["Peptide", "Sequence_loc"]) #header names
                    for key in matched:
                        sequence = matched[key]
                        w.writerow([key, sequence])

            print("\t"+ "writing " + str(len(unmatched)) + " unmatched peptides")
            create_files(unmatched, 'no_match_6ft')

    else:
        print("no peptides to search in 6FT")
    

    #searching for proteosome catalyzed peptide spliced 
    
    print(str(len(unmatched)) + " peptides being searched for proteosome catalyzed peptide spliced variants")

    PCPS("no_match_6ft_2.fasta")
    
    
    #saving output files
    SAAV_out = pd.read_table('categorized_blast_out_SAAV', sep = "\t")
    Sixframe_out = pd.read_csv('matches_to_6ft.csv')
    Sixframe_notmatched = pd.read_table('no_match_6ft_2.fasta', sep = "\t", names = ["Peptides"])
    cis_PCPS = pd.read_table("cis_PCPS", names= ["Peptide", "Spliced_peptide", "Protein_origin"])
    #writing data to excel file
    if (len(mismatched) == 0 & len(matched) == 0 & len(unmatched) == 0):
        print("No significant peptides found")
        
    else:
        
        #creates output directory
        try:
            os.makedirs(output_file_path)
            os.chdir(output_file_path)
            print("Search directory created at " + os.getcwd())
        except OSError as error:
            os.chdir(output_file_path)
            print("Directory already exists at " + os.getcwd())
            
        print("Creating excel file with results")
        with pd.ExcelWriter(file_name + '_immuno_search_out.xlsx', engine='openpyxl') as writer:
            # Write each DataFrame to a different sheet
            if (len(mismatched) != 0):
                SAAV_out.to_excel(writer, sheet_name='Single_AA_variants', index=False)
            if (len(matched) != 0):
                Sixframe_out.to_excel(writer, sheet_name='Matches_to_six_frame', index=False)
            if (len(unmatched) != 0):
                Sixframe_notmatched.to_excel(writer, sheet_name='Six_frame_non_matched', index=False)
            if (len(cis_PCPS) != 0):
                cis_PCPS.to_excel(writer, sheet_name= "cis_PCPS", index= False)

        os.chdir(work_dir)
        print("returning back to work directory at", work_dir)
        print("search and classification done, file saved at ", output_file_path + file_name + '_immuno_search_out.xlsx')
