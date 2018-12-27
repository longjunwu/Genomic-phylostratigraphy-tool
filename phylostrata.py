# This script is used to determine the phylostrata of each gene in a genome
# Input is a collection of protein sequences of a genome in fasta format
# Output is the phylostrata of each gene
# required module: biopython
# required NCBI blast software installed; blast database: NCBI nr database downloaded and indexed locally

# The pipeline is divided into three main steps.

##########Step 1: BLAST SEARCH##########

# The purpose of this step is to blast each protein seq in the query genome file to the target nr database
# the pipeline then appends the blast results of each run to the output file
# results from each query sequence are separated by the characters '!@~'


import os
import sys

# import biopython
sys.path.append("link/to/biopython/")
# import bio.blast from biopython to handle blast results
from Bio.Blast import NCBIXML

# open input file (protein sequences of a genome) and store the input file in "all_seq"
all_seq = open('/link/to/input/file/input_file', "r")
# set up output file "blast_result_file", which will contain blast results
blast_result_file = open("blast_result_file", "a")
# split sequences in the input file, and store each sequence in a list "all_seq_split"
all_seq_split = all_seq.split('>')[1:]

# loop through each sequence in the sequence list
for sequence in all_seq_split:
    # the string lacks the symbol ">", add it here
    sequence = '>' + sequence
    # create temporary file for writing each query sequence
    temp_query_seq = open("temp_query_seq", "w")
    # write one query sequence string into file "temp_query_seq"
    temp_query_seq.write(sequence+'\n')

    # start blastp using temp_query_seq as the query against the nr database, and store results (xml format) in
    # blast_result
    blast_result = os.popen(
        'blastp  -outfmt "6 qseqid sseqid evalue  sgi sacc staxids stitle" -evalue 0.01 -query %s -db '
        '/link/to/nr/database/%s' % (temp_query_seq, nr)).read()
		
    # blast result is stored in "blast_result"
    # the e-value can be set to different values as desired, it is set as 0.01 here

    # store the blast result in the output file "blast_result_file"
    blast_result_file.write(blast_result+'\n')
    # add a "marker" at the end of the blast result of each query sequence
    blast_result_file.write('!@~'+'\n')
	
	

################### STEP 2: EXTRACT TAXONOMY LINEAGE FROM THE BLAST RESULTS #######################


# In this step, the NCBI taxonomic ID [Accession Number] of each hit in the blast results will be extracted
# This ID is then converted into Taxonomy Lineage using NCBI taxonomic information database "gi_taxid_prot.dmp"
# example of taxonomy lineage: "cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;
# Lophotrochozoa;Mollusca;Bivalvia;Pteriomorphia;Ostreoida;Ostreoidea;Ostreidae;Crassostrea;Crassostrea gigas;"
# This pipeline requires gi_taxid_prot.dmp, downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
# gi_taxid_prot.dmp could be downloaded and prepared automatically from a bash shell "download.sh",
# which is attached at the end of this script

# open output file of blast results from last step, and store in "all_result"
all_result = open('blast_result_file', "r")

# create output file for this step, which will contain taxonomy lineages of each hit from each blast result
outfile = open("taxon_info_output_file", "a")

# split the blast results file into a list of blast results, each element containing blast results from each
# query sequence(gene), and store all of them in a list "all_result_split"
all_result_split = all_result.split("!@~")[:-1]

# loop through each blast result query
for result in all_result_split:
    # split blast results by hit, each hit being on its own line
    result_split = result.split("\n")[:-1]
    # extract the NCBI taxonomic ID from each hit
    result_id = result_split[0].split('	')[0]
    # write the taxonomic ID to the output file, each blast query separated by '!@~'
    outfile.write("!@~"+result_id+'\n')

    # loop through each hit in each blast result query
    for hit in result_split:
        # split each hit line into separate column elements, with each column element representing specific information
        # including hit id, e-value, and hit score and so on
        hit_split = hit.split('	')

        # set the e-value cutoff here; cutoff can be set to any number smaller than the e-value cutoff in the blast
        # search, thus when changing the blast cut off, blast does not need to be re-run
        # hit_split [2] stores the e-value of each hit
        if float(hit_split[2]) < 0.00001:
            # extract the taxonomic id of the significant hit
            hit_taxon_ids = hit_split[5]
            # some hits have more than one id separated by ";" ; all id's are processed
            hit_taxon_ids_split = hit_taxon_ids.split(';')
            for hit_taxon_id in hit_taxon_ids_split:
                # the script below uses the shell script 'taxon_tid.sh' to convert each taxonomic id into
                # its corresponding NCBI taxonomy lineage
                # each significant hit's taxonomy lineage is written to the output file "taxon_info_output_file"
                # the shell script "taxon_tid.sh" is attached at the end of this script for reference
                os.system("bash /link/to/taxon_tid.sh %s >> taxon_info_output_file" % hit_taxon_id)
				

##########STEP 3: EXTRACT PHYLOSTRATA FROM TAXON LINEAGE INFORMATION##########


# Here we use Crassostrea gigas (oyster) as an example. If a different query species is used, adjust accordingly
# "oyster_taxon" stores the taxonomy lineage of the Crassostrea gigas from NCBI
oyster_taxon = "cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Lophotrochozoa;" \
               "Mollusca;Bivalvia;Pteriomorphia;Ostreoida;Ostreoidea;Ostreidae;Crassostrea;Crassostrea gigas;"
# split the oyster taxonomy lineage in "oyster_taxon", it will be used later when we loop through each hit
oyster_taxon_split = oyster_taxon.split(';')[:-1]

# open the taxonomy lineage file from last step
taxon_result = open('taxon_info_output_file', "r")

# split taxonomy lineage from each blast query in taxon_result and store in list "taxon_result_split"
taxon_result_split = taxon_result.split('!@~')[1:]

# the taxonomy lineage from each query oyster gene from blast is analysed separately
for query_result in taxon_result_split:
    # extract the query id
    query_id = query_result.split('\n')[0]

    # Each hit in the taxonomy lineage from each query gene is analysed separately
    for hit in query_result.split('\n')[1:]:
        # hit level is used to show the taxonomy level of each hit
        # at the start of each hit, set initial hit level to 16, representing the youngest possible taxon level
        # this initial hit level could be adjusted according to the query species
        hit_level = 16

        ###########below start the processing of the hits################

        # loop through each taxonomic level of the oyster (the oyster_taxon variable), ("Metazoa", for example),
        # starting with the oldest taxon and ending at the youngest (species-specific)
        # after each iteration of a taxon in oyster_taxon, hit_level (phylostrata) of each hit will increase by 1
        # this iterates until the taxon in taxonomy lineage of oyster is not found in the taxonomy lineage of the hit
        # the loop breaks, and hit_level represents the lowest (youngest) taxon of that hit


        # if "cellular organisms" is not found, the hit is not a cellular organism, we do not process these hits any
        # further, similar to other studies
        # therefore before analysis check if "cellular" is found, ie. the hit is a cellular organism
        # for "environmental samples" or "uncultured", the taxonomy is unclear; these hits will not be processed any
        # further similar to other studies

    if hit.find("cellular") >= 0 and hit.find("environmental") < 0 and hit.find("uncultured") < 0:

        # for each level in oyster taxon
        # initialize i
        i = 0

        # start the loop
        for level in oyster_taxon_split:
            i = i + 1

            # Platyhelminthes (flatworm) is not listed as "lophotrochozoa" in NCBI taxonomy, we adjust it here
            if hit.find("Platyhelminthes") >= 0:
                # if hit matches to "Platyhelminthes" then the phylostrata will be adjusted to "8" (the level of
                # "lophotrochozoa")
                hit_level = 8
                # do not need to process further for the hit, stop the loop
                break

            # if the oyster taxon cannot be found in the taxonomy lineage of the hit, the hit is no longer shared with
            # the oyster taxonomy lineage
            if hit.find(level) < 0:
                # for example, if the oyster taxonomy lineage "Metazoa" cannot be found in the taxonomy lineage of
                # the hit, this means the species of the hit is not "Metazoa"
                # since the gene must already exist in one phylostrata level lower: Opisthokonta, which has a level of
                # 3, hit_level = i-1
                hit_level = i - 1
                # stop the loop
                # the current value of i represents the lowest (youngest) taxon of that hit
                break

        # if hit_level of the current hit is smaller than hit_level of all the previous hits processed in the same
        # query gene of the blast results
    if hit_level < hit_level_lowest:
        # assign the value of hit_level to hit_level_lowest
        hit_level_lowest = hit_level

    # after processing all hits in blast results from a single query gene, "hit_level_lowest" represents the
    # phylostrata of that gene
    # assign the id and phylostrata of the query gene to "output_result_phylostrata"; separate id and phylostrata by a
    # ":", separate each query gene result by a ";"
    output_result_phylostrata = query_id + ':' + str(hit_level_lowest) + ';'

    # print result, each query gene on its own line
    print
    output_result_phylostrata

########## APPENDIX ##########

# APPENDIX 1:below is the shell script "download.sh" that download and prepare "gi_taxid_prot.dmp" from NCBI

'''
#!/bin/bash
# from https://www.biostars.org/p/13452/, the author is Frederic Mahe, Kaiserslautern

## Download NCBI's taxonomic data and GI (GenBank ID) taxonomic
## assignation.

## Variables
NCBI="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
TAXDUMP="taxdump.tar.gz"
TAXID="gi_taxid_nucl.dmp.gz"
NAMES="names.dmp"
NODES="nodes.dmp"
DMP=$(echo {citations,division,gencode,merged,delnodes}.dmp)
USELESS_FILES="${TAXDUMP} ${DMP} gc.prt readme.txt"

## Download taxdump
rm -rf ${USELESS_FILES} "${NODES}" "${NAMES}"
wget "${NCBI}${TAXDUMP}" && \
    tar zxvf "${TAXDUMP}" && \
    rm -rf ${USELESS_FILES}

## Limit search space to scientific names
grep "scientific name" "${NAMES}" > "${NAMES/.dmp/_reduced.dmp}" && \
    rm -f "${NAMES}" && \
    mv "${NAMES/.dmp/_reduced.dmp}" "${NAMES}"

## Download gi_taxid_nucl
rm -f "${TAXID/.gz/}*"
wget "${NCBI}${TAXID}" && \
    gunzip "${TAXID}"

exit 0
'''

# APPENDIX 2 below is the shell script of taxon_tid.sh
'''		 

#!/bin/bash
#modified from https://www.biostars.org/p/13452/, the author is Frederic Mahe, Kaiserslautern


NAMES="/scratch/lwu10/hourglass/ncbi_taxon_info/names.dmp"
NODES="/scratch/lwu10/hourglass/ncbi_taxon_info/nodes.dmp"
GI_TO_TAXID="/scratch/lwu10/hourglass/ncbi_taxon_info/gi_taxid_prot.dmp"
TAXONOMY=""
TAXID="${1}"

# Obtain the name corresponding to a taxid or the taxid of the parent taxa
get_name_or_taxid()
{
    grep --max-count=1 "^${1}"$'\t' "${2}" | cut --fields="${3}"
}

# Get the taxid corresponding to the GI number
#TAXID=$(get_name_or_taxid "${GI}" "${GI_TO_TAXID}" "2")

# Loop until you reach the root of the taxonomy (i.e. taxid = 1)
while [[ "${TAXID}" -gt 1 ]] ; do
    # Obtain the scientific name corresponding to a taxid
    NAME=$(get_name_or_taxid "${TAXID}" "${NAMES}" "3")
    # Obtain the parent taxa taxid
    PARENT=$(get_name_or_taxid "${TAXID}" "${NODES}" "3")
    # Build the taxonomy path
    TAXONOMY="${NAME};${TAXONOMY}"
    TAXID="${PARENT}"
done

echo -e "${TAXONOMY}"

exit 0

'''


