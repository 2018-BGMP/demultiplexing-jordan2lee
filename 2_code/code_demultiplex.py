#!/usr/bin/env python3

#use below for talapas
#!/usr/bin/env python

##############################################
# Demultiplex: create and count index hopping
#bi624 assignment demultiplex
##############################################

# Method 1:    argparse

import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='demultiplex program for paired end reads and count instances of index hopping')
    parser.add_argument("-R1", "--read1_fq", help ="path to forward reads fastq", required=True, type=str)
    parser.add_argument("-R2", "--read2_fq", help ="path to reverse reads fastq", required=True, type=str)
    parser.add_argument("-I1", "--index1_fq", help ="path to forward index fastq", required=True, type=str)
    parser.add_argument("-I2", "--index2_fq", help ="path to reverse index fastq", required=True, type=str)
    parser.add_argument("-rds", "--total_reads", help ="total number of reads in a given input file. note this number needs to be the same across all input files", required=True, type=int)
    return parser.parse_args()
    
args = get_arguments()
read1 = args.read1_fq
read2 = args.read2_fq
index1 = args.index1_fq
index2 = args.index2_fq
num_reads = args.total_reads

########################################## or use ##############################################################

# Method 2:    via objects (real Talapas files)
#read1 = "/projects/bgmp/jlee33/bi624_gen_lab/2_demultiplex/data/1294_S1_L008_R1_001.fastq"
#read2 = "/projects/bgmp/jlee33/bi624_gen_lab/2_demultiplex/data/1294_S1_L008_R4_001.fastq"
#index1 = "/projects/bgmp/jlee33/bi624_gen_lab/2_demultiplex/data/1294_S1_L008_R2_001.fastq"
#index2 = "/projects/bgmp/jlee33/bi624_gen_lab/2_demultiplex/data/1294_S1_L008_R3_001.fastq"

#num_reads = 363246735

# Method 2:    via objects (testfiles)
#read1 = "/mnt/c/Users/Jordan/UO_Docs/Genomics_lab/assignments/demultiplexing-jordan2lee/2_code/test_data/R1_testfile.fq"
#read2 = "/mnt/c/Users/Jordan/UO_Docs/Genomics_lab/assignments/demultiplexing-jordan2lee/2_code/test_data/R2_testfile.fq"
#index1 = "/mnt/c/Users/Jordan/UO_Docs/Genomics_lab/assignments/demultiplexing-jordan2lee/2_code/test_data/i1_testfile.fq"
#index2 = "/mnt/c/Users/Jordan/UO_Docs/Genomics_lab/assignments/demultiplexing-jordan2lee/2_code/test_data/i2_testfile.fq"

#num_reads = 10
############## Functions ###################
########################################

# Create reverse compliment function
def reverse_complement(barcode):
    '''Takes string of seq and output string of reverse complementary'''
    comp_bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    return ''.join([comp_bases[base] for base in barcode[::-1]])


# Create function to convert scores to phred scores
def convert_phred(letter):
    """Converts a single character into a phred score based on Illumina 1.8+Phred+33"""
    ASCII = ord(letter) -33
    return(ASCII)


# Create dictionary that has all 24 expected barcoes
#index_dict: keys = barcodes   values= counts (start at 0)
index_dict = {}
b_init = ["GTAGCGTA", "CGATCGAT", "GATCAAGG", "AACAGCGA", "TAGCCATG", "CGGTAATC", "CTCTGGAT", "TACCGGAT", \
    "CTAGCTCA", "CACTTCAC", "GCTACTCT", "ACGATCAG", "TATGGCAC", "TGTTCCGT", "GTCCTAAG", "TCGACAAG", \
    "TCTTCGAC", "ATCATGCG", "ATCGTGGT", "TCGAGAGT", "TCGGATTC", "GATCTTGC", "AGAGTCCA", "AGGATAGC"]
ct_init = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
for i in range(len(b_init)):
    index_dict[b_init[i]]= ct_init[i]

	
# Create function to see if barcode in index_dict
def check_barcode(index_seq):
    '''Take string of index seq and see if in index_dict, returns boolean + input index_seq (True if index seq matches expected - returns False if index seq not expected)'''
    for key in index_dict:
        #print(key)
        if index_seq == key:
            #if input is a index in dict then return true + exit function
            return index_seq == key
        else:
            #if don't find input seq in dict for that index of dict then pass 
            #and continue through function by repeating for next index in dict 
            pass
    #returns False only when else statment is true -- because will exit function on "if index_seq == key:"     
    return False
 
 
 # Create function to see if any bp has a qual score below qual score cutoff
def qs_cutoff (quality_seq):
	'''Take string of a full quality seq from fq/fa + convert to phred score, THEN return True if phred score >30 else false'''
	pos = 0  # start with 1st qual in seq
	#stop after reach end of quality seq
	for pos in range(len(quality_seq)): 
		score = convert_phred(quality_seq[pos])
        
        #if phred score is = or < quality score cuttoff of 30f
		if score <= 30:  
            #return boolean False if any quality score is below or equal quality score cutoff
			return(score>30)
        
		#if phred score is greater than quality score cutoff of 30
		else:   
			continue
			pos+=1
    
    #return boolean True if all quality scores are above quality score cutoff
	return(score>30)

	
########### Demultiplex: body of Code ############################
###################################################

LN = 0

# Initialize counters to 0
total_N = 0
total_unexpected = 0
total_low_qual = 0
total_index_hop = 0
total_correct = 0

#create empty dictionary where key= index seq and values =file handle
fh1_dict = {}
fh2_dict = {}

# OPEN all files
with open(read1, "r") as r1, open(read2, "r") as r2, open(index1, "r") as i1, open(index2, "r") as i2, \
        open("undetermined_r1", "w") as und_r1, open("undetermined_r2", "w") as und_r2:
    
    # Opens file handles for that specific index seq
    for index in b_init: 
        good_r1 = open(index + "_r1.fq", "w")
        good_r2 = open(index + "_r2.fq", "w")
        #create dictionary that has key=index seq and vlaues=file handle (dow this twice)
        fh1_dict[index] = good_r1
        fh2_dict[index] = good_r2

    #tells program to read fastq files 4 lines at a time (header, seq, plus, qual) for ALL input files 
    while i1 and i2 and r1 and r2: # open while have  info 
        
        # read i1 file 4 lines at a time
        header_i1 = i1.readline().strip()
        if not header_i1: #tells to exit loop when reach last line in file
            break
        seq_i1 = i1.readline().strip()
        plus_i1 = i1.readline().strip()
        qual_i1 = i1.readline().strip()
        #lines_i1 = header_i1 + "\n" + seq_i1  + "\n" + plus_i1 + "\n"  + qual_i1  #only for debugging
        #print(lines_i1)

        # read i2 file 4 lines at a time
        header_i2 = i2.readline().strip()
        seq_i2 = i2.readline().strip()
        plus_i2 = i2.readline().strip()
        qual_i2 = i2.readline().strip()
        #lines_i2 = header_i2 + "\n" + seq_i2  + "\n" + plus_i2 + "\n"  + qual_i2  #only for debugging
        #print(lines_i2)

        # read r1 file 4 lines at a time
        header_r1 = r1.readline().strip()
        seq_r1 = r1.readline().strip()
        plus_r1 = r1.readline().strip()
        qual_r1 = r1.readline().strip()
        # Create lines obj so can call all 4 lines (aka all info of one read)
        lines_r1 = header_r1 + ":" + seq_i1 + "\n" + seq_r1  + "\n" + plus_r1 + "\n"  + qual_r1
        #print(lines_r1)
        
        # read r2 file 4 lines at a time
        header_r2 = r2.readline().strip()
        seq_r2 = r2.readline().strip()
        plus_r2 = r2.readline().strip()
        qual_r2 = r2.readline().strip()
        lines_r2 = header_r2 + ":" + seq_i2 + "\n" + seq_r2  + "\n" + plus_r2 + "\n"  + qual_r2
        #print(lines_r2)


        # If either index seq has any N --> write to undeterminded file + add count to counter   
        if "N" in seq_i1 or "N" in seq_i2:  # or reads as "and/or"
            
            # Test
            #print("\nReads that contain N\n", lines_i1, "\n", lines_i2)
            
            # add all 4 lines of r1 read into undetermined_r1 file
            und_r1.write(lines_r1 + "\n") 
            total_N +=1  #add to counter for i1  
            
            # add all 4 lines of r2 read into undetermined_r2 file
            und_r2.write(lines_r2 + "\n")   
            total_N+=1  #add to counter for i2 
            
                    
        # if either index seq is not part of the 24 expected index seq --> write to undetermined file + add count to counter
        elif check_barcode(seq_i1) == False or check_barcode(reverse_complement(seq_i2)) == False:

            # Test
            #print("\n**Index unexpected AND new --> add to und_files**\n", lines_i1, "\n", lines_i2)

            # add all 4 lines of r1 read into undetermined_r1 file
            und_r1.write(lines_r1 + "\n")
            total_unexpected +=1 # add to counter for i1

            # add all 4 lines of r2 to undtermined_r2 file
            und_r2.write(lines_r2 + "\n")
            total_unexpected +=1 # add to counter for i2
                

        # If any either index seq has a bp phred quality score is < quality score cutoff (30) -->  write to undetermined file + add to counter
        elif qs_cutoff(qual_i1) == False or qs_cutoff(qual_i2) == False:

            # add all 4 lines of r1 read into undetermined_r1 file
            und_r1.write(lines_r1 + "\n") 
            total_low_qual +=1  #add to counter for i1

            # add all 4 lines of r2 read into undetermined_r2 file
            und_r2.write(lines_r2 + "\n") 
            total_low_qual+=1  #add to counter for i2
                    
            
            
        # If index 1 and 2 are not reverse complementary to each other --> write to undeter file + add count to counter
        elif reverse_complement(seq_i2) != seq_i1 and "N" not in seq_i1 and "N" not in seq_i2:

            # Test
            #print("\n**Index hop and new --> add to und_files**\n", lines_i1, "\n", lines_i2)

            # add all 4 lines of r1 read into undetermined_r1 file
            und_r1.write(lines_r1 + "\n") 
            total_index_hop +=1  #add to counter for i1

            # add all 4 lines of r2 read into undetermined_r2 file
            und_r2.write(lines_r2 + "\n") 
            total_index_hop+=1  #add to counter for i2


        # If index 1 and 2 are reverse complementary to each other --> write to "good" file
        elif reverse_complement(seq_i2) == seq_i1:

            # test
            #print("\n**Reads that are good**\n", lines_i1, "\n", lines_i2)

			
            #grab value from dict -- this is the fh
            good_r1 = fh1_dict[seq_i1]
            # add all 4 lines of r1 read into undetermined_r1 file
			good_r1.write(lines_r1 + "\n") 
            total_correct +=1  #add to counter for i1

            #create file with <reverse complement of barcode>_r2.fq
            #so that matches with the forward read barcode (e.g. good_r1)
            good_r2 = fh2_dict[reverse_complement(seq_i2)]
            # add all 4 lines of r2 read into undetermined_r2 file
			good_r2.write(lines_r2 + "\n") 
            total_correct+=1  #add to counter for i2 
        
    #close file handles for that specific index seq
    for index in b_init:
        fh1_dict[index].close()
        fh2_dict[index].close()
                    
    
    #Create a useful report for the end user of your code
    str_N =str(total_N)
    str_N_prop = str(total_N/num_reads)
    str_low_qual = str(total_low_qual)
    str_low_qual_prop = str(total_low_qual/num_reads)
    str_hop = str(total_index_hop)
    str_hop_prop = str(total_index_hop/num_reads)
    str_good = str(total_correct)
    str_good_prop = str(total_correct/num_reads)
    str_une = str(total_unexpected)
    str_une_prop = str(total_unexpected/num_reads)
    
    with open("report_demultiplex", "w") as report:
        report.write("Total reads passed (good): " + str_good + "\n" + "Proportion reads passed (good): " + str_good_prop + "\n")
        report.write("Total reads with N: " + str_N + "\n" + "Proportion reads with N: " + str_N_prop + "\n")
        report.write("Total reads with low quality: " + str_low_qual + "\n" + "Proportion reads with low quality: " + str_low_qual_prop + "\n")
        report.write("Total reads with index hopping: " + str_hop + "\n" + "Proportion reads with index hopping: " + str_hop_prop + "\n")
        report.write("Total reads without N and not matching expected barcodes: " + str_une + "\n" + "Proportion reads without N and not matching expected barcodes: " + str_une_prop)
        