#!/usr/bin/env python3
import argparse
import gzip

def get_args(): # Usage for user input
	parser = argparse.ArgumentParser(description='Demultiplex sequencing data based on respective barcodes')
	parser.add_argument("-r1", "--read1", help = "This needs to be a FASTQ file of forward sequence reads", required = True)
	parser.add_argument("-r2", "--read2", help = "This needs to be a FASTQ file of forward barcode reads", required = True)
	parser.add_argument("-r3", "--read3", help = "This needs to be a FASTQ file of reverse barcode reads", required = True)
	parser.add_argument("-r4", "--read4", help = "This needs to be a FASTQ file of reverse sequence reads", required = True)
	parser.add_argument("-i", "--index", help = "This needs to be a .tsv file of all known indexes present in the sequencing outputs", required = True)
	return parser.parse_args()
args = get_args()

read1_file = args.read1 # read in user input and save as a file record
read2_file = args.read2 # read in user input and save as a file record
read3_file = args.read3 # read in user input and save as a file record
read4_file = args.read4 # read in user input and save as a file record
Index_file = args.index # read in user input and save as a file record

Read1_file = gzip.open(read1_file, 'rt') #open file handle for associated read file
Index1_file = gzip.open(read2_file, 'rt') #open file handle for associated read file
Index2_file = gzip.open(read3_file, 'rt') #open file handle for associated read file
Read2_file = gzip.open(read4_file, 'rt') #open file handle for associated read file

def convert_phred(qual_sequence):
	'''Converts a single character into a phred score'''
	mean_scores = [] # initialize empty dictionary to append scores to
	for i in range(len(qual_sequence)): # for each letter in the quality sequence
	    mean_scores.append(0.0) #initialize an empty array with the length of the sequence
	index = 0 # counter to keep track of which index the following loop is on
	for letter in qual_sequence: # itterate throug each letter of the quality sequence
		phred_score = ord(letter)-33 # convert the letter to a numerical score
		mean_scores[index] += phred_score # insert the score into the mean score list for the appropriate position
		index+=1 # incriment counter
	total = 0 # initialize summation of all known scores
	for value in mean_scores: # iterate through each score in the list of quality scores
		total += value # add to the total score
	average_score = total/len(qual_sequence) # take the average score for a quality sequence
	return average_score

def rev_comp(DNA_seq):
	''' Take DNA string and return its reverse compliment)'''
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} # dictionary of bases as keys and their complement as values
	sequence_list = list(DNA_seq) # make the DNA string into a list object
	complement_list = [] # initialize an empty list to append complements to
	for pos in sequence_list: # iterate over the sequence list
		sequence_list = complement_list.append(complement[pos]) # append to the complement list the complement of each position in sequence list
	comp = ''.join(complement_list) # convert complement list into a string
	DNA_rev_comp = comp[::-1] # reverse the complement
	return DNA_rev_comp

def get_indexes(Known_Index_File):
	''' Take the entire .tsv file of known indexes and put all indexes  and their reverse complements into a dictionary'''
	file_of_indexes  = open(Known_Index_File, "r") # open the index file with a handle
	index_dictionary = {} # initialize an empty dictionary to store indexes as keys and their reverse complement as values
	for line in file_of_indexes: # iterate over the lines in the index file
		line = line.strip() # strip the current line of trailing white space
		sample, barcode = line.split() # split the line into columns
		index_dictionary[barcode] = rev_comp(barcode) # put the barcode into the dictionary as a key and its rev comp as a value
		index_dictionary[rev_comp(barcode)] = barcode # put the rev comp into the dictionary as a key and the barcode as a value
	file_of_indexes.close() # close the file being read
	return index_dictionary

def fastq_records_iterator(file_handle): # from help with Jared Galloway
    ''' Take in a file handle and return a generator (iterator) object which yeilds four \n stripped lines of a file as a list of strings. '''
    while True: # always do this for the file given
        record = [] # initialize an empty record
        for i in range(4): # iterate over the lines in the file
            record.append(file_handle.readline().strip()) # strip the line of trailing whitespace, and append it to the record
        if record[0] == '': # check if the first item in the record is empty (this will happen once the file is read completely)
            break # end the creation of records
        else: # while there are records to generate 4 lines at a time
            yield record # return a generator object to iterate over containing 4 lines of a complete fastq record
    return None

index_dict = get_indexes(Index_file) # create a dictionary of indexes and their reverse complements as keys containing the rev comp as values

total_counter = 0 # keep a count of all the records written out to files
total_hopped_counter = 0 # keep a count of only all the hopped records written out to files
hopped_counter = {} # initialize a dictionary to count the number of each index that is encountered (both forward and reverse)
handle_dict = {} # create a dictionary to store unique file names for each key file handle (for unique write files)
barcode_counter = {} # initialize a dictionary to count the total number of barcodes encountered (to report number of records in each codition)
encountered_counter = 0


# read in four lines at a time from all files at once
for record in zip(fastq_records_iterator(Read1_file), fastq_records_iterator(Index1_file), fastq_records_iterator(Index2_file), fastq_records_iterator(Read2_file)):
    encountered_counter += 1
    Read1,Index1,Index2,Read2 = record # break up the record list into its sublists
    index1_quality = convert_phred(Index1[3]) # get average wuality of the Index 1 sequence
    index2_quality = convert_phred(Index2[3]) # get the average quality of the Index 2 sequence
    if (Index1[1] not in index_dict) or (Index2[1] not in index_dict): # check the indexes doesn't match any dictionary key and (already doesnt contain 'N') (LOW QUALITY)
        if "undetermined" not in barcode_counter: # if the query is not in the dict, put it in there
            barcode_counter['undetermined'] = 1 # # create an undetermined counter if it doesnt exist
            R1_undetermined_out = open(f"R1_undetermined.fastq", "w") # open out write files for undetermined case
            R2_undetermined_out = open(f"R2_undetermined.fastq", "w") # open out write files for undetermined case
            R1_undetermined_out.write(f"{Read1[0]} {Index1[1]}-{Index2[1]}\n{Read1[1]}\n{Read1[2]}\n{Read1[3]}\n") # write complete R1 record out
            R2_undetermined_out.write(f"{Read2[0]} {Index1[1]}-{Index2[1]}\n{Read2[1]}\n{Read2[2]}\n{Read2[3]}\n") # write complete R2 record out
            total_counter += 1 # increment total counter
        else:
            barcode_counter['undetermined'] += 1 # Increment unknown counter in dictionary
            R1_undetermined_out.write(f"{Read1[0]} {Index1[1]}-{Index2[1]}\n{Read1[1]}\n{Read1[2]}\n{Read1[3]}\n") # write complete R1 record out
            R2_undetermined_out.write(f"{Read2[0]} {Index1[1]}-{Index2[1]}\n{Read2[1]}\n{Read2[2]}\n{Read2[3]}\n") # write complete R2 record out
            total_counter += 1 # increment total counter
    elif index1_quality<30 or index2_quality<30: # check quality of an index is below a determnined threshold (LOW QUALITY)
        if "undetermined" not in barcode_counter: # if the query is not in the dict, put it in there
            barcode_counter['undetermined'] = 1 # create an undetermined counter if it doesnt exist
            R1_undetermined_out = open(f"R1_undetermined.fastq", "w") # open out write files for undetermined case
            R2_undetermined_out = open(f"R2_undetermined.fastq", "w") # open out write files for undetermined case
            R1_undetermined_out.write(f"{Read1[0]} {Index1[1]}-{Index2[1]}\n{Read1[1]}\n{Read1[2]}\n{Read1[3]}\n") # write complete R1 record out
            R2_undetermined_out.write(f"{Read2[0]} {Index1[1]}-{Index2[1]}\n{Read2[1]}\n{Read2[2]}\n{Read2[3]}\n") # write complete R2 record out
            total_counter += 1 # increment total counter
        else: # if undetermined already exists in the dictionary keeping track of conditions
            barcode_counter['undetermined'] += 1 # incriment counter of undetermined records
            R1_undetermined_out.write(f"{Read1[0]} {Index1[1]}-{Index2[1]}\n{Read1[1]}\n{Read1[2]}\n{Read1[3]}\n") # write complete R1 record out
            R2_undetermined_out.write(f"{Read2[0]} {Index1[1]}-{Index2[1]}\n{Read2[1]}\n{Read2[2]}\n{Read2[3]}\n") # write complete R2 record out
            total_counter += 1 # increment total counter
    elif Index1[1] in index_dict: # checks if the index read from index 1 is in the index dictionary
        if Index2[1] == index_dict[Index1[1]]: # check two indexes match a key,value pair in the dictionary (MATCHED)
            if f"{Index1[1]}-{Index2[1]}" not in barcode_counter: # if the unique barcode pair query is not in the dict, put it in there
                barcode_counter[f"{Index1[1]}-{Index2[1]}"] = 1 # create a unique barcode pair counter if it doesnt exist
                handle_dict[f"R1_{Index1[1]}_{Index2[1]}"] = f"R1_{Index1[1]}-{Index2[1]}_out" # store a unique file handle in a dictionary for R1 unique barcode pair
                handle_dict[f"R2_{Index1[1]}_{Index2[1]}"] = f"R2_{Index1[1]}-{Index2[1]}_out" # store a unique file handle in a dictionary for R2 unique barcode pair
                handle_dict[f"R1_{Index1[1]}_{Index2[1]}"] = open(f"R1_{Index1[1]}-{Index2[1]}.fastq", "w") # open a unique R1 barcode pair out file
                handle_dict[f"R2_{Index1[1]}_{Index2[1]}"] = open(f"R2_{Index1[1]}-{Index2[1]}.fastq", "w") # open a unique R1 barcode pair out file
                handle_dict[f"R1_{Index1[1]}_{Index2[1]}"].write(f"{Read1[0]} {Index1[1]}-{Index2[1]}\n{Read1[1]}\n{Read1[2]}\n{Read1[3]}\n") # write complete R1 record out
                handle_dict[f"R2_{Index1[1]}_{Index2[1]}"].write(f"{Read2[0]} {Index1[1]}-{Index2[1]}\n{Read2[1]}\n{Read2[2]}\n{Read2[3]}\n") # write complete R2 record out
                total_counter += 1 # increment total counter
            else: # if a unique barcode-pair already exists in the dictionary keeping track of conditions
                barcode_counter[f"{Index1[1]}-{Index2[1]}"] += 1 # incriment counter of the unique barcode-pair records
                handle_dict[f"R1_{Index1[1]}_{Index2[1]}"].write(f"{Read1[0]} {Index1[1]}-{Index2[1]}\n{Read1[1]}\n{Read1[2]}\n{Read1[3]}\n") # write complete R1 record out
                handle_dict[f"R2_{Index1[1]}_{Index2[1]}"].write(f"{Read2[0]} {Index1[1]}-{Index2[1]}\n{Read2[1]}\n{Read2[2]}\n{Read2[3]}\n") # write complete R2 record out
                total_counter += 1 # increment total counter
        elif (Index2[1] != index_dict[Index1[1]]) and (Index2[1] in index_dict): # both of the indexes matches a key in the dictionary and they are not each others direct dictionary pairs (INDEX HOPPED)
            if "index_hopped" not in barcode_counter: # if the index hopped query is not in the dict, put it in there
                barcode_counter["index_hopped"] = 1 # create an index hopped counter if it doesnt exist
                R1_hopped_out = open(f"R1_index_hopped.fastq", "w") # open out write files for hopped case
                R2_hopped_out = open(f"R2_index_hopped.fastq", "w") # open out write files for hopped case
                R1_hopped_out.write(f"{Read1[0]} {Index1[1]}-{Index2[1]}\n{Read1[1]}\n{Read1[2]}\n{Read1[3]}\n") # write complete R1 record out
                R2_hopped_out.write(f"{Read2[0]} {Index1[1]}-{Index2[1]}\n{Read2[1]}\n{Read2[2]}\n{Read2[3]}\n") # write complete R2 record out
                if Index1[1] not in hopped_counter: # check if index 1 is in the hopped counter dictionary
                    hopped_counter[Index1[1]] = 1 # create it if its not there
                if Index2[1] not in hopped_counter: # also check if index 2 is in the hopped counter dictionary
                    hopped_counter[Index2[1]] = 1 # create it if its not there
                total_counter += 1 # increment total counter
                total_hopped_counter += 1 # increment total hopped counter
            else: # if index hopped is already a key in the barcode dictionary
                barcode_counter["index_hopped"] += 1 # increment the index hopped key
                R1_hopped_out.write(f"{Read1[0]} {Index1[1]}-{Index2[1]}\n{Read1[1]}\n{Read1[2]}\n{Read1[3]}\n") # write complete R1 record out
                R2_hopped_out.write(f"{Read2[0]} {Index1[1]}-{Index2[1]}\n{Read2[1]}\n{Read2[2]}\n{Read2[3]}\n") # write complete R2 record out
                if Index1[1] not in hopped_counter: # check if index 1 is in the hopped counter dictionary
                    hopped_counter[Index1[1]] = 1 # create it if its not there
                else: # if the current index1 is in the dictionary of hopped indexes already
                    hopped_counter[Index1[1]] += 1 # increment present index 1
                if Index2[1] not in hopped_counter: # also check if index 2 is in the hopped counter dictionary
                    hopped_counter[Index2[1]] = 1 # create it if its not there
                else: # if the current index2 is in the dictionary of hopped indexes already
                    hopped_counter[Index2[1]] += 1 # increment present index 1
                total_counter += 1 # increment total counter
                total_hopped_counter +=1 # increment only hopped index counter

Index1_file.close() # close all files being read
Read1_file.close() # close all files being read
Read2_file.close() # close all files being read
Index2_file.close() # close all files being read

out_record = open("demultiplex_record.out", "w") # open a file to write records out to

out_record.write(f"Files Processed: \n {read1_file},\n {read2_file},\n {read3_file},\n {read4_file},\n {Index_file}\n\n") # Write out records in a report
out_record.write(f"Total number of files written: {len(barcode_counter)*2}\n")
out_record.write(f"Total number of records written: {total_counter}\n")
out_record.write(f"Total number of records encountered: {encountered_counter}\n\n")
for key in barcode_counter:
	out_record.write(f"{key} read-pair files contain {barcode_counter[key]} records\n ({float((barcode_counter[key]/total_counter)*100)}% of total records)\n")
out_record.write("\nOverall amount of index swapping: \n")
out_record.write(f"Barcode\tObservations\tPercent of total hopped indexes\n")
for key in hopped_counter:
	value = hopped_counter[key] # associate the number of index hopped with the index
	out_record.write(f"{key}\t{value}\t{(value/total_hopped_counter)*100}\n") # write out average observation for each index hopped
out_record.close()
