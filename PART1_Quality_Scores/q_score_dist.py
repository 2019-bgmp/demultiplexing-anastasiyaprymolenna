#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import gzip

def get_args(): # Usage
	parser = argparse.ArgumentParser(description='Parse Fasta file for the longest protein sequence')
	parser.add_argument("-f", "--file", help = "This needs to be a FASTA file of contigs. Sequences must be contained on one line.", required = True)
	parser.add_argument("-p", "--figure_name", help = "Specify the name you want the figure to be")
	return parser.parse_args()
args = get_args()

file = args.file
figure = args.figure_name

leng = 0 #get the length of the line of quality/seq scores
with gzip.open(file, "rt") as lineFile:
	LC = 0
	for line in lineFile:
		line = line.strip()
		LC += 1
		if LC == 2:
			leng = len(line)
			break

def convert_phred(letter):
	"""Converts a single character into a phred score"""
	phred_score = ord(letter)-33
	return phred_score

mean_scores = []
for i in range(0,leng):
    mean_scores.append(0.0) #initialize an empty array

LN = 0 #start counter
with gzip.open(file, 'rt') as fastq:
    for line in fastq:
        line = line.strip()
        if LN%4==3: #only work with quality score lines
            index = 0
            for position in line:
                score = convert_phred(position) #use def function to convert the score to a number
                mean_scores[index] += score
                index+=1
        LN+=1

number_of_reads = LN/4
print(f"# Base Pair\tMean Quality Score")
for i in range(0,leng):
    average = mean_scores[i]/number_of_reads
    mean_scores.pop(i)
    mean_scores.insert(i,average)
    print("{0}\t{1}".format(i, round(average,1)))

x = list(range(0,leng))
y = mean_scores
plt.bar(x,y)
plt.title("Means of Phred Quality Scores at Respective Sequence Positions")
plt.xlabel("Position in Sequence (0-based)")
plt.ylabel("Mean Quality Score")
plt.show()
plt.savefig(f"{figure}_plot.png")
