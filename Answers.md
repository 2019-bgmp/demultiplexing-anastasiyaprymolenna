Table of File Labels

Part 1:

| File Name  | Category |
| ------------- | ------------- |
| 1294_S1_L008_R1_001.fastq.gz  | Read 1 |
| 1294_S1_L008_R2_001.fastq.gz  | Index 1 |
| 1294_S1_L008_R3_001.fastq.gz  | Index 2 |
| 1294_S1_L008_R4_001.fastq.gz  | Read 2 |

2b) What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively?

 After a slight existential crisis, I would say that a good quality score for the indexes and read pairs depends on the data that came off the sequencer. In general a score greater than 10 will catch 90% of our data and 20 will catch 99.9%. I would generally lean more towards 20 because I want to get quality data. However, a lot of parameters must be considered when determinging a cut off, sometimes with trial and error. I know that now. It is important to consider the conditions of the sequencer and whether an entire lane could be messed up or just a couple reads could be bad that can be filtered. Essentially, we don't want to throw away the baby with the bathwater if all of our reads are bad but it was on the fault of the sequencer. Then again there are a lot of things that could go wrong so we need to consider case by case each filter we do and how that impacts out downstream assembly or analysis and try maybe a couple cut offs that gives us the most reproducible data in the end. Reproducibility is always the goal.

2c)How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command you used. CHALLENGE: use a one-line command

 Command: $ zcat /path/file | awk 'NR%4==2' | grep 'N' | wc -l
 
 Index 1 undetermined count = 3976613
 
 Index 2 undetermined count = 3328051
