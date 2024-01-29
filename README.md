import pysam
import pandas
import numpy

input_bam = "/Users/meggielam/Desktop/test/SC-1_1_SC-1_2.hg38plusAkata_inverted-Aligned.sortedByCoord.out.bam"



def get_var(path_to_bamfile, variant_chromosome, variant_position):
    
    bam = pysam.AlignmentFile(path_to_bamfile, "rb")
    reads = bam.fetch(variant_chromosome, variant_position, variant_position + 1)
    base_counts = {'A':0, 'C':0, "G":0, "T":0}

    for read in reads:
        difference = variant_position - read.pos
        seq = read.seq
        if difference < len(seq):
            base = seq[difference-1]
            base_counts[base] += 1
    return base_counts

df = pd.read_csv("/Users/meggielam/Desktop/test/filtered_final_2.csv", usecols = ["Chrom", "Start"])
rows = len(df.index)
columns = 4 
m = np.zeros([rows,columns])

i = 0 

for chrom, start in zip(df['Chrom'], df['Start']):
    base_counts = get_var(bam, chrom, start)
    A = base_counts['A']
    C = base_counts['C']
    G = base_counts['G']
    T = base_counts['T']
    m[i, 0] = A 
    m[i, 1] = C
    m[i, 2] = G
    m[i, 3] = T

    i = i + 1
    if i == 1000:
        break 

#this will tell me what is at row 50 -- m[50, :]
pd.DataFrame(m, columns=['A','C','G', 'T'])
