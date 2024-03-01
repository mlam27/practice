import pysam
import pandas as pd
import numpy as np

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

df = pd.read_csv("/Users/meggielam/Desktop/test/filtered_final_2.csv")
df = df[df['Chrom']!= 'chrM']
df.index=range(len(df.index))

rows = len(df.index)
columns = 4 
m = np.zeros([rows,columns])

i = 0 

for chrom, start, alt,cline in zip(df['Chrom'], df['Start'], df['Alt'],df['CellLineName']):
    base_counts = get_var(input_bam, chrom, start)
    called_base = max(base_counts, key=base_counts.get)
    depth = sum(val for val in base_counts.values())
    max_base = base_counts[called_base]
    if depth < 100:
        continue
    elif max_base / depth < .8:
        continue
   
    elif alt == called_base :
        print(cline)

    A = base_counts['A']
    C = base_counts['C']
    G = base_counts['G']
    T = base_counts['T']
    m[i, 0] = A 
    m[i, 1] = C
    m[i, 2] = G
    m[i, 3] = T

    i = i + 1


minimum_reads = 30
md = pd.DataFrame(m, columns=['A','C','G', 'T']) #this will tell me what is at row 50 -- m[50, :]
md["Ref"] = list(df.Ref) #for each row want to see what the base was suppose to be before the mutation
md['called_nuc'] = [md.columns[i] for i in np.argmax(m,1)]  # Get most freq base
md['CellLineName']=df['CellLineName']
md['expected_alt'] = df['Alt']

md = md[np.sum(md[md.columns[:4]],1) > minimum_reads]
t_f = []
for ref, called_nuc in zip(md['Ref'], md['called_nuc']):
    t_f.append(ref == called_nuc)
md['match']=t_f
md['alt_match']=md['called_nuc'] == md['expected_alt']
md[(~md['match']) & (md['alt_match'])& (md['Ref'].str.len() ==1)]['CellLineName'].value_counts()







