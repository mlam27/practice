import pandas as pd
#filtering file to keep values >0.4 for column AF
new = pd.read_table("/Users/meggielam/Desktop/test/OmicsSomaticMutations.csv",sep=",")
new ["AF"] #will show what is in column AF
new ["AF"] >0.4 #filters values that are greater than 0.4
new[new ["AF"] > 0.4] #filtered into original file
filtered = new[new ["AF"] > 0.4] #called it 'filtered' in terminal
filtered.to_csv ("filtered_list.csv") #saved list onto computer 

#filtering previous file (filtered_list) and only keeping certain columns
filtered = new[new ["AF"] > 0.4]
filtered_keep = ["Chrom", "Pos", "Ref", "Alt", "ModelID"] #will only keep columns that are called those
filtered = filtered[filtered_keep] 
filtered.to_csv ("filtered_keep_only.csv")

#adding a column to establish "stop" position
filtered = pd.read_csv('/Users/meggielam/Desktop/test/filtered_keep_only.csv')
filtered['Stop'] = filtered['Pos'] + 1
filtered = filtered[['Chrom', 'Pos', 'Stop', 'Ref', 'Alt', 'ModelID']] #reordering the columns
filtered = filtered.rename(columns={'Pos': 'Start'}) #renames 'Pos' to "Start"
filtered.to_csv ("filtered_new.csv")

#merged columns chromo and start and made a new column of result
filtered['Chrom:Start'] = filtered['Chrom'].astype(str) + ':' + filtered['Start'].astype(str)
filtered.to_csv ("filtered_merged.csv")

#looked for duplications of chrom:start
duplicated = filtered[filtered.duplicated(subset='Chrom:Start', keep=False)]
duplicated.to_csv ("filtered_duplicated.csv")

#removed the duplicated rows based off chrom:start
mask_duplicated = filtered.duplicated(subset='Chrom_Start', keep=False) #Create a boolean mask for duplicated rows based on 'Chrom_Start'

duplicated_rows = filtered[mask_duplicated] #Created a DataFrame with the duplicated rows
duplicated_rows.to_csv ("dup_removal_from_filtered.csv")

filtered_no_duplicates = filtered[~mask_duplicated] # Created a DataFrame without the duplicated rows
filtered_no_duplicates.to_csv ("filtered_no_dups.csv")

#matched ModelID to CellLineName from depmap database to my list**
model = pd.read_csv('/Users/meggielam/Desktop/test/Model.csv')
filtered_no_dups = pd.read_csv('/Users/meggielam/Desktop/test/filtered_no_dups.csv')
merged_df = pd.merge(filtered_no_dups, model[['ModelID', 'CellLineName']], how='left', on='ModelID')

#can also use this loop to match ModelID to CellLineName**
for model_id in filtered_no_dups['ModelID'].unique():
    if model_id in model['ModelID'].values:
        model_subset = model[model['ModelID'] == model_id][['ModelID', 'CellLineName']]
# Merge the subset with filtered_no_dups using a left merge
        cell_line = pd.merge(filtered_no_dups[filtered_no_dups['ModelID'] == model_id], model_subset, how='left', on='ModelID')
# Update the 'CellLineName' column in filtered_no_dups
        filtered_no_dups.loc[filtered_no_dups['ModelID'] == model_id, 'CellLineName'] = cell_line['CellLineName_y'].tolist()
        #CellLineName_y: 'CellLineName' column from model_subset