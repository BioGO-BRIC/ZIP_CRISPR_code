import pysam
import pandas as pd
import os

# List all the barcodes to analyze
bam_files = [
    './Barcode/barcode01.sorted.bam',
    './Barcode/barcode02.sorted.bam',
    ]


# Input and output files
excel_file_path = './model.xlsx'
output_excel_path = './resultsBam2Excel.xlsx'

# Charge excel file
df_excel = pd.read_excel(excel_file_path)  # Add skiprows if necessary


print(df_excel.columns)


# Positions to extract
start_position = 125815045 #genomic position (for instance 125815045)
end_position = 125815075 #genomic position (for instance 125815075)
reference_name = 'chr10'  # Chromosom number (for instance 'chr10')

# To prepare the excel ouput file
with pd.ExcelWriter(output_excel_path, engine='openpyxl') as writer:
    for bam_file in bam_files:
        bam_basename = os.path.basename(bam_file)
        sheet_name = bam_basename[:9] 

        # To open bam file
        bam = pysam.AlignmentFile(bam_file, "rb")

        # To prepare a list to collect data
        data = []

        # To extract data with filters simalar to IGV'ones
        for pileupcolumn in bam.pileup(reference_name, start_position-1, end_position, truncate=True, min_base_quality=0, min_mapping_quality=0):
            pos_data = [pileupcolumn.pos + 1]  # +1 for indexing
            base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Del': 0, 'Ins_A': 0, 'Ins_T': 0, 'Ins_C': 0, 'Ins_G': 0}
            delins_details = []

            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    base_counts[base.upper()] += 1
                if pileupread.is_del:
                    base_counts['Del'] += 1
                elif pileupread.indel > 0:
                    inserted_seq = pileupread.alignment.query_sequence[pileupread.query_position + 1:pileupread.query_position + 1 + pileupread.indel]
                    for ins_base in inserted_seq:
                        base_counts[f'Ins_{ins_base.upper()}'] += 1

            pos_data.extend(base_counts.values())
            data.append(pos_data)

        # To create a dataframe with extracted data
        columns = ['Position', 'A', 'C', 'G', 'T', 'Del', 'Ins_A', 'Ins_T', 'Ins_C', 'Ins_G']
        df_pileup = pd.DataFrame(data, columns=columns)

        # To add reversed columns 
        reversed_df = df_pileup.iloc[::-1].reset_index(drop=True)
        reversed_columns = [f'{col}_inv' for col in columns[1:]]  # Creating inverted column names
        reversed_df.columns = ['Position'] + reversed_columns
        
        # To fuse dataframe
        df_merged = pd.merge(df_excel, df_pileup, on='Position', how='left')
        df_merged = pd.concat([df_merged, reversed_df.iloc[:, 1:]], axis=1)  # Horizontal concatenation

        # To save in a new excel sheet
        df_merged.to_excel(writer, sheet_name=sheet_name, index=False)

        # To close bam file
        bam.close()