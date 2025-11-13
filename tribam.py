import pysam
from difflib import SequenceMatcher
import matplotlib.pyplot as plt
import os
import pandas as pd
from openpyxl import Workbook
from openpyxl.drawing.image import Image
from openpyxl.utils.dataframe import dataframe_to_rows

# To read bam file
def read_bam(bam_file):
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    reads = [read for read in bamfile.fetch()]
    bamfile.close()
    return reads

# To read sequence.txt file
def read_target_sequences(sequences_file):
    target_sequences = {}
    with open(sequences_file, 'r') as file:
        for line in file:
            title, sequence = line.strip().split('\t')
            target_sequences[sequence] = title
    return target_sequences

# To construct read sequence including CIGAR strings
def rebuild_sequences(read):
    sequence = read.query_sequence
    if sequence is None:
        return None

    cigar = read.cigar
    expanded_sequence = ""
    seq_pos = 0
    
    for operation, length in cigar:
        if operation == 0:  # same
            expanded_sequence += sequence[seq_pos:seq_pos + length]
            seq_pos += length
        elif operation == 1:  # Insertion (I)
            expanded_sequence += sequence[seq_pos:seq_pos + length]
            seq_pos += length
        elif operation == 2:  # Deletion (D)
            expanded_sequence += '-' * length
        elif operation == 4:  # Soft clip (S)
            seq_pos += length
        elif operation == 5:  # Hard clip (H)
            pass
    
    return expanded_sequence

# To determine reverse complement
def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(sequence))

# To check read quality
def satisfactory_quality_reads(read, quality_threshold=40):
    return read.mapping_quality >= quality_threshold

# To find the best match
def best_match(read_sequence, target_sequences, tolerance=0):
    if read_sequence is None:
        return None
    
    for sequence, title in target_sequences.items():
        for seq in [sequence, reverse_complement(sequence)]:
            if seq in read_sequence:
                return sequence
            # For partial match
            for i in range(len(read_sequence) - len(seq) + 1):
                sub_sequence = read_sequence[i:i+len(seq)]
                ratio = SequenceMatcher(None, sub_sequence, seq).ratio()
                if ratio >= (1 - tolerance / len(seq)):
                    return sequence
    return None

# To sort the reads
def read_sorting(reads, target_sequences, quality_threshold=40):
    groups = {title: [] for title in target_sequences.values()}
    others = []
    bad_quality = []
    
    for read in reads:
        if not satisfactory_quality_reads(read, quality_threshold):
            bad_quality.append(read)
            continue
        
        sequence = rebuild_sequences(read)
        match = best_match(sequence, target_sequences)
        if match:
            title = target_sequences[match]
            groups[title].append(read)
        else:
            others.append(read)
    
    return groups, others, bad_quality

# To add reads in bam files
def bam_writing(groups, others, bad_quality, output_folder, bamfile_template):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for groupe, reads in groups.items():
        bam_file = os.path.join(output_folder, f"{groupe}.bam")
        with pysam.AlignmentFile(bam_file, "wb", template=bamfile_template) as outfile:
            for read in reads:
                outfile.write(read)
        # To index bam file
        pysam.index(bam_file)
    
    bam_file_others = os.path.join(output_folder, "others.bam")
    with pysam.AlignmentFile(bam_file_others, "wb", template=bamfile_template) as outfile:
        for read in others:
            outfile.write(read)
    # To index bam file with other reads
    pysam.index(bam_file_others)
    
    bam_file_bad_quality = os.path.join(output_folder, "bad_quality.bam")
    with pysam.AlignmentFile(bam_file_bad_quality, "wb", template=bamfile_template) as outfile:
        for read in bad_quality:
            outfile.write(read)
    # To index bam file with low quality reads
    pysam.index(bam_file_bad_quality)

# To generate a chart
def generate_distribution_chart(groups, output_folder, file_name):
    labels = list(groups.keys())
    values = [len(reads) for reads in groups.values()]

    # To define chart colors
    colors = ['#aec6cf', '#cfc9e9', '#f4c2c2', '#f49ac2', '#c1c1ff', '#d6cadd', '#ffb3de']
    colors = colors[:len(labels)]  # To adjust to sequence number
    
    plt.figure(figsize=(10, 6))
    plt.pie(values, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
    plt.axis('equal')
    plt.title('read_distribution')
    
    # To save chart
    chart_file = os.path.join(output_folder, f"{file_name}.png")
    plt.savefig(chart_file)
    plt.close()

    return chart_file

# To deal with bam file
def bam_file_treatment(bam_file, target_sequences, output_folder):
    # To read bam files
    reads = read_bam(bam_file)
    
    # To sort reads
    groups, others, bad_quality = read_sorting(reads, target_sequences)
    
    # To whrite sorted reads in bam files and index files
    with pysam.AlignmentFile(bam_file, "rb") as bamfile_template:
        bam_writing(groups, others, bad_quality, output_folder, bamfile_template)
    
    # To create a chart and save it
    file_name = os.path.splitext(os.path.basename(bam_file))[0]
    chart_file = generate_distribution_chart(groups, output_folder, file_name)
    
    # To create metrics
    metrics = {
        "Total reads": len(reads),
        "bad_quality_reads": len(bad_quality),
        "reads in groups": {groupe: len(reads) for groupe, reads in groups.items()},
        "others reads": len(others)
    }
    
    return metrics, chart_file

# To create excel file with the results
def excel_file_creation(results, excel_file):
    wb = Workbook()
    ws_summary = wb.active
    ws_summary.title = "abstract"
    
    for i, (file_name, (metrics, chart_file)) in enumerate(results.items()):
        ws = wb.create_sheet(title=file_name)
        
        # To add metrics to the sheet
        df_metrics = pd.DataFrame.from_dict(metrics["reads in groups"], orient='index', columns=['number of reads'])
        df_metrics.loc["Total reads"] = metrics["Total reads"]
        df_metrics.loc["bad_quality_reads"] = metrics["bad_quality_reads"]
        df_metrics.loc["others reads"] = metrics["others reads"]
        
        for r in dataframe_to_rows(df_metrics, index=True, header=True):
            ws.append(r)
        
        # To add the chart
        img = Image(chart_file)
        ws.add_image(img, 'E5')
        
        # To add data to the summary
        ws_summary.append([file_name, metrics["Total reads"], metrics["bad_quality_reads"], metrics["others reads"]])
        for groupe, count in metrics["reads in groups"].items():
            ws_summary.append([f"{file_name} - {groupe}", count])
    
    wb.save(excel_file)

# Main function to deal with bam files in a folder
def main(bam_folder, sequences_file, output_folder, excel_file):
    target_sequences = read_target_sequences(sequences_file)
    
    results = {}
    for bam_file in os.listdir(bam_folder):
        if bam_file.endswith(".bam"):
            bam_file_path = os.path.join(bam_folder, bam_file)
            metrics, chart_file = bam_file_treatment(bam_file_path, target_sequences, output_folder)
            results[bam_file] = (metrics, chart_file)
    
    excel_file_creation(results, excel_file)

# To use main function with input and output files
if __name__ == "__main__":
    main("bam_folder", "sequences.txt", "sortieTriBam", "results_Tribam")
    