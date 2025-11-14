# ZIP_CRISPR_code
Repository for code that was used in the publication "Improving efficiency of homology-directed repair with ZIP CRISPR"
 
---

## 📂 Repository structure
```
├── Barcode
│ └── barcode01.sorted.bam
│ └── barcode02.sorted.bam
│ └──...
├── bam2Excel.py
├── model.xlxs
├── tribam.py
├── sequences.txt


```


## ⚙️ Analyze

- Python ≥ 3.7  
- [pysam](https://pysam.readthedocs.io/en/latest/)  
- [pandas](https://pandas.pydata.org/)
  
Tested with:  
pysam 0.22.0   
pandas 2.3.1  
Python 3.13.2  

1. **bam2Excel.py**
  
This script processes .sorted.bam files and calculates the proportion of each nucleotide for each selected position to determine the HDR editing rate per base, returning these data in xlsx format.

Run from project root:
```bash
python3 bam2excel.py
```

2. **tribam.py**  
  
This script searches specific sequences in the BAM corresponding to unedited, HDR edited sequences or eventually predominant InDels seen with IGV (Integrative Genomics Viewer) visualization to determine the HDR editing rate per read.

Run from project root:

```bash
python3 tribam.py
```

Output:
Excel file with proportion of CRISPR outome scenarios in % + sub files with BAMs. 




Date of last update: 13 november 2025

📄 License

This project is licensed under the Creative Commons CC0 1.0 Universal license.
You are free to copy, modify, distribute and use the code, even for commercial purposes, without asking permission.

Citation or acknowledgment of the original author is appreciated.

