"""Generate real-world test datasets for BioAgent."""
import pandas as pd
import numpy as np
from pathlib import Path

Path("data/sample").mkdir(parents=True, exist_ok=True)

# --- Dataset 1: Breast Cancer RNA-seq ---
genes = [
    'ESR1','ERBB2','MKI67','BRCA1','BRCA2','TP53','PIK3CA',
    'PTEN','CDH1','GATA3','FOXA1','AR','PGR','EGFR','MYC',
    'CCND1','RB1','ATM','CHEK2','PALB2','RAD51','AURKA',
    'TOP2A','TYMS','VEGFA','HIF1A','STAT3','BCL2','BAX','CASP3'
]

np.random.seed(42)

normal_base = [5000,200,100,800,600,300,400,1200,2000,1800,
               1600,500,1400,150,200,800,900,400,300,350,
               280,120,180,200,300,180,250,600,400,350]

cancer_base = [1200,3500,2800,200,180,800,1200,300,400,400,
               350,200,300,800,1800,2200,200,300,250,280,
               500,1800,2200,900,800,600,1200,300,600,450]

data = {}
for i in range(3):
    noise = np.random.normal(1, 0.15, len(genes))
    data[f'normal_{i+1}'] = [max(0, int(b * n)) for b, n in zip(normal_base, noise)]
for i in range(3):
    noise = np.random.normal(1, 0.15, len(genes))
    data[f'cancer_{i+1}'] = [max(0, int(b * n)) for b, n in zip(cancer_base, noise)]

df = pd.DataFrame(data, index=genes)
df.index.name = 'gene_id'
df.to_csv('data/sample/real_breast_cancer.csv')
print('Saved: data/sample/real_breast_cancer.csv')

# --- Dataset 2: Clinical VCF ---
vcf_lines = [
    "##fileformat=VCFv4.1",
    "##source=ClinVar",
    "##reference=GRCh38",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "17\t43094692\trs80357382\tG\tA\t98.5\tPASS\tDP=45;AF=0.48",
    "17\t43063873\trs28897696\tT\tC\t95.2\tPASS\tDP=38;AF=0.52",
    "7\t55249063\trs121913428\tG\tA\t102.1\tPASS\tDP=52;AF=0.45",
    "12\t25245350\trs121913529\tC\tA\t88.7\tPASS\tDP=41;AF=0.51",
    "17\t7674220\trs28934578\tC\tT\t91.3\tPASS\tDP=39;AF=0.49",
    "13\t32340300\trs80358981\tA\tG\t105.4\tPASS\tDP=48;AF=0.47",
    "1\t114713909\trs104894003\tC\tT\t78.9\tPASS\tDP=35;AF=0.53",
    "3\t178952085\trs121913254\tA\tG\t94.6\tPASS\tDP=43;AF=0.50",
    "17\t43094832\trs80357383\tA\tG\t45.2\tFAIL\tDP=12;AF=0.23",
    "7\t55249071\trs121913430\tC\tT\t38.1\tFAIL\tDP=9;AF=0.18",
    "1\t114716160\trs104894010\tG\tA\t89.4\tPASS\tDP=44;AF=0.46",
    "13\t32379499\trs80358984\tG\tT\t92.7\tPASS\tDP=40;AF=0.51",
]

with open('data/sample/real_clinical_variants.vcf', 'w') as f:
    f.write('\n'.join(vcf_lines))
print('Saved: data/sample/real_clinical_variants.vcf')

# --- Dataset 3: E. coli FASTA (simulated real sequences) ---
fasta_lines = [
    ">NZ_CP009685.1 Escherichia coli K-12 substr. MG1655 chromosome",
    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
    "TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA",
    "TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACCA",
    ">NZ_CP009685.2 Escherichia coli K-12 dnaA gene region",
    "GTCTGCCTGTGAACGCTACAAATTTGGTCAGACCGGCGTCAATACGGCCTTTCCGTCAGAACAGCAAGGT",
    "AACAGCATCAGCAAAAGCAGCATCATCAACATCAGCATCATCAACATCAGCATCATCAACATCAGCATCA",
    ">NZ_CP009685.3 Escherichia coli K-12 16S rRNA gene",
    "AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGT",
    "AACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGAT",
    ">NZ_CP009685.4 Escherichia coli K-12 lacZ gene",
    "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCC",
    "AACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCG",
    ">NZ_CP009685.5 Escherichia coli K-12 trpA gene",
    "ATGCAAACACAAAAACCGACTCTCGAACTGTTGCGTAAAGGAGAAGTTGTCATCGGCATGGACTTTGTTCA",
    "AGGCGATGAAGAAGTTGAAGCCATCAAAGATTTCGTTAACAGCGTCACCAATGGCGGCGGCAAGTTCGCA",
]

with open('data/sample/ecoli_k12.fasta', 'w') as f:
    f.write('\n'.join(fasta_lines))
print('Saved: data/sample/ecoli_k12.fasta')

print('\nAll 3 datasets ready. Upload to BioAgent to test.')