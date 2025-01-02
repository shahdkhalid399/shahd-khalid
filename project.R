

# rnaseq_pipeline.sh
#
# This script demonstrates a full RNA-seq workflow in Bash:
# 1. Download SRA data with sra-tools (prefetch, fastq-dump).
# 2. Quality check with FastQC.
# 3. Trimming with Trimmomatic.
# 4. Reference genome index (HISAT2).
# 5. Alignment with HISAT2.
# 6. SAM -> BAM -> sorted BAM using SAMtools.
# 7. Counting reads per gene with featureCounts (Subread).
# 0. Setup Project Directories
echo "Creating directory structure..."
mkdir -p ~/rna_project/{raw_data,qc_reports,trimmed_data,alignment,counts,reference}
# Explanation:
#  - raw_data:      will store the downloaded FASTQ files
#  - qc_reports:    will store FastQC outputs
#  - trimmed_data:  will store trimmed reads
#  - alignment:     will store .sam/.bam alignment files
#  - counts:        will store the featureCounts results
#  - reference:     will store the reference genome + index
echo "Directories created (if they did not exist
# 1. Download Data from SRA (Optional)
# Example SRA run IDs: SRR000111, SRR000112
# Adjust these IDs to real ones from your study.
echo "Starting SRA download..."
cd ~/rna_project/raw_data
prefetch SRR000111  # Download .sra file
prefetch SRR000112
# Convert .sra to FASTQ
fastq-dump --split-files SRR000111
fastq-dump --split-files SRR000112
# Explanation:
#  - 'prefetch' obtains the .sra file from NCBI.
#  - 'fastq-dump --split-files' splits paired-end reads into *_1.fastq and *_2.fastq.
echo "SRA data downloaded and converted to FASTQ."
# 2. Quality Control with FastQC
echo "Running FastQC on raw data..."

fastqc SRR000111_1.fastq SRR000111_2.fastq \
       SRR000112_1.fastq SRR000112_2.fastq \
       -o ../qc_reports

# Explanation:
#  - FastQC checks read quality, adapter contamination, etc.
#  - The '-o' flag specifies the output directory for the HTML and .zip files.

echo "Raw reads QC done. Check the HTML reports in qc_reports folder."
# 3. Trimming with Trimmomatic
echo "Trimming reads with Trimmomatic..."

cd ~/rna_project
trimmomatic PE -threads 4 \
  raw_data/SRR000111_1.fastq raw_data/SRR000111_2.fastq \
  trimmed_data/SRR000111_1_paired.fq.gz trimmed_data/SRR000111_1_unpaired.fq.gz \
  trimmed_data/SRR000111_2_paired.fq.gz trimmed_data/SRR000111_2_unpaired.fq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

trimmomatic PE -threads 4 \
  raw_data/SRR000112_1.fastq raw_data/SRR000112_2.fastq \
  trimmed_data/SRR000112_1_paired.fq.gz trimmed_data/SRR000112_1_unpaired.fq.gz \
  trimmed_data/SRR000112_2_paired.fq.gz trimmed_data/SRR000112_2_unpaired.fq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

echo "Trimming completed. Paired reads are in trimmed_data/."

# 4. Download & Prepare Reference Genome
# Example for human GRCh38 from Ensembl. If you already have a reference genome,
# skip downloading and just build the index. Adjust URLs for your organism.
echo "Downloading and decompressing reference genome..."

cd ~/rna_project/reference
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

echo "Reference genome downloaded and unzipped."
# 5. Build the Index for HISAT2
echo "Building HISAT2 index..."

hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38_index

echo "HISAT2 index built. Files: GRCh38_index.*"
# 6. Alignment with HISAT2
echo "Aligning reads..."

cd ~/rna_project
mkdir -p alignment

hisat2 -p 4 -x reference/GRCh38_index \
  -1 trimmed_data/SRR000111_1_paired.fq.gz \
  -2 trimmed_data/SRR000111_2_paired.fq.gz \
  -S alignment/SRR000111.sam

hisat2 -p 4 -x reference/GRCh38_index \
  -1 trimmed_data/SRR000112_1_paired.fq.gz \
  -2 trimmed_data/SRR000112_2_paired.fq.gz \
  -S alignment/SRR000112.sam

echo "Alignment done. SAM files are in alignment/."
# 7. Convert SAM -> Sorted BAM -> Index with SAMtools
echo "Converting, sorting, and indexing alignments..."

cd alignment

samtools view -bS SRR000111.sam | samtools sort -o SRR000111.sorted.bam
samtools index SRR000111.sorted.bam

samtools view -bS SRR000112.sam | samtools sort -o SRR000112.sorted.bam
samtools index SRR000112.sorted.bam

echo "Sorted BAM files ready. Index files created."

# 8. Gene-Level Counting (featureCounts)
echo "Downloading annotation (GTF) for featureCounts..."

cd ~/rna_project/reference
wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz

cd ~/rna_project
mkdir -p counts

echo "Running featureCounts..."

featureCounts -T 4 \
  -a reference/Homo_sapiens.GRCh38.109.gtf \
  -o counts/all_samples_counts.txt \
  alignment/SRR000111.sorted.bam \
  alignment/SRR000112.sorted.bam

# Explanation:
#  - '-T 4': use 4 threads
#  - '-a': path to GTF annotation
#  - '-o': output file
#  - Provide sorted BAM files from the alignment step

echo "featureCounts completed. Output saved to counts/all_samples_counts.txt."


echo "RNA-seq pipeline completed successfully!"
# End of rnaseq_pipeline.sh
# deseq2_analysis.R
#
# This script demonstrates how to:
# 1. Load featureCounts output into R
# 2. Run a basic DESeq2 analysis
# 3. (Optional) Visualize a KEGG pathway via pathview

# 0. Install DESeq2 and pathview if needed:
# (Uncomment if you haven't already installed them in your conda environment)
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("DESeq2")
# BiocManager::install("pathview")

library(DESeq2)
library(pathview)
# 1. Load the featureCounts file
counts_file <- "counts/all_samples_counts.txt"

# featureCounts output has comment lines at the top that start with '#'.
# We'll skip them by using comment.char="#"
raw_counts <- read.table(counts_file, header=TRUE, row.names=1, comment.char="#")

# The actual count columns typically start at column 7 in featureCounts output
# but let's confirm how the columns are arranged. The last columns are the samples.
countData <- raw_counts[,6:ncol(raw_counts)]

# Check the first few rows/columns
print("Preview of countData:")
print(head(countData))
# 2. Create Sample Metadata (colData)
# Let's assume we have 2 samples: SRR000111 and SRR000112.
# We must match the column names from 'countData' exactly.
# Example column headers might be "alignment/SRR000111.sorted.bam"
# and "alignment/SRR000112.sorted.bam".
#
# We define conditions: e.g., SRR000111 is "control", SRR000112 is "treated".
sampleNames <- colnames(countData)
condition   <- c("control", "treated")  # Adjust as needed

colData <- data.frame(
  row.names = sampleNames,
  condition = factor(condition, levels=c("control","treated"))
)

print("colData:")
print(colData)
# 3. Create DESeqDataSet and Run DE Analysis
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design=~condition)

# Run the DESeq2 pipeline
dds <- DESeq(dds)

# Extract results (contrast: treated vs. control)
res <- results(dds, contrast=c("condition","treated","control"))
res <- res[order(res$padj), ]  # sort by adjusted p-value ascending
head(res)
# 4. Export the DE Results to a CSV file
write.csv(as.data.frame(res), file="deseq2_results.csv")
cat("DESeq2 results saved to 'deseq2_results.csv'.\
# 5. Optional: Basic Visualization
# Plot an MA-plot
png("MA_plot.png", width=600, height=600)
plotMA(res, main="DESeq2 MA-Plot", ylim=c(-5,5))
dev.off()
cat("MA-plot saved to 'MA_plot.png'
# 6. Optional: Pathway Visualization with pathview

# Example: let's say we want to visualize "Alzheimer’s disease" KEGG pathway (hsa05010).
# We'll create a named numeric vector of log2 fold changes.
gene_fc <- res$log2FoldChange
# We need the gene symbols or Entrez IDs as names. If rownames(res)
# are Ensembl IDs, you might have to map them to Entrez or symbol for pathview.
# For demonstration, we assume rownames(res) are gene symbols:
names(gene_fc) <- rownames(res)

# pathview call
# "hsa05010" is the KEGG ID for "Alzheimer’s disease" in humans (hsa = homo sapiens).
pv_out <- pathview(
  gene.data   = gene_fc,
  pathway.id  = "hsa05010",
  species     = "hsa",
  out.suffix  = "alz",
  kegg.native = TRUE
)
# This generates PNG and/or PDF in your working directory.

cat("Pathview: see generated files for 'hsa05010_alz.png' or 'hsa05010_alz.pdf'.\n")

# End of deseq2_analysis.R
