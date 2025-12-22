# 16S rRNA V3–V4 (341f/805r) ASV Analysis  
## Execution Manual for the NIG Supercomputer

**Target:**  
16S rRNA (V3–V4 region) ASV analysis using **SILVA v138.2** and **DADA2**

---

## Table of Contents

1. [Overview of the Analysis](#1-overview-of-the-analysis)  
2. [Prerequisites (Required)](#2-prerequisites-required)  
3. [Project Structure](#3-project-structure)  
4. [Script Preparation and Configuration](#4-script-preparation-and-configuration)  
5. [Job Submission and Execution](#5-job-submission-and-execution)  
6. [Result Inspection and Post-processing](#6-result-inspection-and-post-processing)  
7. [Parameter Tuning Guide](#7-parameter-tuning-guide)  
8. [List of Output Files](#8-list-of-output-files)  
9. [Downloading the SILVA Database (For Instructors)](#9-downloading-the-silva-database-for-instructors)

---

## 1. Overview of the Analysis

### 1-1. DADA2 Analysis Pipeline

```

FASTQ placement → Main analysis (js1) → Result inspection → Post-processing (js2) → Completion

````

**Processes executed in the main analysis:**

1. Loading FASTQ files (automatic detection of R1/R2)
2. Generation of quality plots (PDF output)
3. Filtering (removal of low-quality regions and primers)
4. Error model learning (`learnErrors`)
5. ASV inference (`dada`)
6. Paired-end merging (`mergePairs`)
7. Chimera removal (`removeBimeraDenovo`)
8. ASV table generation (saved as RDS)
9. Taxonomic assignment (`assignTaxonomy`)
10. Read tracking (recording read counts at all steps)

**Processes executed in post-processing:**

- Output of ASV sequences in FASTA format
- Generation of count tables and taxonomy tables
- Calculation of relative abundance (%)
- Aggregation by taxonomic rank (Kingdom to Genus)

---

### 1-2. Tools Used

- **Analysis method:** DADA2  
- **Database:** SILVA v138.2  
- **Execution environment:** NIG supercomputer  

---

## 2. Prerequisites (Required)

### 2-1. Placement of FASTQ Files (Most Important)

**Procedure:**

1. On your local computer, open the **“16S analysis result folder”** returned from the sequencing facility
2. Open the **`raw_fastq/`** folder inside it
3. Copy **all `*.fastq.gz` files** in that folder to the following directory on the NIG supercomputer:

```bash
~/16S_SILVA/raw_fastq/
````

**Important notes:**

* Place **all FASTQ files in the same `raw_fastq/` directory**
* Do not separate files by sample or experiment into different folders
* When additional data arrive, **add them to the same directory**
* File names can be kept as provided by the sequencing facility
  (e.g., `*_R1_001.fastq.gz`, `*_R2_001.fastq.gz`)

---

### 2-2. Logging in to the NIG Supercomputer

```bash
# Log in to the gateway node
↓

# Log in to a compute node
e.g., ssh a001
```

---

### 2-3. Creating Required Directories

```bash
# Create project root directories
mkdir -p ~/16S_SILVA/raw_fastq
mkdir -p ~/16S_SILVA/scripts

# Create log directory
mkdir -p ~/log
```

**Note:**
If `~/log/` is not created, job submission will fail with an error.

---

## 3. Project Structure

### 3-1. Directory Structure

**Example structure:**

```

~/16S_SILVA/                    ← Project root
├── raw_fastq/                  ← FASTQ file directory
│   ├── Sample1_R1_001.fastq.gz
│   ├── Sample1_R2_001.fastq.gz
│   ├── Sample2_R1_001.fastq.gz
│   └── Sample2_R2_001.fastq.gz
├── scripts/                    ← Scripts
│   ├── s1_dada2_pipeline.r
│   ├── s2_asv_postprocess.r
│   ├── js1_dada2_pipeline.sh
│   └── js2_asv_postprocess.sh
└── output/                     ← Analysis results (automatically created)
├── filtered/               ← Filtered FASTQ files
├── ASV_nochim.rds
├── taxonomy.rds
└── ... (other output files)

```

**Location of the SILVA database (example; update the path once the instructor’s SILVA database location is finalized):**
```

/home/ryohonda/db/dada2/SILVA_138.2/silva_nr99_v138.2_toGenus_trainset.fa.gz

````
(inside Prof. Honda’s directory)

---

### 3-2. What Is the Project Root?

The **project root** refers to the top-level directory that contains all materials required for the analysis.

In the scripts, modify the following line to match your own project root:

```r
PROJECT_DIR <- "/home/username/16S_SILVA"
````

---

## 4. Script Preparation and Configuration

### 4-1. Placement of Script Files

Place the following four files in `~/16S_SILVA/scripts/`:

1. **s1_dada2_pipeline.r** – R script for the main analysis
2. **s2_asv_postprocess.r** – R script for post-processing
3. **js1_dada2_pipeline.sh** – Shell script for submitting the main analysis job
4. **js2_asv_postprocess.sh** – Shell script for submitting the post-processing job

---

### 4-2. Editing the R Scripts

#### Sections to edit in `s1_dada2_pipeline.r`

**1. Specify the project directory:**

```r
# ★ Replace with your own username
PROJECT_DIR <- "/home/username/16S_SILVA"
FASTQ_DIR   <- file.path(PROJECT_DIR, "raw_fastq")
OUT_DIR     <- file.path(PROJECT_DIR, "output")
```

**2. Specify the SILVA database path:**

```r
# ★ Replace with Prof. Honda’s path (or the location where you downloaded it)
silva_file <- "/home/ryohonda/db/dada2/SILVA_138.2/silva_nr99_v138.2_toGenus_trainset.fa.gz"   # ★ Change to the instructor’s DB path
```

**3. Filtering parameters (adjust as needed):**

```r
out <- filterAndTrim(
    fwd = fnFs, filt = filtFs,
    rev = fnRs, filt.rev = filtRs,
    trimLeft = c(17, 21),      # Primer length (V3–V4: 341F = 17 bp, 805R = 21 bp)
    truncLen = c(280, 220),    # Adjust based on quality plots
    maxN = 0, maxEE = c(2, 2),
    truncQ = 2, compress = TRUE, multithread = TRUE
)
```

---

#### Sections to edit in `s2_asv_postprocess.r`

**Specify the output directory:**

```r
# ★ Replace with your own username
OUT_DIR <- "/home/username/16S_SILVA/output"
```

---

### 4-3. Creating Job Submission Scripts

#### `js1_dada2_pipeline.sh` (Main analysis)

This step requires a moderate amount of computational resources.
Using the **medium** or **rome** partition is recommended.

* Request **8 CPU cores** and **32 GB memory**
* If the number of samples exceeds 100, request **64 GB memory**
  (note that queue waiting time will be longer)

```bash
#!/bin/bash
#SBATCH --output=/home/username/log/%x_%j.out.log    # ★ Replace username
#SBATCH --error=/home/username/log/%x_%j.err.log     # ★ Replace username
#SBATCH -p medium       # Change if necessary
#SBATCH -N 1-1
#SBATCH -n 8
#SBATCH --mem=32G       # Sufficient for ~45 samples (adjust if needed)

# Execute the script (★ update the path)
Rscript /home/username/16S_SILVA/scripts/s1_dada2_pipeline.r
```

---

#### `js2_asv_postprocess.sh` (Post-processing)

This step is lightweight, so the **short** partition is sufficient.
It can also be run on a local PC.

```bash
#!/bin/bash

#SBATCH --output=/home/username/log/%x_%j.out.log    # ★ Replace username
#SBATCH --error=/home/username/log/%x_%j.err.log     # ★ Replace username
#SBATCH -p short

# Execute the script (★ update the path)
Rscript /home/username/16S_SILVA/scripts/s2_asv_postprocess.r
```

---

### 4-4. Granting Execute Permission

```bash
cd ~/16S_SILVA/scripts
chmod u+x js1_dada2_pipeline.sh
chmod u+x js2_asv_postprocess.sh
```

---

## 5. Job Submission and Execution

### 5-1. Submitting the Main Analysis Job

```bash
cd ~/16S_SILVA/scripts
sbatch js1_dada2_pipeline.sh
````

**Message displayed upon successful submission:**

```
Submitted batch job 123456
```

The number `123456` is the **job ID**.

---

### 5-2. Checking Job Status

```bash
# Check the status of your jobs
squeue -u username
```

**Example output:**

```
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
123456    medium dada2_me username  R       1:23      1 node01
```

**Meaning of job status codes:**

* `R` (Running): Job is currently running
* `PD` (Pending): Job is waiting for available resources
* No output: Job has completed or terminated

---

### 5-3. Monitoring Logs in Real Time

```bash
# Standard output log (analysis progress)
tail -f ~/log/dada2_medium_123456.out.log

# Error log (error messages)
tail -f ~/log/dada2_medium_123456.err.log
```

**To stop monitoring:** press `Ctrl + C`

---

### 5-4. Cancelling a Job (If Necessary)

```bash
scancel 123456
```

---

## 6. Result Inspection and Execution of Post-processing

### 6-1. Mandatory Checks After Completion of the Main Analysis

**Before running the post-processing step (`js2`), be sure to check the following items:**

---

#### Check 1: Quality Plot Inspection

Quality plots generated before filtering:

* `R1_quality_before.pdf`
* `R2_quality_before.pdf`

Quality plots generated after filtering (to evaluate filtering effects):

* `R1_quality_after.pdf`
* `R2_quality_after.pdf`

**Download the PDFs and inspect them carefully.**

**Key points to check:**

* Identify regions where the quality score (green band) for R1 and R2 remains above Q30
* Assess whether the positions specified by `truncLen` are appropriate

---

#### Check 2: Read Tracking Inspection

```bash
cat read_tracking.tsv
```

**Key points to check:**

* Whether the number of reads drops excessively at any step
* In particular, confirm that the number of reads in `merged` and `nonchim` is not drastically reduced compared with `input`
* **`final_perc_reads_retained` (final retention rate) is calculated**

  * This value indicates the percentage of reads retained relative to the initial read count
  * Formula: `(nonchim / input) × 100`

---

**Good example:**

```
           input filtered denoisedF denoisedR merged nonchim final_perc_reads_retained
Sample1   150000   140000    135000    130000  90000   88000  58.7
Sample2   120000   110000    105000    100000  70000   68000  56.7
```

→ Retention rate is above 50%, with no extreme drop at any step

---

**Bad example (requires parameter re-tuning):**

```
           input filtered denoisedF denoisedR merged nonchim final_perc_reads_retained
Sample1   150000   140000    135000    130000   5000    4800   3.2
Sample2   120000   110000    105000    100000   3000    2900   2.4
```

---

```
Actual workflow in practice:

First analysis run:
  (1) Inspect "before" quality plots → decide truncLen
  (2) Run the script (including filterAndTrim)
  (3) Inspect "after" quality plots → evaluate filtering effects
  (4) Check read_tracking.tsv → confirm read counts at each step
     ↓ If results are unsatisfactory (e.g., too few merged reads)

Second analysis run:
  (1) Modify truncLen and re-run the s1 script
  (2) Re-check "after" quality plots
     ↓ If acceptable
  Proceed to the next step
```

#### Check 3: Inspect Log Files

```bash
# Confirm the "analysis completed" message
grep "Done/解析完了" ~/log/dada2_medium_123456.out.log

# Check for error messages
cat ~/log/dada2_medium_123456.err.log
````

---

### 6-2. Troubleshooting When Problems Occur

**Problem: Low merging rate (few merged reads)**

→ Adjust `truncLen` and **re-run the main analysis**

For details, see “[7. Parameter Tuning Guide](#7-parameter-tuning-guide)”.

---

### 6-3. Submitting the Post-processing Job

**If all checks are completed and no problems are found:**

```bash
cd ~/16S_SILVA/scripts
sbatch js2_asv_postprocess.sh
```

**Check job status:**

```bash
squeue -u username
tail -f ~/log/dada2_post_*.out.log
```

**Confirm completion:**

```bash
# Confirm the "ASV post-processing completed" message
grep "Done/後処理完了" ~/log/dada2_post_*.out.log
```

---

## 7. Parameter Tuning Guide

### 7-1. How to Determine `truncLen` (Read Truncation Length)

```r
truncLen = c(280, 220)  # R1 = 280 bp, R2 = 220 bp
```

**Decision procedure:**

1. **Inspect quality plots**

   * Truncate reads just before the quality (green band) drops below Q30

2. **Ensure sufficient overlap**

   * The overlapping region between R1 and R2 should be **~50 bp, which is sufficient**
   * Too much or too little overlap is undesirable (as low as ~20 bp can still work)
   * Formula:
     `overlap = truncLen_R1 + truncLen_R2 - amplicon_length`
   * Example:
     `280 + 220 - 460 = 40 bp` → OK
   * A rough guideline is `R1 + R2 ≈ 510 bp` (460 + 50)
   * For example, R1 = 280 and R2 = 230
   * R1 tends to be longer than R2

---

**Examples of parameter choices:**

| Pattern     | R1  | R2  | Overlap | Evaluation                      |
| ----------- | --- | --- | ------- | ------------------------------- |
| Recommended | 280 | 220 | 40 bp   | Good                            |
| Too short   | 240 | 180 | -40 bp  | Merge failure                   |
| Too long    | 300 | 260 | 100 bp  | May include low-quality regions |

---

### 7-2. How to Determine `trimLeft` (Primer Removal)

```r
trimLeft = c(17, 21)  # 341F = 17 bp, 805R = 21 bp
```

**Decision procedure:**

Count the number of bases **after the hyphen (`-`)** in the primer sequences used.

**V3–V4 primers from the sequencing facility (from the work report):**

```
1st-341f: ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNN-CCTACGGGNGGCWGCAG
                                                  └─────────┬─────────┘
                                                         17 bp

1st-805r: GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNN-GACTACHVGGGTATCTAATCC
                                                  └──────────┬──────────┘
                                                          21 bp
```

**Explanation:**

* **Before** the hyphen (`-`): adapter sequence (already removed prior to DADA2)
* **After** the hyphen (`-`): actual primer sequence (**must be removed by DADA2**)

**Calculation:**

* 341F: `CCTACGGGNGGCWGCAG` → 17 bases
* 805R: `GACTACHVGGGTATCTAATCC` → 21 bases

→ `trimLeft = c(17, 21)`

**Important:**
If different primers are used, be sure to modify these values accordingly.

---

### 7-3. Other Parameters

```r
maxN = 0           # Remove reads containing ambiguous bases (N)
maxEE = c(2, 2)    # Maximum expected errors (R1, R2)
truncQ = 2         # Truncate at positions with quality score ≤ 2
```

**In most cases, these parameters do not need to be changed.**

---
## 8. List of Output Files

### 8-1. Outputs from the Main Analysis (after running js1)

| File name | Description |
|----------|-------------|
| `ASV_seqtab.rds` | ASV table before chimera removal |
| `ASV_nochim.rds` | ASV table after chimera removal |
| `taxonomy.rds` | ASV × taxonomy table |
| `filter_stats.tsv` | Filtering statistics |
| `read_tracking.tsv` | Read count tracking across all steps |
| `R1_quality_before.pdf` | R1 quality plot before filtering (★ must be checked) |
| `R2_quality_before.pdf` | R2 quality plot before filtering (★ must be checked) |
| `R1_quality_after.pdf` | R1 quality plot after filtering (★ must be checked) |
| `R2_quality_after.pdf` | R2 quality plot after filtering (★ must be checked) |
| `filtered/` | Filtered FASTQ files (intermediate files) |

---

### 8-2. Outputs from Post-processing (after running js2)

| File name | Description |
|----------|-------------|
| `ASV_sequences.fasta` | DNA sequences of all ASVs (FASTA format) |
| `ASV_counts.tsv` | ASV × sample count table |
| `ASV_taxonomy.tsv` | ASV × taxonomic annotation table |
| `ASV_results.tsv` | Integrated table of counts and taxonomy |
| `ASV_results_relabund.tsv` | Relative abundance (%) table |
| `ASV_Kingdom.tsv` | Aggregated results at Kingdom level |
| `ASV_Phylum.tsv` | Aggregated results at Phylum level |
| `ASV_Class.tsv` | Aggregated results at Class level |
| `ASV_Order.tsv` | Aggregated results at Order level |
| `ASV_Family.tsv` | Aggregated results at Family level |
| `ASV_Genus.tsv` | Aggregated results at Genus level |
| `dada2_post_workspace.RData` | R workspace file (all data saved) |

---

## 9. Downloading the SILVA Database (For Instructors)

**Download source:**  
https://www.arb-silva.de/current-release/DADA2/1.36.0/SSU

**Required file:**
```

silva_nr99_v138.2_toGenus_trainset.fa.gz

```

---

## Revision History

- 2025/12/15: Initial version created (Miho Kobayashi)

---

## Personal Notes (Kobayashi)

# Illustrated Explanation of Paired-end Sequencing and Overlap

## 1. Basics of Paired-end Sequencing

### 1-1. What the Illumina Sequencer Does

```

[Actual DNA fragment (V3–V4 region)]
5'==================================================3'  ~460 bp
↑                                              ↑
341F                                          805R
Primer                                        Primer

[How Illumina reads]
R1 (Forward) →→→→→→→→→
           5'==================================================3'
           3'==================================================5'
                                                     ←←←←←←←←← R2 (Reverse)

* R1: Read forward from the 5' end (typically 280–300 bp)
* R2: Read backward from the 3' end (typically 220–280 bp)

```

---

### 1-2. Why Reads Are Generated from Both Ends

**Reason 1:** Illumina sequencers cannot read long sequences in a single run  
(approximately 300 bp is the practical limit)

**Reason 2:** The V3–V4 region is ~460 bp long, so it must be read from both ends and merged in the middle

---

## 2. What Is an Overlap?

### 2-1. Visual Explanation

```

[Case with overlap (normal)]

  R1: 280 bp read
5'━━━━━━━━━━━━━━━━━━━━━━━━━━━→
                      ■■■■■■■  ← overlapping region
                     ←━━━━━━━━━━━━━━━━━━━━━━3'
                      R2: 220 bp read

Overlap = 280 + 220 − 460 = 40 bp

[Case without overlap (merge failure)]

  R1: 240 bp read
5'━━━━━━━━━━━━━━━━→
                   ?????? ← unread region (gap)
                         ←━━━━━━━━━━━━━3'
                          R2: 180 bp read

Overlap = 240 + 180 − 460 = −40 bp (negative!)
→ Cannot be merged!

```

---

### 2-2. Formula for Calculating Overlap

```

Overlap (bp) = length of R1 + length of R2 − actual length of the DNA fragment

Example 1: truncLen = c(280, 220), V3–V4 = 460 bp
→ 280 + 220 − 460 = 40 bp  OK

Example 2: truncLen = c(240, 180), V3–V4 = 460 bp
→ 240 + 180 − 460 = −40 bp  NG (gap occurs)

Example 3: truncLen = c(300, 260), V3–V4 = 460 bp
→ 300 + 260 − 460 = 100 bp  Too much overlap (may include low-quality regions)

```

---
## 3. What Is Merging (`merge`)?

### 3-1. Meaning of Merging

**Merging = connecting R1 and R2 into a single continuous sequence**

```

[Before merging]
R1: ATCGATCGATCG... (Forward read)
R2: CGTAGCTAGCTA... (Reverse read)

[After merging]
One complete sequence: ATCGATCGATCG...CGTAGCTAGCTA

```

---

### 3-2. DADA2 Merging Process

```

1. Identify the overlapping region between R1 and R2

2. Check whether the sequences in the overlapping region match
   ┌─────┐
   │ATCGC│ ← end of R1
   │ATCGC│ ← start of R2 (after reverse-complement conversion)
   └─────┘
   Match! → Merge successful

3. If they do not match → Merge fails (the read is discarded)

```

---

End of document.
```




