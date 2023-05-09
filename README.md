# customscripts

A dirty repo with my custom scripts for any project.

## BASH

### GR_rnaseq_pipe_launcher.sh

A bash script to perform a RNAseq analysis on flamingo, using the pipelines from Thibault DAYRIS :

* "Manual" QC on raw reads (fastqc, fastq_screen, multiqc). This step is out of the pipeline despite it could perform fastqc/multiqc, due to a current impossibility to get the fastqc results from the pipeline (the pipeline attempts to write in an unauthorized directory). Should be fixed someday !
* "Manual" trimming (fastp) when needed. This step is out of the pipeline as intended during its coopted creation with the team.
* "Manual" QC on trimmed reads.
* Pseudo-mapping and quantification (salmon) using the 'rna-count-salmon' pipeline.
* Immune infiltration estimation (6 tools) using the 'immune-decov' pipelone.
* Differential gene expression analysis using (deseq2) the 'rna-dge-salmon-deseq2' pipeline.

At each step, results are zipped and automatically sent to my NextCloud shared directory.

### get_data_irods_byrun.sh

A bash script to retrieve data from the warehouse using iRODS. Data will be written in the default data_input directory of the project, with a sub-directory corresponding to the sequencing run, then sub-directories corresponding to dataset_id (to avoid smashing older data with the same name, from an older run by a newer one)

### get_data_irods_simple.sh

Same as get_data_irods_byrun.sh, but without taking care of the sequencing run.


## R

### AnnovarTSV_reformat.R

Aggregation then filtering of Annovar-annotated TSV files  tables (from the Agilent SureSelect XTHS / HaloPlex HS pipelines, corresponding to tabular outputs from Varscan2/FreeBayes+Annovar), and additional filtering on variant frequency, ExAC_ALL frequency, and refGene functions. 

### bedGC.R

Computes GC% from a bed file

### ChAMP_wrapper_v2.20.1.R

A R wrapper to use ChAMP on illumina methylation microarrays (450K or EPIC designs).
Tested on R v4.0.4 with ChAMP v2.20.1

### ChiFisher.R

Performs a series of Fisher / X2  or Wilcoxon / Kruskal-Wallis tests on selected query / target columns from a tab-separated annotation file

### DEseqObj2normalizedmatrix.R

Normalization of a RAW COUNTS DEseq2 object (DESeqDataSet) with vst, returns the normalized counts matrix.
Extra parameters (...) are passe to DESeq2::vst()
Requires DESeq2 and clusterProfiler packages, and functions from customscripts/R/diffexp2gsea_design.R

### diffexp.R

Set of functions to perform differential gene expression using DESeq2.
Additional function also allows to assess covariate, in a way to evaluate which to regress towards inceasing the expected biological signal.

### diffexp_design.R (NEWER)

Same as diffexp.R, but using a design table (more convenient, allows a fine control of compared entities).

### diffexp2gsea.R

Automate GSEA/ORA functional annotation and analysis, based on clusterProfiler/DOSE.
Compatible with MSigDb (thanks to the msigdbr package), KEGG, GO, DO, WikiPathways, Reactome, KEGG/MKEGG, Mesh

### EaCoN_TCN_GIS_autoscorer.R

Generates a table containing GIS (genomic instability scores) from an EaCoN TCN output using the ASCAT segmenter.
This is an executable script that requires arguments.

### HTG_analysis_functions.R

Set of functions to perform the analysis of HTG EdgeSeq target RNAseq data. Requires diffexp or diffextp_design.

### immune-deconv_difftest.R

Differential analysis on immune cellularity prediction results (from the 'immune-decov' pipeline) on (clinical) sample annotations, using a Wilcoxon sum-rank test (T-test, optionally).

### maelstrom.R

A series of functions :

* Device
    * rasterpdf.open() : Opens a "connection" to output multiple plots to a multipage TIFF with rasters
    * rasterpdf.close() : Closes the "connexion" opened with rasterpdf.open()
    * multipng2pdf() : Concatenates PNG files to a multipage (raster) pdf
* Parsing
    * write.table.fast() : Fast file writer using iotools::write.csv.raw but corrected for header handling (~ 5x faster)
    * read.table.fast() : Fast file reader using data.table::fread
    * db.load() : Load data from a SQLite ".db" file, using DBI and RSQLite. Returns a list of tables.
    * hdf5.load() : Load data from a HDF5 file, using rhdf5. Returns a list of tables.
    * read.horiz.csv() : Read a csv/tsv file with data ordered HORIZONTALLY (to a df)
* Conversion
    * factors2char.df() : Converts any factor column to a character column in a dataframe
    * factors2num.df() : Converts any factor column to a numeric column in a dataframe
    * chrConv() : Converts chrom <-> chr (ie, "chr1" <-> 1) with support to alphanumerical output. If alpha = TRUE, "chrX" -> "X", else "chrX" -> 23 (for homo sapiens)
* Matrix
    * rotate.matrix.clockwise() : Rotates a matrix (90 degrees, clockwise)
    * matrix.ranks.scaler() : Scales the columns of a matrix using a given vector of values. This is used to scale additional samples to another quantiles-normalized matrix, using its ranks.
    * matrix.rows.merger() : Merges values from multiple lines of a numerical matrix according to their same rowname. Supported methods are : median, mean, min, max
    * matrix.rows.aggregator() : Merges values from multiple lines of a numerical matrix according to their same rowname (aggregate version). Supported methods are any R function that can coerce a vector to a single value (such as : median, mean, min, max, ...)
    * matrix.na.replace.byrow() : Replace missing values in a numerical matrix using the median / mean / min or max of corresponding line. Please only use this if knn imputation failed (see package "impute")
* Vector
    * interleave.numeric() : Interleaves two NUMERIC vectors into one
* Plots
    * pca2d() : Biplots for PCA results, with centroids
    * pca3d() : 3D-plot of PCA with optional classes and projections
* System
    * get.os() : A more robust way to get machine OS type
* Genomics
    * bed2gc() : Computes GC% fro ma bed file

### NMF_run.R

A wrapper to ease the use of the NMF clusteringpackage for R. This allows NMF clustering using various methods, along with capturing the features/variables contributing to the different clusters, and assess the clusters to a null distribution by shuffling data.

### Recount3.R

Very preliminary script to mess with ReCount3 data, thanks to the recount3 R package.

### skmeans_run.R

A wrapper to ease the use of the skmeans package for R. This allows spherical-kmeans clustering using various methods.

### SNV_Panel_F2_aggregator.R

Performs the aggregation of "F2" tables (from the Agilent SureSelect XTHS / HaloPlex HS pipelines, corresponding to tabular outputs from Varscan2+Annovar), and additional filtering on variant frequency, ExAC_ALL frequency, and refGene functions.

### SoupX_usage.R

Performs removal of ambient RNA expression for a single cell counts matrix using SoupX. Usage of the raw (ie, without discarding the empty barcodes) is recommended. Without, lower efficiency is expected. TESTED WITH SoupX v1.6.2

### WTF.R

Performs a Wilcoxon rank test (W), a students' T-test (T), and/or a Fisher's exact test (F) on a continuous variable and a class.
