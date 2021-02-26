# customscripts

A dirty repo with my custom scripts for any project.

* GR_rnaseq_pipe_launcher.sh

A bash script to perform a RNAseq analysis on flamingo, using the pipelines from Thibault DAYRIS :

    * "Manual" QC on raw reads (fastqc, fastq_screen, multiqc). This step is out of the pipeline despite it could perform fastqc/multiqc, due to a current impossibility to get the fastqc results from the pipeline (the pipeline attempts to write in an unauthorized directory). Should be fixed someday !
    * "Manual" trimming (fastp) when needed. This step is out of the pipeline as intended during its coopted creation with the team.
    * "Manual" QC on trimmed reads.
    * Pseudo-mapping and quantification (salmon) using the 'rna-count-salmon' pipeline.
    * Immune infiltration estimation (6 tools) using the 'immune-decov' pipelone.
    * Differential gene expression analysis using (deseq2) the 'rna-dge-salmon-deseq2' pipeline.

At each step, results are zipped and automatically sent to my NextCloud shared directory.

* GR_wes_pipe_launcher.sh
