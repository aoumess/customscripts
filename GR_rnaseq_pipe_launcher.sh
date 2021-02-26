## NEXTCLOUD
NCUSER='B_JOB'
NCCRED=${NCUSER}':Z03F4nny!0104'
NCROOTDIR='https://nextcloud.gustaveroussy.fr/remote.php/dav/files/B_JOB/'

## PARAMETERS
R1CODE='R1'
R2CODE='R2'
MINBASEQ=20  # MINIMUM BASE QUALITY FOR FILTERING
ZCOMP=6
ADAPTERSFASTA='/mnt/beegfs/scratch/bioinfo_core/B20029_SOGA_02/data_input/RESOURCES/ADAPTERS/ADAPTERS_CONCAT.fa'
FASTQSCREENCONF='/mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/fastq_screen.conf'
MYTMP='/mnt/beegfs/userdata/b_job/mytmp/'
MYGITSREP='/mnt/beegfs/userdata/b_job/gits/'
NCDIR='https://nextcloud.gustaveroussy.fr/remote.php/dav/files/B_JOB/B20067_FRDA_05/'



## PATHS

### P29_AUMA
PFBPROJECTNAME='B20067_FRDA_05'
PFGPROJECTNAME='P29_AUMA'
RUNNAME='CONCATENATED'
SNAMEBLOCKS="1,2"  # PART(S) OF THE FASTQ FILENAMES THAT WILL FORM THE SAMPLENAME

### P30_AUMA
PFBPROJECTNAME='B20067_FRDA_05'
PFGPROJECTNAME='P30_AUMA'
# RUNNAME='191114_A00461_0089_AHM7J3DMXX'
# RUNNAME='191202_A00461_0091_AHHYJHDRXX'
# RUNNAME='191216_A00461_0093_BHJGCTDRXX'
# RUNNAME='200109_A00461_0099_BHJM2GDRXX'
# RUNNAME='200130_A00461_0102_BHJNG3DRXX'
RUNNAME='CONCATENATED'
SNAMEBLOCKS="1,2"  # PART(S) OF THE FASTQ FILENAMES THAT WILL FORM THE SAMPLENAME

### P31_DELE
PFBPROJECTNAME='B21002_DELE_01'
PFGPROJECTNAME='P31_DELE'
RUNNAME='210211_A00461_0161_AH2KTKDRXY'
SNAMEBLOCKS="1,2"

## COMMON
INDIR='/mnt/beegfs/scratch/bioinfo_core/'${PFBPROJECTNAME}'/data_input/'
OUTDIR='/mnt/beegfs/scratch/bioinfo_core/'${PFBPROJECTNAME}'/data_output/'
RAWDATADIR=${INDIR}'/'${PFGPROJECTNAME}'/'${RUNNAME}
WORKDIR=${OUTDIR}'/'${PFGPROJECTNAME}'/'${RUNNAME}
RAWQCDIR=${WORKDIR}'/RAW_QC'
TRIMDIR=${WORKDIR}'/TRIM'
SALMONDIR=${WORKDIR}'/SALMON'
IMMUNEDIR=${WORKDIR}'/IMMUNE'


## Creating project dir
NCOUTDIR=${NCROOTDIR}'/'${PFBPROJECTNAME}
curl --user ${NCCRED} -X MKCOL $NCOUTDIR
NCOUTDIR=${NCOUTDIR}'/'${PFGPROJECTNAME}
curl --user ${NCCRED} -X MKCOL $NCOUTDIR
NCOUTDIR=${NCOUTDIR}'/'${RUNNAME}'/'
curl --user ${NCCRED} -X MKCOL $NCOUTDIR

NCOUTDIR=${NCROOTDIR}'/'${PFBPROJECTNAME}'/'${PFGPROJECTNAME}'/'${RUNNAME}'/'

## GENERATING SAMPLENAMES
SNARRAY=`ls ${RAWDATADIR}/*.f*q* | xargs -n 1 basename | cut -d "_" -f ${SNAMEBLOCKS}`
SNAMES=`printf "%s\n" "${SNARRAY[@]}" | sort -u`

## PER-RUN RAW QC
NTHREADS=2
MEM=5
conda activate readsqc
mkdir -p ${RAWQCDIR}'/fastqc' ${RAWQCDIR}'/fastqscreen' ${RAWQCDIR}'/multiqc' ${RAWQCDIR}'/logs'
cd ${RAWQCDIR}
for SAMPLENAME in ${SNAMES[*]}; do {
	echo ${SAMPLENAME};
	TOOLEO='rawqc'
	sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J ${TOOLEO}'.'${SAMPLENAME} \
	-e ${RAWQCDIR}'/logs/'${SAMPLENAME}'_'${TOOLEO}'.e' -o ${RAWQCDIR}'/logs/'${SAMPLENAME}'_'${TOOLEO}'.o' --wrap \
	"fastqc -o ${RAWQCDIR}/fastqc -t ${NTHREADS} ${RAWDATADIR}/*${SAMPLENAME}*_${R1CODE}*.f*gz \
	${RAWDATADIR}/*${SAMPLENAME}*_${R2CODE}*.f*gz && \
	fastq_screen --threads ${NTHREADS} --force --outdir ${RAWQCDIR}/fastqscreen --subset 100000 \
	--conf ${FASTQSCREENCONF} ${RAWDATADIR}/*${SAMPLENAME}*_${R1CODE}*.f*gz && \
	multiqc -n RAW_${SAMPLENAME} -i ${PROJECTNAME}' '$SAMPLENAME' RAW FASTQ' -z -f -o ${RAWQCDIR}'/multiqc' ${RAWQCDIR}/fastqscreen/*${SAMPLENAME}* ${RAWQCDIR}/fastqc/*${SAMPLENAME}*.zip"
}
done

NTHREADS=3
MEM=5
sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'mqc_RAWQC' -e ${RAWQCDIR}'/logs/'${PFBPROJECTNAME}'_mqc_raw.e' -o ${RAWQCDIR}'/logs/'${PFBPROJECTNAME}'_mqc_raw.o' --wrap \
"multiqc -n ${PFBPROJECTNAME}_${RUNNAME}_RAW_QC -i ${PFBPROJECTNAME}' ALL FASTQ RAW' -z -f -o ${RAWQCDIR} ${RAWQCDIR}/fastqc/*${R1CODE}*.zip ${RAWQCDIR} ${RAWQCDIR}/fastqc/*${R2CODE}*.zip ${RAWQCDIR}/fastqscreen/*.txt"
conda deactivate

## Sending report to NC
### Zipping
NTHREADS=1
MEM=2
sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'zip_RAWQC' -e ${WORKDIR}'/zip_rawqc.e' -o ${WORKDIR}'/zip_rawqc.o' --wrap \
"cd ${WORKDIR} && echo ${PWD} && \
zip -r -9 --exclude='*.zip' ${WORKDIR}'/'${PFBPROJECTNAME}'_'${PFGPROJECTNAME}'_'${RUNNAME}'_RAW_QC.zip' `basename ${RAWQCDIR}`"
### Sending zip to NC
curl -u ${NCCRED} -T ${WORKDIR}'/'${PFBPROJECTNAME}'_'${PFGPROJECTNAME}'_'${RUNNAME}'_RAW_QC.zip' ${NCOUTDIR}

## TRIMMING
NTHREADS=5
MEM=10
mkdir -p ${TRIMDIR}'/fastp_report/' ${TRIMDIR}'/logs/' ${TRIMDIR}'/failed_reads/'
conda activate fastp
cd ${TRIMDIR}
for SAMPLENAME in ${SNAMES[*]}; do {
  echo ${SAMPLENAME};
  sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'trim.'${SAMPLENAME} -e ${TRIMDIR}'/logs/'${SAMPLENAME}'_trim.e' -o ${TRIMDIR}'/logs/'${SAMPLENAME}'_trim.o' --wrap \
  "fastp -i ${RAWDATADIR}/${SAMPLENAME}*_${R1CODE}*.f*.gz -I ${RAWDATADIR}/${SAMPLENAME}*_${R2CODE}*.f*.gz \
  -o ${TRIMDIR}/${SAMPLENAME}'_'${R1CODE}'_trim.fq.gz' -O ${TRIMDIR}/${SAMPLENAME}'_'${R2CODE}'_trim.fq.gz' \
  --unpaired1 ${TRIMDIR}/failed_reads/${SAMPLENAME}'_'${R1CODE}'_unp.fq.gz' \
	--unpaired2 ${TRIMDIR}/failed_reads/${SAMPLENAME}'_'${R2CODE}'_unp.fq.gz' \
  --failed_out ${TRIMDIR}/failed_reads/${SAMPLENAME}'_fail.fq.gz' \
	-z ${ZCOMP} --adapter_fasta ${ADAPTERSFASTA} --trim_poly_g \
  --cut_tail --cut_tail_mean_quality ${MINBASEQ} --length_required 36 --n_base_limit 3 --qualified_quality_phred ${MINBASEQ} -p \
  --thread ${NTHREADS} --html ${TRIMDIR}/fastp_report/${SAMPLENAME}'_fastp.html' --json ${TRIMDIR}/fastp_report/${SAMPLENAME}'_fastp.json'"
}
done
conda deactivate

## QC post-trim
NTHREADS=3
MEM=5
conda activate readsqc
mkdir -p ${TRIMDIR}'/fastqc' ${TRIMDIR}'/fastqscreen' ${TRIMDIR}'/multiqc' ${TRIMDIR}'/logs'
for SAMPLENAME in ${SNAMES[*]}; do {
	echo ${SAMPLENAME};
	TOOLEO='trimqc'
	sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J ${TOOLEO}'.'${SAMPLENAME} \
	-e ${TRIMDIR}'/logs/'${SAMPLENAME}'_'${TOOLEO}'.e' -o ${TRIMDIR}'/logs/'${SAMPLENAME}'_'${TOOLEO}'.o' --wrap \
	"fastqc -o ${TRIMDIR}/fastqc -t ${NTHREADS} ${TRIMDIR}/*${SAMPLENAME}*_${R1CODE}*_trim.fq.gz \
	${TRIMDIR}/*${SAMPLENAME}*_${R2CODE}*_trim*.fq.gz && \
  fastq_screen --threads ${NTHREADS} --force --outdir ${TRIMDIR}/fastqscreen --subset 100000 \
	--conf ${FASTQSCREENCONF} ${TRIMDIR}/*${SAMPLENAME}*_${R1CODE}*_trim.f*gz && \
	multiqc -n TRIM_${SAMPLENAME} -i ${PROJECTNAME}' '$SAMPLENAME' FASTQ TRIMMED' -z -f -o ${TRIMDIR}'/multiqc' ${TRIMDIR}/fastqc/*${SAMPLENAME}*.zip ${TRIMDIR}/fastqscreen/*${SAMPLENAME}* ${TRIMDIR}/fastp_report/*${SAMPLENAME}*.json"
}
done

NTHREADS=3
MEM=5
sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'TRIMQC_ALL' \
-e ${TRIMDIR}'/logs/'${PFBPROJECTNAME}'_mqc_trim.e' -o ${TRIMDIR}'/logs/'${PFBPROJECTNAME}'_mqc_trim.o' --wrap \
"multiqc -n ${PFBPROJECTNAME}_${RUNNAME}_TRIM_QC -i ${PFBPROJECTNAME}' ALL FASTQ TRIMMED' -z -f -o ${TRIMDIR} ${TRIMDIR}/fastqc/*.zip ${TRIMDIR}/fastp_report/*.json ${TRIMDIR}/fastqscreen/*.txt"
conda deactivate

### Zipping
NTHREADS=1
MEM=2
sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'zip_TRIMQC' -e ${WORKDIR}'/zip_trimqc.e' -o ${WORKDIR}'/zip_trimqc.o' --wrap \
"cd ${WORKDIR} && echo $PWD && \
zip -r -9 --exclude='*.zip' --exclude='*.gz' ${WORKDIR}'/'${PFBPROJECTNAME}'_'${PFGPROJECTNAME}'_'${RUNNAME}'_TRIM_QC.zip' `basename ${TRIMDIR}`"
### Sending zip to NC
curl -u ${NCCRED} -T $WORKDIR'/'${PFBPROJECTNAME}'_'${PFGPROJECTNAME}'_'${RUNNAME}'_TRIM_QC.zip' ${NCOUTDIR}


######

## SALMON COUNT through rna-count-salmon

### 01. Creating a user-specific tmp directory (needed at the moment, should be fixed later)
mkdir ${MYTMP}
### 02. Loading the env
#### The global one won't work due to trying to write data without permission, despite the --conda-prefix $MYTMP...
# conda activate /mnt/beegfs/pipelines/rna-count-salmon/env/
#### Thus, using my own git clone
conda activate rna-count-salmon
### 03. Creating and getting into the output directory
mkdir ${SALMONDIR}
cd ${SALMONDIR}
### 04. Creating the config and design files
NTHREADS=10
python3 ${RNA_COUNT_LAUNCHER} config --no-fastqc --no-multiqc --aggregate --threads ${NTHREADS} --design ${SALMONDIR}'/design.tsv' ${FASTA} ${GTF}
python3 ${RNA_COUNT_LAUNCHER} design -o ${SALMONDIR}'/design.tsv' ${TRIMDIR}
### 05. MANUAL EDITION of the design file (renaming samples)
### 06. Calling the pipeline in manual mode to use our custom tmp dir. **FASTQC WILL FAIL**, but we can run it manually afterwards
python3 ${RNA_COUNT_LAUNCHER} snakemake --snakemake-args ' --conda-prefix '${MYTMP}' -n --dag | dot -Tpng > dag.png '
python3 ${RNA_COUNT_LAUNCHER} snakemake --snakemake-args ' --conda-prefix '${MYTMP}' --attempt 5 '
### 06. Generating a multiqc report from salmon results
conda deactivate

conda activate readsqc
NTHREADS=3
MEM=2
sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'SALMON_ALL' \
-e ${SALMONDIR}'/logs/'${PFBPROJECTNAME}'_mqc_salmon.e' -o ${TRIMDIR}'/logs/'${PFBPROJECTNAME}'_mqc_salmon.o' --wrap \
"multiqc -n ${PFBPROJECTNAME}_${PFGPROJECTNAME}_${RUNNAME}_SALMON -i ${PROJECTNAME}' '${RUNNAME}' ALL SALMON COUNT' -z -f -o ${SALMONDIR} ${SALMONDIR}/pseudo_mapping/*"
conda deactivate
curl -u ${NCCRED} -T ${SALMONDIR}'/'${PFBPROJECTNAME}'_'${PFGPROJECTNAME}'_'${RUNNAME}'_SALMON.html' ${NCOUTDIR}


## IMMUNE INFILTRATE ESTIMATION through tpm-immunedecov (default for human tumor with all methods. Edit the YAML files for any other data type)

### 01. Cloning the pipeline locally
git clone 'https://github.com/tdayris/tpm-immunedecov' ${MYGITSREP}'/tpm-immunedecov'
### 02. Creating the output dir
mkdir -p ${IMMUNEDIR}'/logs'
cd ${IMMUNEDIR}
### 03. Converting the aggregated salmon counts table **WARNING** : USE THE "GENES" TABLE !
conda activate mamba
NTHREADS=3
MEM=2
sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'decovsalmon' \
-e ${IMMUNEDIR}'/logs/'${PFBPROJECTNAME}'_mqc_salmon.e' -o ${IMMUNEDIR}'/logs/'${PFBPROJECTNAME}'_mqc_salmon.o' --wrap \
"python3.9 ${MYGITSREP}/tpm-immunedecov/scripts/from_salmon_aggregate.py \
--output ${IMMUNEDIR}/${PFBPROJECTNAME}_salmon_immune_input.tsv ${SALMONDIR}/pseudo_mapping/aggregation/TPM.genes.tsv"
### 04. Running the pipeline
for CFILE in `ls ${IMMUNEDIR}/config*.yaml`; do {
	echo $CFILE;
	snakemake --use-conda --conda-frontend mamba --cores 1 --jobs ${NTHREADS} --configfile ${CFILE} \
	--snakefile ${MYGITSREP}'/tpm-immunedecov/Snakefile';
}
done
conda deactivate
### 04. Creating an archive and sending it to NextCloud
cd ${IMMUNEDIR}
NTHREADS=1
MEM=2
sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'DECOV_zip' \
-e ${IMMUNEDIR}'/logs/'${PFBPROJECTNAME}'_decov_zip.e' -o ${IMMUNEDIR}'/logs/'${PFBPROJECTNAME}'_decov_zip.o' --wrap \
"zip -r -9 ${WORKDIR}'/'${PFBPROJECTNAME}_immune.zip CIBERSORT/ CIBERSORT_ABS/ EPIC/ MCPcounter/ quanTIseq/ xCell/ ${PFBPROJECTNAME}_salmon_immune_input.tsv"
curl -u ${NCCRED} -T ${WORKDIR}'/'${PFBPROJECTNAME}'_immune.zip' ${NCOUTDIR}



## DGE
DGEDIR=${WORKDIR}'/DGE'
mkdir -p ${DGEDIR}
cd ${DGEDIR}
conda activate /mnt/beegfs/pipelines/rna-dge-salmon-deseq2/env
## Building config
python3 ${DGE_LAUNCHER} config --models 'Macrophage_type,Proinf,Basal,~Macrophage_type' 'F2020,Treated,Untreated,~F2020' 'Macrophage_type_F2020,Proinf_Treated,Basal_Treated,~Macrophage_type_F2020' 'Macrophage_type_F2020,Proinf_Untreated,Basal_Untreated,~Macrophage_type_F2020' --workdir ${DGEDIR} --threads 5 --columns Macrophage_type F2020 Macrophage_type_F2020 --no-multiqc ${GTF}
## Manually create the design file
## Generate the DAG plot
python3 ${DGE_LAUNCHER} snakemake --snakemake-args ' -n --dag | dot -Tpng > dag.png '
## Run DGE !
python3 ${DGE_LAUNCHER} snakemake --snakemake-args ' --attempt 5 '

### Creating an archive and sending it to NextCloud
NTHREADS=1
MEM=2
cd ${WORKDIR}
sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'DGE_zip' \
-e ${DGEDIR}'/logs/'${PFBPROJECTNAME}'_dge_zip.e' -o ${DGEDIR}'/logs/'${PFBPROJECTNAME}'_dge_zip.o' --wrap \
"zip -r -9 --exclude='*.RDS' --exclude=*.gtf ${WORKDIR}'/'${PFBPROJECTNAME}_dge.zip `basename ${DGEDIR}`"

curl -u ${NCCRED} -T ${WORKDIR}'/'${PFBPROJECTNAME}'_dge.zip' ${NCOUTDIR}









### Building an analysis for pairs from P29 & P30
# PROJECTNAME='C2J1_screening_PAIRS'
# WORKDIR=${OUTDIR}'/'${PROJECTNAME}
# SNAMEBLOCKS="1,2"  # PART(S) OF THE FASTQ FILENAMES THAT WILL FORM THE SAMPLENAME
# SALMON30='/mnt/beegfs/scratch/bioinfo_core/B20067_FRDA_05/data_output/P30_AUMA/CONCATENATED/SALMON/'
# SALMON29='/mnt/beegfs/scratch/bioinfo_core/B20067_FRDA_05/data_output/P29_AUMA/CONCATENATED/SALMON/'
# ls ${SALMON30}/pseudo_mapping/*creening*/quant.sf ${SALMON30}/pseudo_mapping/*C2J1*/quant.sf \
# ${SALMON29}/pseudo_mapping/13*creening*/quant.sf ${SALMON29}/pseudo_mapping/13*C2J1*/quant.sf > ${WORKDIR}'/FQpairs.txt'
#
# /mnt/beegfs/scratch/bioinfo_core/B20067_FRDA_05/data_output/P30_AUMA/CONCATENATED/SALMON//pseudo_mapping/*creening*/quant.sf
#
# ls ${FASTQDIR30}/*creening*trim.fq.gz ${FASTQDIR30}/*C2J1*trim.fq.gz > $WORKDIR'/FQlist.txt'
# ls ${FASTQDIR29}/13*creening*trim.fq.gz ${FASTQDIR29}/13*C2J1*trim.fq.gz >> $WORKDIR'/FQlist.txt'
#
# ## Manual editting of the FQlist to remove samples that do not belong to a pair
#
#
# FASTQDIR30='/mnt/beegfs/scratch/bioinfo_core/B20067_FRDA_05/data_output/P30_AUMA/CONCATENATED/TRIM'
# FASTQDIR29='/mnt/beegfs/scratch/bioinfo_core/B20067_FRDA_05/data_output/P29_AUMA/CONCATENATED/TRIM'
# TRIMFILES30=`ls ${FASTQDIR30}/*creening*trim.fq.gz ${FASTQDIR30}/*C2J1*trim.fq.gz`
# SNARRAY30=`echo ${TRIMFILES30} | xargs -n 1 basename | cut -d "_" -f ${SNAMEBLOCKS}`
# SNAMES30=`printf "%s\n" "${SNARRAY30[@]}" | sort -u`
# TRIMFILES29=`ls ${FASTQDIR29}/13*creening*trim.fq.gz ${FASTQDIR29}/13*C2J1*trim.fq.gz`
# SNARRAY29=`echo ${TRIMFILES29} | xargs -n 1 basename | cut -d "_" -f ${SNAMEBLOCKS}`
# SNAMES29=`printf "%s\n" "${SNARRAY29[@]}" | sort -u`
#
# SNAMES=("${SNAMES30[@]}" "${SNAMES29[@]}")
