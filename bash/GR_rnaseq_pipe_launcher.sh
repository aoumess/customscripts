## NEXTCLOUD
# FLUSER='b_job'
FLUSER=${USER}

## PARAMETERS
FRCODE='R1'
RRCODE='R2'
MINBASEQ=20  # MINIMUM BASE QUALITY FOR FILTERING
ZCOMP=6
## FastqScreen
FASTQSCREENCONF='/mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/fastq_screen.conf'
## Fastp
ADAPTERSFASTA="/mnt/beegfs/userdata/${FLUSER}/RESOURCES/ADAPTERS/CONCAT.fa"
MYTMP="/mnt/beegfs/userdata/${FLUSER}/mytmp/"
MYGITSREP="/mnt/beegfs/userdata/${FLUSER}/gits/"

## PATHS

### PROJECT (BIGR)
PFBPROJECTNAME='B24033_LULA_01'

### RUN01
PFGPROJECTNAME='RT18822_RNASeq_Batch62'
RUNNAME='240426_A00461_0558_BH7FV3DRX5'

## MAIN DIRS
### WHERE THE RAW FASTQ FILES SHOULD END to init the analysis
INDIR="/mnt/beegfs/scratch/bioinfo_core/${PFBPROJECTNAME}/data_input/"
RAWDATADIR="${INDIR}/${PFGPROJECTNAME}/${RUNNAME}/"
### WORKING DIRECTORY
OUTDIR="/mnt/beegfs/scratch/bioinfo_core/${PFBPROJECTNAME}/data_output/"
WORKDIR="${OUTDIR}/${PFGPROJECTNAME}/"

RAWQCDIR=${WORKDIR}'/RAW_QC'
SALMONDIR=${WORKDIR}'/SALMON'
IMMUNEDIR=${WORKDIR}'/IMMUNE'

## Nextcloud (requires a 'nextcloud_credentials.txt' file in the user home, set as chmod 700, which contains the user password)
NCUSER=`echo ${FLUSER} | tr '[:lower:]' '[:upper:]'`
source "/home/${FLUSER}@intra.igr.fr/nextcloud_credentials.txt"
NCCRED="${NCUSER}:${NCPWD}"
NCROOTDIR="https://nextcloud.gustaveroussy.fr/remote.php/dav/files/${NCUSER}/transfert_bigr/"
NCOUTDIRA=${NCROOTDIR}'/'${PFBPROJECTNAME}
curl --user ${NCCRED} -X MKCOL $NCOUTDIRA
NCOUTDIR=${NCOUTDIRA}'/'${RUNNAME}
curl --user ${NCCRED} -X MKCOL $NCOUTDIRB


## GENERATING SAMPLENAMES
SNARRAY=(`ls ${RAWDATADIR}/*.f*q* | xargs -n 1 basename | sed  's/_S.*.f*q*//'`)
SNAMES=(`printf "%s\n" "${SNARRAY[@]}" | sort -u`)
echo ${SNAMES[@]}

## PER-RUN RAW QC
NTHREADS=2
MEM=5
conda activate readsqc
mkdir -p ${RAWQCDIR}'/fastqc' ${RAWQCDIR}'/multiqc' ${RAWQCDIR}'/logs'
cd ${RAWQCDIR}
for SAMPLENAME in ${SNAMES[@]}; do {
	echo ${SAMPLENAME};
	TOOLEO='rawqc'
	sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J ${TOOLEO}'.'${SAMPLENAME} -e ${RAWQCDIR}'/logs/'${SAMPLENAME}'_'${TOOLEO}'.e' -o ${RAWQCDIR}'/logs/'${SAMPLENAME}'_'${TOOLEO}'.o' --wrap "fastqc -o ${RAWQCDIR}/fastqc -t ${NTHREADS} ${RAWDATADIR}/*${SAMPLENAME}*_${FRCODE}*.f*q* ${RAWDATADIR}/*${SAMPLENAME}*_${RRCODE}*.f*q*"
}
done

NTHREADS=3
MEM=5
sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J 'mqc_RAWQC' -e ${RAWQCDIR}'/logs/'${PFBPROJECTNAME}'_mqc_raw.e' -o ${RAWQCDIR}'/logs/'${PFBPROJECTNAME}'_mqc_raw.o' --wrap \
"multiqc -n ${PFBPROJECTNAME}_${RUNNAME}_RAW_QC -i ${PFBPROJECTNAME}' ALL FASTQ RAW' -z -f -o ${RAWQCDIR} ${RAWQCDIR}/fastqc/*${FRCODE}*.zip ${RAWQCDIR} ${RAWQCDIR}/fastqc/*${RRCODE}*.zip ${RAWQCDIR}/fastqscreen/*.txt"
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

## QC post-trim
NTHREADS=3
MEM=5
conda activate readsqc
mkdir -p ${TRIMDIR}'/fastqc' ${TRIMDIR}'/fastqscreen' ${TRIMDIR}'/multiqc' ${TRIMDIR}'/logs'
for SAMPLENAME in ${SNAMES[@]}; do {
	echo ${SAMPLENAME};
	TOOLEO='trimqc'
	sbatch --mem=${MEM}'G' --cpus-per-task=${NTHREADS} --tasks=1 --nodes=1 -J ${TOOLEO}'.'${SAMPLENAME} \
	-e ${TRIMDIR}'/logs/'${SAMPLENAME}'_'${TOOLEO}'.e' -o ${TRIMDIR}'/logs/'${SAMPLENAME}'_'${TOOLEO}'.o' --wrap \
	"fastqc -o ${TRIMDIR}/fastqc -t ${NTHREADS} ${TRIMDIR}/*${SAMPLENAME}*_${FRCODE}*_trim.fq.gz \
	${TRIMDIR}/*${SAMPLENAME}*_${RRCODE}*_trim*.fq.gz && \
  fastq_screen --threads ${NTHREADS} --force --outdir ${TRIMDIR}/fastqscreen --subset 100000 \
	--conf ${FASTQSCREENCONF} ${TRIMDIR}/*${SAMPLENAME}*_${FRCODE}*_trim.f*gz && \
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
"zip -r -9 --exclude='*.RDS' --exclude=*.gtf --exclude=*.out --exclude=*.err ${WORKDIR}'/'${PFBPROJECTNAME}_dge.zip `basename ${DGEDIR}`"

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
# ls ${INDIR30}/*creening*trim.fq.gz ${INDIR30}/*C2J1*trim.fq.gz > $WORKDIR'/FQlist.txt'
# ls ${INDIR29}/13*creening*trim.fq.gz ${INDIR29}/13*C2J1*trim.fq.gz >> $WORKDIR'/FQlist.txt'
#
# ## Manual editting of the FQlist to remove samples that do not belong to a pair
#
#
# INDIR30='/mnt/beegfs/scratch/bioinfo_core/B20067_FRDA_05/data_output/P30_AUMA/CONCATENATED/TRIM'
# INDIR29='/mnt/beegfs/scratch/bioinfo_core/B20067_FRDA_05/data_output/P29_AUMA/CONCATENATED/TRIM'
# TRIMFILES30=`ls ${INDIR30}/*creening*trim.fq.gz ${INDIR30}/*C2J1*trim.fq.gz`
# SNARRAY30=`echo ${TRIMFILES30} | xargs -n 1 basename | cut -d "_" -f ${SNAMEBLOCKS}`
# SNAMES30=`printf "%s\n" "${SNARRAY30[@]}" | sort -u`
# TRIMFILES29=`ls ${INDIR29}/13*creening*trim.fq.gz ${INDIR29}/13*C2J1*trim.fq.gz`
# SNARRAY29=`echo ${TRIMFILES29} | xargs -n 1 basename | cut -d "_" -f ${SNAMEBLOCKS}`
# SNAMES29=`printf "%s\n" "${SNARRAY29[@]}" | sort -u`
#
# SNAMES=("${SNAMES30[@]}" "${SNAMES29[@]}")
