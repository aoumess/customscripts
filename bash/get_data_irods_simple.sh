PFBID='B22070_PAHO_01'
DATADIR="/mnt/beegfs/scratch/bioinfo_core/${PFBID}/data_input/${RUN}/"

## Initialize the iRods connexion
iinit
## Get list of available datasets
COLLECLIST=(`imeta qu -C 'projectName' like ${PFBID} | grep -v '\-\-\-\-' | awk -F ': ' '{print $2}'`)
## Get list of files per dataset
mkdir -p ${DATADIR}
for CL in ${COLLECLIST[@]}; do {
  echo ${CL}
  DSETID=`basename ${CL}`
  echo ${DSETID}
  echo ${DSETID} >> "${DATADIR}/filelist"
  mkdir -p "${DATADIR}/${DSETID}"
  FILELIST=(`ils "${CL}/archive" | grep -v ':' | cut -d' ' -f 3`)
  for FL in ${FILELIST[@]}; do {
    echo "Retrieving ${FL} ...";
    echo ${FL} >> "${DATADIR}/filelist"
    iget -fvK "${CL}/archive/${FL}" "${DATADIR}/${DSETID}/";
  }
  done;
}
done;
