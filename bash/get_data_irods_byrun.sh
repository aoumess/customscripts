## Retrieve data from the warehouse using iRODS.
## Data are stored in the default folder on flamingo, and organised in subfolders by run name, then by dataset_id

PFBID='B23019_BAPI_01'
PFGID='ICARUS'
# PFGRUN='230324_A00461_0426_AHTNTVDSX5'
# PFGRUN='230329_A00461_0427_BHK3C5DMXY'
PFGRUN='230405_A00461_0430_AH25J7DRX3'
DATADIR="/mnt/beegfs/scratch/bioinfo_core/${PFBID}/data_input/${PFGRUN}/"
PROTOCOL="%RNA%"

## Initialize the iRods connexion
iinit
## Get list of available datasets [PFB TYPE]
COLLECLIST=(`imeta qu -C 'projectName' like ${PFGID} and 'protocol' like ${PROTOCOL} | grep -v '\-\-\-\-' | awk -F ': ' '{print $2}'`)
## Get list of files per dataset_id
mkdir -p ${DATADIR}
echo -n > "${DATADIR}/${PFGRUN}_filelist"
for CL in ${COLLECLIST[@]}; do {
  DSETID=`basename ${CL}`
  echo ${DSETID}
  echo ${DSETID} >> "${DATADIR}/${PFGRUN}_filelist"
  mkdir -p "${DATADIR}/${DSETID}"
  FILELIST=(`ils "${CL}/archive" | grep -v ':' | cut -d' ' -f 3`)
  for FL in ${FILELIST[@]}; do {
    echo "Retrieving ${FL} ...";
    echo ${FL} >> "${DATADIR}/${PFGRUN}_filelist"
    iget -fvK "${CL}/archive/${FL}" "${DATADIR}/${DSETID}/";
  }
  done;
}
done;

for CL in ${COLLECLIST[@]}; do {
  DSETID=`basename ${CL}`
  echo ${DSETID}
  echo ${DSETID} >> "${DATADIR}/${PFGRUN}_filelist"
  mkdir -p "${DATADIR}/${DSETID}"
  FILELIST=(`ils "${CL}/archive" | grep -v ':' | cut -d' ' -f 3`)
  for FL in ${FILELIST[@]}; do {
    echo "Retrieving ${FL} ...";
    echo ${FL} >> "${DATADIR}/${PFGRUN}_filelist"
    # iget -fvK "${CL}/archive/${FL}" "${DATADIR}/${DSETID}/";
  }
  done;
}
done;
