PFGID='P32_NAIB'
PFBID='B21031_JEMI_01'
RUNNAME='210715_A00461_0214_BHG3KKDRXY'
DATADIR='/mnt/beegfs/scratch/bioinfo_core/'${PFBID}'/data_input/'${PFGID}'/'${RUNNAME}

## Initialize the iRods connexion
iinit
## Get list of available datasets
COLLECLIST=(`imeta qu -C 'projectName' like ${PFGID} and 'datasetName' = ${RUNNAME} | grep -v '\-\-\-\-' | awk -F ': ' '{print $2}'`)
## Get list of files per dataset
mkdir -p ${DATADIR}
for CL in ${COLLECLIST[@]}; do {
  echo ${CL}
  FILELIST=(`ils ${CL}'/archive' | grep -v ':' | cut -d' ' -f 3`)
  for FL in ${FILELIST[@]}; do {
    iget -vK ${CL}'/archive/'${FL} ${DATADIR};
  }
  done;
}
done;
