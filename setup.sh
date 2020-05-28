APPDIR=/afs/cern.ch/work/g/gsaha/public
WORKDIR=$APPDIR/projectHiggsDelphes
cd $APPDIR/CMSSW_8_0_29/src
echo 'CMSSW_8_0_29'
eval `scramv1 runtime -sh`
echo 'All in CMSSW Environment'
cd $WORKDIR
export LD_LIBRARY_PATH=/afs/cern.ch/work/g/gsaha/public/Delphes-3.4.1/:$LD_LIBRARY_PATH
echo libPath=$LD_LIBRARY_PATH