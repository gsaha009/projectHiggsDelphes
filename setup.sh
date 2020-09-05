echo " ... Setting Delphes-3.4.2 environment ... "
APPDIR=/home/gsaha
WORKDIR=$APPDIR/Work/HiggsProject_Pheno
cd $WORKDIR
export LD_LIBRARY_PATH=/home/gsaha/Packages/Delphes-3.4.2/:$LD_LIBRARY_PATH
echo libPath=$LD_LIBRARY_PATH
