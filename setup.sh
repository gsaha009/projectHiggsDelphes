echo "Setting Delphes-3.4.2 environment ===>>>"
APPDIR=/home/gsaha
WORKDIR=$APPDIR/DelphesHiggsProject
cd $WORKDIR
#export LD_LIBRARY_PATH=/home/sinpcms/DEBABRATA/Packages/Pheno_Packages/Delphes-3.4.1/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/sinpcms/GOURAB/Delphes-3.4.2/:$LD_LIBRARY_PATH
echo libPath=$LD_LIBRARY_PATH
