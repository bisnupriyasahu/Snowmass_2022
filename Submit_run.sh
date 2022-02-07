#!/bin/sh                                                                                                                                                                                                
# 

# untar your Python installation. Make sure you are using the right version!
tar -xzf python##.tar.gz
# (optional) if you have a set of packages (created in Part 1), untar them also
tar -xzf packages.tar.gz

export PATH=$PWD/python/bin:$PATH
export PYTHONPATH=$PWD/packages
export HOME=$PWD

echo "input parameters: cluster, process" $1 $2

export pwd=/nfs_scratch/bsahu/Snowmass_2022/Post_Analyzr/CMSSW_12_2_0/delphes
source /cvmfs/cms.cern.ch/cmsset_default.sh
CLUSTER=$1
PROCESS=$2
echo "PROCESS=    " $2
export $SCRAM_ARCH=slc7_amd64_gcc700
cd CMSSW_12_2_0/src/delphes

# make sure the script will use your Python installation, 
# and the working directory as its home location

#RUNPATH=/nfs_scratch/bsahu/Snowmass_2022/Post_Analyzr/CMSSW_12_2_0/delphes
#RUNPATH=/nfs_scratch/bsahu/Snowmass_2022/Post_Analyzr/CMSSW_12_2_0/delphes#
#cd $RUNPATH


#########  Smaple/Job splitting                                                                                                                                                                           
echo "path is PWD " $PWD
echo "path is pwd " $pwd

SplitingNumber=10
DataSetArray=($(cat InputSamples_DATA.txt)) # array of the input datadets                                                                                                                                
echo "DataSetArray is " $DataSetArray
echo "process" $PROCESS                                                                                                                                                                                  
#echo ${DataSetArray[$PROCESS / $SplitingNumber]}                                                                                                                                                         
#echo "DATASETNAME:" ${DataSetArray[$PROCESS / $SplitingNumber]}
DataSetName=${DataSetArray[$PROCESS / $SplitingNumber]}
rootNumber=$(($PROCESS % $SplitingNumber))
echo ${DataSetArray[$PROCESS / $SplitingNumber]}
echo "root Number is " $(($PROCESS % $SplitingNumber))

########### complie the Skimmer                                                                                                                                                                           
#make


########### loop over all root file in a dataset directory                                                                                                                                                
#xrdfs root://cmsxrootd.hep.wisc.edu ls "/hdfs/store/"$DataSetName | grep .root | while read FullDataSetName
#ls "/hdfs/store/user/bsahu/"$DataSetName | grep $rootNumber.root | while read FullDataSetName
#ls "/hdfs/store/user/varuns/"$DataSetName | grep $rootNumber.root | while read FullDataSetName
#xrdfs "root://cmsxrootd.hep.wisc.edu/store/user/varuns/"$DataSetName | grep $rootNumber.root | while read FullDataSetName
#xrdfs "root://eoscms.cern.ch//"$DataSetArray | grep $rootNumber.root | while read FullDataSetArray

#ls "$DataSetName" | grep .root

#for FILE in $DataSetName;
#for FILE in *; do echo $DataSetName/* ; while read FullDataSetName
#for FILE in *; do cat $FILE; done
xrdfs "root://cmseos.fnal.gov/" ls $DataSetArray | grep $rootNumber.root

#xrdfs ls $DataSetArray | grep $rootNumber.root | while read FullDataSetArray

xrdfs "root://cmseos.fnal.gov/" ls $DataSetArray | grep $rootNumber.root | while read FullDataSetName

############  Here is where the Skimmer is running     ############                                                                                                                                       
do
 echo "full data set name "$FullDataSetName
 file=`echo root://cmseos.fnal.gov/$DataSetArray$FullDataSetName`
 echo "File name is "$file
 ShortName=$file
 echo "ShortName= "$ShortName
 # This removes all the string before Moriond17 (including Moriond17)                                                                                                                                     
 #ShortName=${file##*Moriond17}  # This removes all the string before Moriond17 (including Moriond17)                                                                                                     
 echo "Here is the short Name   ------>" $ShortName

 python3 PostAN.py $ShortName  "Out_"$OutName$rootNumber.root
done
############  Here is where the Skimmer ends          ############

IFS="/"                                                                                                                                                                                                   
set $DataSetName
echo "argument:1 :" $1 
echo "argument:2 :" $2 
echo "argument: 3 :" $3
echo "argument: 4 :" $4 
echo "argument: 5 :" $5 
echo "argument: 6 :" $6 
echo "argument:7 :" $7 
echo "argument:8 :" $8
echo "argument:9 :" $9 
echo "argument:10 :" $10 
echo "argument:11 :" $11
#OutName=$4"_"$5"_"$6"_"$rootNumber".root"  # this makes the 4th and 6th pieces of the                                                                                                         
OutName="Snowmass_"$7"_"$rootNumber".root"  # this makes the 4th and 6th pieces of the                                                                                                         
echo "Outname is " $OutName

hadd -f $OutName "Out_"$OutName$rootNumber*.root


##########  remove the unneccesat files                                                                                                                                                                   
rm Out_*root
echo "Done execution ..."
cd ${_CONDOR_SCRATCH_DIR}
rm -rf ${3}

exit $exitcode
echo "DONE!"
