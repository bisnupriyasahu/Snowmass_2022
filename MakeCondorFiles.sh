#!/bin/sh
#export pwd=/nfs_scratch/sacharya/UL_studies/CMSSW_10_6_20/src
cat>Job_${4}.sh<<EOF
#!/bin/sh  
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_10_6_26/src
eval `scramv1 runtime -sh`
cd \${_CONDOR_SCRATCH_DIR}                                                                                                                                                                                                                   
python3 ${1} ${2} ${3}                                                                                                                                                                                                               
EOF

chmod 775 Job_${4}.sh

cat>condor_${4}<<EOF
x509userproxy =/tmp/x509up_u10044
executable = ./Job_${4}.sh
notification         = never
aswhenToTransferOutput = On_Exit
shouldTransferFiles  = yes
#requirements = (TARGET.UidDomain == "hep.wisc.edu" && TARGET.HAS_CMS_HDFS)  
requirements = (OpSysAndVer == "CENTOS7" && TARGET.Arch == "X86_64" && (MY.RequiresSharedFS=!=true || TARGET.HasAFS_OSG) && (TARGET.HAS_OSG_WN_CLIENT =?= TRUE || TARGET.IS_GLIDEIN=?=true) && IsSlowSlot=!=true)
on_exit_remove       = (ExitBySignal == FALSE && (ExitCode == 0 || ExitCode == 42 || NumJobStarts>3))
+IsFastQueueJob      = True
getenv = true 
request_memory       = 1992
request_disk         = 2048000
transfer_input_files = ${1}, ${5}, ${6}
output               = \$(Cluster)_\$(Process)_${4}.out 
error                = \$(Cluster)_\$(Process)_${4}.err
Log                  = \$(Cluster)_\$(Process)_${4}.log
Queue
EOF

condor_submit condor_${4}
