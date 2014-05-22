#!/bin/sh

runall=true

if [[ $runall  == "true" ]];
then
    source functions.sh
    cycle="ElectronCF"
    skinput="True"
    useskim="DiLep"
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    
    
    declare -a input_samples=("DY10to50" "DY50plus" "ttbar")
        
    stream="egamma"
    outputdir=$LQANALYZER_DIR"/data/output/CFElectron/"
    
 
    ### submit this configured job (uses bin/submit.sh)
    source submit.sh $1
fi




echo ""
echo "End of example_submit.sh script."