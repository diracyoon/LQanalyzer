#!/bin/sh
declare -a input_samples=("")
stream=""
cycle=""
njobs=1
data_lumi=""
LibList="NULL"
outputdir=${LQANALYZER_DIR}/data/output/
loglevel="INFO"
remove_workspace="True"
logstep=-1
queue="NONE" 
skipevent=-1
nevents=-1
totalev=-1
xsec=-1.
targetlumi=-1.
efflumi=-1.
remove="True"
skinput="True"
useCATv742ntuples="True"
runevent=-1
DEBUG="False"
useskim="Lepton"
runnp="False"
runcf="False"
runtau="False"
skflag="NULL"
catversion="NULL"
usebatch="NULL"
run_in_bkg="False"
drawhists="False"
setnjobs="False"
submitallfiles="False"
function makeParseVariable() {
    if ! [ -z ${2} ] 
        then
        local myresult="-"${1}" "${2}
        echo "${myresult}"
    fi
}

function checkDir() {
    ls ${1} > /dev/null
    if [ $? -ne 0 ]
        then
        echo "ERROR:   "${1}" not found!"
    fi
}

function checkAndMakeDir() {
    echo -ne "Checking for directory "${1}
    ls ${1} > /dev/null
    if [ $? -ne 0 ]
        then
        mkdir ${1}
        echo "Made directory "${1}
        echo ""
    else
        echo "... found."
    fi  
}

function setupLQanalyzer() {
    if [ -z "$LQANALYZER_DIR" ]; then
	cd ../
	source setup.sh
	cd -
    fi
}