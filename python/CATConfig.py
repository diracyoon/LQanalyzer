###################################################################                                                                         
###configure Job                                                                                                                                         ####################################################################                                                                                                           
###timeWait=1#                                                                                                                                           ###################################################                                                                                                      ### Make Input File                         

import os, getpass, sys,ROOT,time
from functions import *

from optparse import OptionParser


path_jobpre="/data1/"
if "tamsa2.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    path_jobpre="/data2/"




#Import parser to get options                                                                                                                                                  
parser = OptionParser()
parser.add_option("-p", "--period", dest="period", default="A",help="which data period or mc sample")
parser.add_option("-s", "--stream", dest="stream", default="NULL", help="Which data channel- ee,or mumu?")
parser.add_option("-j", "--jobs", dest="jobs", default=1, help="Name of Job")
parser.add_option("-J", "--setnjobs", dest="setnjobs", default=1, help="Name of Job")
parser.add_option("-c", "--cycle", dest="cycle", default="Analyzer", help="which cycle")
parser.add_option("-t", "--tree", dest="tree", default="ntuple/event", help="What is input tree name?")
#parser.add_option("-t", "--tree", dest="tree", default="event", help="What is input tree name?")                                                                              
parser.add_option("-o", "--logstep", dest="logstep", default=-1, help="How many events betwene log messages")
parser.add_option("-d", "--data_lumi", dest="data_lumi", default="A", help="How much data are you running on/ needed to weight mc?")
parser.add_option("-l", "--loglevel", dest="loglevel", default="INFO", help="Set Log output level")
parser.add_option("-n", "--nevents", dest="nevents", default=-1, help="Set number of events to process")
parser.add_option("-k", "--skipevent", dest="skipevent", default=-1, help="Set number of events to skip")
parser.add_option("-a", "--datatype", dest="datatype", default="NULL", help="Is data or mc?")
parser.add_option("-e", "--totalev", dest="totalev", default=-1, help="How many events in sample?")
parser.add_option("-x", "--xsec", dest="xsec", default=-1., help="How many events in sample?")
parser.add_option("-X", "--tagger", dest="tagger", default="123", help="random number string?")
parser.add_option("-T", "--targetlumi", dest="targetlumi", default=-1., help="How many events in sample?")
parser.add_option("-E", "--efflumi", dest="efflumi", default=-1., help="How many events in sample?")
parser.add_option("-O", "--outputdir", dest="outputdir", default="${LQANALYZER_DIR}/data/output/", help="Where do you like output to go?")
parser.add_option("-w", "--remove", dest="remove", default=True, help="Remove the work space?")
parser.add_option("-S", "--skinput", dest="skinput", default=True, help="Use SKTree as input?")
parser.add_option("-R", "--runevent", dest="runevent", default=True, help="Run Specific Event?")
parser.add_option("-N", "--useCATv742ntuples", dest="useCATv742ntuples", default=True, help="' to run on these samples")
parser.add_option("-L", "--LibList", dest="LibList", default="NULL", help="Add extra lib files to load")
parser.add_option("-D", "--debug", dest="debug", default=False, help="Run submit script in debug mode?")
parser.add_option("-m", "--useskim", dest="useskim", default="Lepton", help="Run submit script in debug mode?")
parser.add_option("-G", "--runtau", dest="runtau", default="runtau", help="Run fake mode for np bkg?")
parser.add_option("-P", "--runnp", dest="runnp", default="runnp", help="Run fake mode for np bkg?")
parser.add_option("-Q", "--runcf", dest="runcf", default="runcf", help="Run fake mode for np bkg?")
parser.add_option("-q", "--queue", dest="queue", default="queue", help="what queue")
parser.add_option("-v", "--catversion", dest="catversion", default="NULL", help="What cat version?")
parser.add_option("-f", "--skflag", dest="skflag", default="NULL", help="add input flag?")
parser.add_option("-b", "--usebatch", dest="usebatch", default="usebatch", help="Run in batch queue?")
parser.add_option("-F","--submitallfiles",dest="submitallfiles",default="False", help="force n=1000")


###################################################                                                                                                                            
#set the local variables using options                                                                                                                                         
###################################################                                                                                                                            
(options, args) = parser.parse_args()
number_of_cores = int(options.jobs)
sample = options.period
channel = options.stream
cycle = options.cycle
logstep = int(options.logstep)
loglevel = options.loglevel
runnp = options.runnp
runcf = options.runcf
tagger= options.tagger
### THESE ARE OPTIONS THAT CAN BE INCLUDED but not in example                                                                                                                  
tree = options.tree
number_of_events_per_job= int(options.nevents)
skipev = int(options.skipevent)
dataType = options.datatype
totalev = int(options.totalev)
xsec = float(options.xsec)
tar_lumi = float(options.targetlumi)
eff_lumi = float(options.efflumi)
data_lumi = options.data_lumi
catversion = options.catversion
Finaloutputdir = options.outputdir
remove_workspace=options.remove
useskinput=options.skinput
runevent= options.runevent
useCATv742ntuples = options.useCATv742ntuples
tmplist_of_extra_lib=options.LibList
DEBUG = options.debug
useskim = options.useskim
skflag = options.skflag
usebatch =options.usebatch

new_channel = channel.replace(":", "")
original_channel = new_channel
queuename = options.queue

if useskim == "SKTree_NoSkim":
    useskim ="NoCut"
elif useskim == "SKTree_LeptonSkim":
    useskim ="Lepton"
elif useskim == "SKTree_DiLepSkim":
    useskim="DiLep"
elif useskim == "SKTree_TriLepSkim":
    useskim="TriLep"

original_sample=sample
sample = sample.replace(":", " ")
datatype=""
splitsample  = sample.split()
if not len(splitsample)==1:
    sample = splitsample[0]
    for conf in range(1,len(splitsample)-1):
        if "nevents" in splitsample[conf]:
            conf+=1
            number_of_events_per_job = splitsample[conf]
        if "remove" in splitsample[conf]:
            conf+=1
            remove_workspace  = splitsample[conf]
        if "loglevel" in splitsample[conf]:
            conf+=1
            loglevel = splitsample[conf]
        if "cycle"  in splitsample[conf]:
            conf+=1
            cycle = splitsample[conf]
        if "njobs" in splitsample[conf]:
            conf+=1
            number_of_cores = splitsample[conf]
        if "skipevent" in splitsample[conf]:
            conf+=1
            skipev = splitsample[conf]
        if "skinput" in splitsample[conf]:
            conf+=1
            useskinput = splitsample[conf]
        if "runevent" in splitsample[conf]:
            conf+=1
            runevent = splitsample[conf]
        if "useCATv742ntuples" in splitsample[conf]:
            conf+=1
            useCATv742ntuples = splitsample[conf]
        if "catversion" in splitsample[conf]:
            conf+=1
            catversion = splitsample[conf]

running_batch=True

printToScreen=False
if str(usebatch) == "NULL":
    if str(running_batch) == "True":
        if printToScreen:
            print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            print "Running batch job:"
            print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    else:
        if printToScreen:
            print "Running standard root interactive job:"
else:
    if str(usebatch) == "False":
        running_batch=False
    elif  str(usebatch) == "false":
        running_batch=False



if not cycle == "SKTreeMaker":
    if not cycle == "SKTreeMakerNoCut":
        if not cycle == "SKTreeMakerDiLep":
            if not useskinput == "True":
                if not useskinput == "true":
                    if printToScreen:
                        print "You are running on FlatCATntuples. This will be more cpu extensive. This is only advisable if you are testing some new branches NOT in SKTrees."


output_mounted="/data2"
if "cmscluster.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    output_mounted="/data4"

##########################################################                                                                                                                     
### Make tmp directory for job                                                                                                                                                 
############################################################                                                                                                                   

tmpwork = output_mounted+"/CAT_SKTreeOutput/"+ getpass.getuser() + "/"
if not (os.path.exists(tmpwork)):
    os.system("mkdir " + tmpwork)

mc = len(sample)>1
if sample == "H_v2" or sample == "H_v3":
    mc= False

if mc:
    datatype="mc"
    if "D1" in sample:
        mc = 0
        datatype="data"
    if "D2" in sample:
        mc = 0
        datatype="data"
else:
    datatype="data"


if sample == "ALL":
    datatype="data"

if datatype == "mc":
    timeWait=10

if not dataType =="":
    datatype=dataType

##################################################################################################################                                                             
### Make a list of input samples: at the moment this is useless. Will add code to include * options in input                                                                   
##################################################################################################################                                                             
list = []
import re
if ("*" in sample) and mc:
    print "ADD code"
else:
    list.append(sample)

original_sample = sample
if useskinput == "true":
    if not mc:
        if useskim == "Lepton":
            new_channel="SK" + new_channel
        else:
            if useskim == "NoCut":
                new_channel="SK" + new_channel + "_nocut"
            else:
                if useskim == "DiLep":
                    new_channel="SK" + new_channel + "_dilep"
                else:
                    if useskim == "TriLep":
                        new_channel="SK" + new_channel + "_trilep"

    else:
        if useskim == "Lepton":
            sample="SK" + sample
        else:
            if useskim == "NoCut":
                sample="SK" + sample + "_nocut"
            else:
                if useskim == "DiLep":
                    sample="SK" + sample + "_dilep"
                else:
                    if useskim == "TriLep":
                        sample="SK" + sample + "_trilep"

elif useskinput == "True":

    if not mc:
        if useskim == "Lepton":
            new_channel="SK" + new_channel
        else:
            if useskim == "NoCut":
                new_channel="SK" + new_channel + "_nocut"
            else:
                if useskim == "DiLep":
                    new_channel="SK" + new_channel + "_dilep"
                else:
                    if useskim == "TriLep":
                        new_channel="SK" + new_channel + "_trilep"

    else:
        if useskim == "Lepton":
            sample="SK" + sample
        else:
            if useskim == "NoCut":
                sample="SK" + sample + "_nocut"
            else:
                if useskim == "DiLep":
                    sample="SK" + sample + "_dilep"
                else:
                    if useskim == "TriLep":
                        sample="SK" + sample + "_trilep"

if printToScreen:
    print "Input sample = " + sample
if not mc:
    if printToScreen:
        print "Input channel = " + new_channel



##############################################################################################                                                                                 
#### Check if sktrees are located on current machines  (not used when running on cmsX at snu)                                                                                  
#############################################################################################                                                                                  

isfile = os.path.isfile
join = os.path.join
import platform

if platform.system() != "Linux":

    localDir = os.getenv("LQANALYZER_DIR")+ "/data/input/"
    if not mc:
        localDir = os.getenv("LQANALYZER_DIR")+ "/data/input/data/" + new_channel  + sample
    else:
        localDir = os.getenv("LQANALYZER_DIR")+ "/data/input/mc/"  + sample


    if not os.path.exists(localDir):
        print "No files in current location: Will copy them over"
        CopySKTrees(new_channel,sample,mc,"True")
    elif  sum(1 for item in os.listdir(localDir) if isfile(join(localDir, item))) == 0:
        print "No files are located locally: Will copy from cms21 machine"
        CopySKTrees(new_channel,sample,mc,"True")
    else:
        update = raw_input("Files already located on current machine. Do you want these updating? Yes/No")
        if update == "Yes":
            print "Updating local sktree"
            CopySKTrees(new_channel,sample,mc,"True")
        elif update == "yes":
            print "Updating local sktree"
            CopySKTrees(new_channel,sample,mc,"True")
        else:
            CheckPathInFile(new_channel,sample,mc)


##################################################################################################################                                                             
#Find the DS name (and lumi if MC) from txt/datasets.txt                                                                                                                       
##################################################################################################################                                                             
inDS = ""
mcLumi = 1.0
filechannel=""

catversions = ["v7-6-4",
               "v7-6-3",
               "v7-6-2",
               "v7-4-5",
               "v7-4-4"]

sample_catversion = ""
output_catversion=os.getenv("CATVERSION")

#### Check latest tag/version for DS.                                                                                                                                          
iversion=0
while inDS == "":
    if platform.system() == "Linux":
        version="_CAT"
        sample_catversion = catversions[iversion]


        if catversion != "":
            sample_catversion = catversion
            output_catversion = catversion

        if printToScreen:
            print "Using CAT " +sample_catversion + " ntuples"
        if mc:
            filename = os.getenv("LQANALYZER_RUN_PATH") + '/txt/datasets_snu_CAT_mc_' +sample_catversion +  '.txt'
            if "cmscluster.snu.ac.kr" in str(os.getenv("HOSTNAME")):
                filename = os.getenv("LQANALYZER_RUN_PATH") + '/txt/Cluster/datasets_snu_cluster_CAT_mc_' +sample_catversion +  '.txt'

        else:
            filename = os.getenv("LQANALYZER_RUN_PATH") + '/txt/datasets_snu_CAT_data_'  +sample_catversion +'.txt'
            if "cmscluster.snu.ac.kr" in str(os.getenv("HOSTNAME")):
                filename = os.getenv("LQANALYZER_RUN_PATH") + '/txt/Cluster/datasets_snu_cluster_CAT_data_'  +sample_catversion +'.txt'
    else:
        filename = os.getenv("LQANALYZER_RUN_PATH") + 'txt/datasets_mac.txt'

    if printToScreen:
        print "Using " + filename
    if not mc:
        if printToScreen:
            print "Running on data "
            print new_channel + " " + sample
        if not (os.path.exists(filename)):
            iversion = iversion +1
            continue
        for line in open(filename, 'r'):
            if not line.startswith("#"):
                entries = line.split()
                if len(entries)==3:
                    if new_channel ==entries[0] and sample == entries[1]:
                        inDS = entries[2]
        sample = "period"+sample
        eff_lumi=1.
        tar_lumi=1.
        filechannel = new_channel+"_"
    else:
        if printToScreen:
            print "Running on MC"
        if not (os.path.exists(filename)):
            iversion = iversion +1
            continue
        for line in open(filename, 'r'):
            if not line.startswith("#"):
                entries = line.split()
                if len(entries)==6:
                    if sample == entries[0]:
                        eff_lumi = entries[4]
                        inDS = entries[5]
    iversion = iversion +1
    if inDS == "":
        if catversion != "":
            print "Input dataset is not available in specifies catversion: Exiting"
            sys.exit()
        print "LQAnalyzer :: WARNING :: Sample is not available in " + filename + "."
        print "Will look in previous compatable version. Need input from user if this is ok or a mistake."
        update = raw_input("This is likely because you have not changed the name of the input file. Since CATVERSION v7-6-3 these were changed. Is using anolder version of tuples ok for what you are doing? If you wish to use an older sample type Y. If not change input. run 'sktree -l' for options" )
        if not  update == "Y":
            sys.exit()
        if iversion == len(catversions):
            print "LQAnalyzer :: ERROR :: Input dataset is not available: Exiting"
            sys.exit()


InputDir = inDS

if  not os.path.exists(InputDir):
    print "No directory " + InputDir
    sys.exit()


isfile = os.path.isfile
join = os.path.join
number_of_files = sum(1 for item in os.listdir(InputDir) if isfile(join(InputDir, item)))

path_log=path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/Cluster_"+original_sample + tagger + ".log"

time.sleep(5.)
os.system("qstat -u '*' > " + path_log)
file_log=open(path_log,"r")

njobs=0

for line in file_log:
    splitline  = line.split()
    if len(splitline) < 6:
        continue
    if not "prior" in line:
        njobs=njobs+1
file_log.close()

if os.path.exists(path_log):
    os.system("rm " + path_log)

if printToScreen:
    print "number_of_files = " + str(number_of_files) + " njobs running in batch = " + str(njobs)

njobs_max_allowed = 1750
if queuename == "long":
    njobs_max_allowed = 2750

if number_of_files < (njobs_max_allowed- njobs):
    if printToScreen:
        print "Job is running in background............"
    os.system("mkdir " + tmpwork + "/" + tagger)



