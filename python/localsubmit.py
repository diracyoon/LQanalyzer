##################################################################
###configure Job
####################################################################
###timeWait=1#

###################################################
### Make Input File
####################################################
import os, getpass, sys,ROOT
from functions import *

from optparse import OptionParser

#Import parser to get options
parser = OptionParser()
parser.add_option("-p", "--period", dest="period", default="A",help="which data period or mc sample")
parser.add_option("-s", "--stream", dest="stream", default="NULL", help="Which data channel- ee,or mumu?")
parser.add_option("-j", "--jobs", dest="jobs", default=1, help="Name of Job")
parser.add_option("-c", "--cycle", dest="cycle", default="Analyzer", help="which cycle")
parser.add_option("-t", "--tree", dest="tree", default="ntuple/event", help="What is input tree name?")
#parser.add_option("-t", "--tree", dest="tree", default="event", help="What is input tree name?")
parser.add_option("-o", "--logstep", dest="logstep", default=-1, help="How many events betwene log messages")
parser.add_option("-d", "--data_lumi", dest="data_lumi", default="A", help="How much data are you running on/ needed to weight mc?")
parser.add_option("-l", "--loglevel", dest="loglevel", default="INFO", help="Set Log output level")
parser.add_option("-n", "--nevents", dest="nevents", default=-1, help="Set number of events to process")
parser.add_option("-k", "--skipevent", dest="skipevent", default=-1, help="Set number of events to skip")
parser.add_option("-a", "--datatype", dest="datatype", default="", help="Is data or mc?")
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
parser.add_option("-P", "--runnp", dest="runnp", default="runnp", help="Run fake mode for np bkg?")
parser.add_option("-G", "--runtau", dest="runtau", default="runtau", help="Run fake mode for np bkg?")
parser.add_option("-Q", "--runcf", dest="runcf", default="runcf", help="Run fake mode for np bkg?")
parser.add_option("-q", "--queue", dest="queue", default="", help="Which queue to use?")
parser.add_option("-v", "--catversion", dest="catversion", default="NULL", help="What cat version?")
parser.add_option("-f", "--skflag", dest="skflag", default="NULL", help="add input flag?")
parser.add_option("-b", "--usebatch", dest="usebatch", default="usebatch", help="Run in batch queue?")
parser.add_option("-J", "--setnjobs", dest="setnjobs", default="False", help="user sets njobs?")
parser.add_option("-F","--submitallfiles",dest="submitallfiles",default="False", help="force n=1000")
parser.add_option("-g","--tmpfilename",dest="tmpfilename",default="", help="tmpname")


###################################################
#set the local variables using options
###################################################
(options, args) = parser.parse_args()
number_of_cores = int(options.jobs)
setjobs = options.setnjobs
setnumber_of_cores=False
if setjobs== "true":
    setnumber_of_cores=True

sample = options.period
channel = options.stream
cycle = options.cycle
logstep = int(options.logstep)
loglevel = options.loglevel
runnp = options.runnp
runtau = options.runtau
runcf = options.runcf
queue = options.queue
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

tmp_filename=options.tmpfilename
if tmp_filename == "None":
    tmp_filename=""
    

tmpsubmit_allfiles=options.submitallfiles
submit_allfiles=False
if tmpsubmit_allfiles == "true":
    submit_allfiles=True




###### New for 801.7 tag  
ClusterStatFile = os.getenv("LQANALYZER_DIR")+ "/python/CheckCluster.py"
os.system("python " + ClusterStatFile + " -x " + tagger + " -s " + sample)


if not skflag == "":
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    print "skflag = " + skflag 

new_channel = channel.replace(":", "")

original_channel = new_channel

path_jobpre="/data1/"
if "tamsa2.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    path_jobpre="/data2/"
    
    
    
queuepath=path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CattupleConfig/QUEUE/ForceQueue.txt"
file_queuepath = open(queuepath,"r")
for line in file_queuepath:
    if "#" in line:
        continue
    sline = line.split()
    for s in sline:
        if "@" in s or "fast" in s or "all" in s or "long" in s:
            queue= s
file_queuepath.close()


if queue == "None":
    queue = "fastq"

queue_command = ''
qlist=[]

queueok=False
if queue:
    queue_command=' -q '
if queue == "allq":
    queue_command=''
    queueok=True
if queue == "fastq":
    queue_command=queue_command+ ' fastq '
    queueok=True
if queue == "longq":
    queue_command=queue_command+ ' longq '
    queueok=True
queue=queue.replace("all.q@cms-0-","node")
queue=queue.replace(".local","")
if "node1" in queue:
    qlist.append(1)
if "node2" in queue:
    qlist.append(2)
if "node3" in queue:
    qlist.append(3)    
if "node4" in queue:
    qlist.append(4)
if "node5" in queue:
    qlist.append(5)
if "node6" in queue:
    qlist.append(6)
if "node" in queue:
    queue_command=queue_command+ ' "' 
for q in qlist:
    queue_command=queue_command+"all.q@cms-0-"+str(q)+".local,"
if len(qlist) > 0    :
    queueok=True
    queue_command=queue_command[:-1]
    queue_command=queue_command+'" '

if queue == "None":
    queue=""
if not  queue:
    queue=""

memoryusage_v=0
memoryusage_p=0
list_of_extra_lib=[]
libname=''
for lib in tmplist_of_extra_lib:
    if '"' in libname:
        libname=""
    if ',' in libname:
        libname=""
        
    libname+=lib
    if ".so" in libname:
        list_of_extra_lib.append(libname)
        libname=""
    
if libname != "NULL":
    if len(list_of_extra_lib) ==0:
        print "Name of library has to contain .so."
   
for lib in list_of_extra_lib:
    print "Adding " + lib + " to list of Libraies to Load"


if DEBUG == "True":
    print "In debug mode"

print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print "Running : " + cycle

if useskim == "SKTree_NoSkim":
    useskim ="NoCut"
elif useskim == "SKTree_LeptonSkim":
    useskim ="Lepton"
elif useskim == "SKTree_DiLepSkim":
    useskim="DiLep"
elif useskim == "SKTree_HNDiLepSkim":
    useskim="HNDiLep"
elif useskim == "SKTree_HNFakeSkim":
    useskim="HNFake"
elif useskim == "SKTree_TriLepSkim":
    useskim="TriLep"
        

if useskinput == "True": 
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    print "Using SKTrees as input."
elif useskinput == "true":
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    print "Using SKTrees as input."
else:
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    print "Using CATntuples as input"    
########  Sample specific configuration ###############
## set the job conguration set for a specific sample###
#######################################################
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
print ""

####################
####
####################



###########################################################################################
###########################################################################################
### DEFAULT  settings for runnning in batch mode. True for cmscluster.snu.ac.kr
running_batch=False

if number_of_cores == 1:
    print "running_batch = False  (default if njobs = 1)"
else:
    running_batch=True
    
if str(usebatch) == "NULL":
    if str(running_batch) == "True":
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print "Running batch job:"
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    else:
        print "Running standard root interactive job:"
else:
    if str(usebatch) == "False":
        running_batch=False
    elif  str(usebatch) == "false":
        running_batch=False
        


if not cycle == "SKTreeMaker":
    if not cycle == "SKTreeMakerNoCut":
        if not cycle == "SKTreeMakerDiLep":
            if not cycle == "SKTreeMakerFake":
                if not cycle == "SKTreeMakerHNDiLep":
                    if not useskinput == "True":
                        if not useskinput == "true":
                            print "You are running on FlatCATntuples. This will be more cpu extensive. This is only advisable if you are testing some new branches NOT in SKTrees."
                        

output_mounted="/data2"
workoutput_mounted="/data2"
merge_mounted="/data8/DATA"
if len(sample)>1:
    if  sample == "H_v2" or sample == "H_v3":
        workoutput_mounted="/data7/DATA"

else:
    workoutput_mounted="/data7/DATA"


if "cmscluster.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    output_mounted="/data5"

                    
##########################################################
### Make tmp directory for job
############################################################

tmpwork = workoutput_mounted+"/CAT_SKTreeOutput/"+ getpass.getuser() + "/"
mergetmpwork = merge_mounted+"/CAT_SKTreeOutput/"+ getpass.getuser() + "/"
if not (os.path.exists(tmpwork)):
    os.system("mkdir " + tmpwork)

if not (os.path.exists(mergetmpwork)):
    os.system("mkdir " + mergetmpwork)



timestamp_dir=tmpwork + "/" + cycle + "_joboutput_" +now() +"_" +os.getenv("HOSTNAME")+"_"+os.getenv("CATVERSION")
#timestamp_dir=tmpwork + "/" + cycle + "_joboutput_" +now() +"_" +os.getenv("HOSTNAME")    
while  (os.path.exists(timestamp_dir)):
    app_dir=1
    timestamp_dir = timestamp_dir + "_v" + str(app_dir)
    app_dir = app_dir+1
if not (os.path.exists(timestamp_dir)):
    os.system("mkdir " + timestamp_dir)


if not os.path.exists(timestamp_dir+"/job_output/"):
    os.system("mkdir " + timestamp_dir+"/job_output/")

local_sub_dir=  timestamp_dir + "/job_output/"  + sample + '_' + new_channel + '_' + now()
  
if not os.path.exists(local_sub_dir):
    os.system("mkdir " + local_sub_dir)
        
##################################################################################################################
#### HARD CODE THE MAXIMUM number of subjobs
##################################################################################################################

large_sample=False
if "TT" in sample:
    large_sample=True
if "DY50"  in sample:
    large_sample=True


ncore_def=number_of_cores
import platform
BusyMachine=False
username = str(os.getenv("USER"))

nj_def=1000    
debug_sktree = False
if number_of_cores > 0:
    print "Setting ncores to max if sktreemaker"
    if cycle == "SKTreeMaker":
        number_of_cores=nj_def
    if cycle == "SKTreeMakerNoCut":
        number_of_cores=nj_def
    if cycle == "SKTreeMakerDiLep":
        number_of_cores=nj_def
    if cycle == "SKTreeMakerHNDiLep":
        number_of_cores=nj_def
    if cycle == "SKTreeMakerFakeHN":
        number_of_cores=nj_def

    if cycle == "SKTreeMakerTriLep":
        number_of_cores=nj_def
    print "number_of_cores  = " + str(number_of_cores)
else:
    debug_sktree=True
    number_of_cores=1


##################################################################################################################            
##### FINISHED CONFIGURATION
##################################################################################################################
singlejob = number_of_cores==1            

if "SKTreeMaker" in cycle:
    if debug_sktree:
        singlejob=True
        running_batch=False
    else:
        singlejob=False
        running_batch=True

#### determine if input is data/mc
mc = len(sample)>1
if sample == "H_v2" or sample == "H_v3":
    mc= False

if mc:
    datatype="mc"

else:
    datatype="data"

if datatype == "mc":
    timeWait=10


##################################################################################################################
### Make a list of input samples: at the moment this is useless. Will add code to include * options in input
##################################################################################################################
list = []
import re
if ("*" in sample) and mc:
    print "ADD code"
else:
    list.append(sample)
#list has only size ==1 currently

##################################################################################################################
##### Specify if the job is running on SKTrees or CATNtuples
##################################################################################################################
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
                    if useskim == "HNDiLep":
                        new_channel="SK" + new_channel + "_hndilep"
                    else:
                        if useskim == "HNFake":
                            new_channel="SK" + new_channel + "_hnfake"
                            
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
                    if useskim == "HNDiLep":
                        sample="SK" + sample + "_hndilep"
                    else:
                        if useskim == "HNFake":
                            sample="SK" + sample + "_hnfake"
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
                    if useskim == "HNDiLep":
                        new_channel="SK" + new_channel + "_hndilep"
                    else:
                        if useskim == "HNFake":
                            new_channel="SK" + new_channel + "_hnfake"
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
                    if useskim == "HNDiLep":
                        sample="SK" + sample + "_hndilep"
                    else:
                        if useskim == "HNFake":
                            sample="SK" + sample + "_hnfake"
                        else:
                            if useskim == "TriLep":
                                sample="SK" + sample + "_trilep"
                            
                
print "Input sample = " + sample
if not mc:
    print "Input channel = " + new_channel

##############################################################################################
#### Check if sktrees are located on current machines  (not used when running on cmsX at snu)                        
#############################################################################################

isfile = os.path.isfile
join = os.path.join
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

if "HN" in  sample or "CHT" in sample:
    datasetfile="datasets_snu_sig_CAT_mc_"
else:
    datasetfile="datasets_snu_nonsig_CAT_mc_"

while inDS == "":
    if platform.system() == "Linux":
        version="_CAT"
        sample_catversion = catversions[iversion]


        if catversion != "":
            sample_catversion = catversion
            output_catversion = catversion

        print "Using CAT " +sample_catversion + " ntuples"
        if mc:
            filename = os.getenv("LQANALYZER_RUN_PATH") + '/txt/'+datasetfile +sample_catversion +  '.txt'

        else:
            filename = os.getenv("LQANALYZER_RUN_PATH") + '/txt/datasets_snu_CAT_data_'  +sample_catversion +'.txt'
    else:
        filename = os.getenv("LQANALYZER_RUN_PATH") + 'txt/datasets_mac.txt'
             

    print "Using " + filename    
    if not mc:
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
        update = raw_input("This is likely because you have not changed the name of the input file. Since CATVERSION v7-6-3 these were changed. Is using anolder version of catuples ok for what you are doing? If you wish to use an older sample type Y. If not change input. run 'sktree -l' for options" )
        if not  update == "Y":
            sys.exit()
        if iversion == len(catversions):
            print "LQAnalyzer :: ERROR :: Input dataset is not available: Exiting"
            sys.exit()

if "tamsa2.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    inDS=inDS.replace("/data2/","/data4")

InputDir = inDS    


##################################################################################################################
print "Input directory= " + inDS    ## now have defined what dur contains input files
##################################################################################################################                    

############################################################
############################################################
###### RUN JOB
############################################################
############################################################

os.system("ls " + InputDir + "/*.root > " + local_sub_dir + "/inputlist.txt")

############################################################
## Get number of files in Input directory
############################################################
isfile = os.path.isfile
join = os.path.join
number_of_files = sum(1 for item in os.listdir(InputDir) if isfile(join(InputDir, item)))

print "number_of_files = " + str(number_of_files)
if number_of_files == 1 and not setnumber_of_cores:
    singlejob=False
    running_batch=True
if number_of_files == 1 and submit_allfiles:
    singlejob=False
    running_batch=True



if DEBUG == "True":
    print "Job has " + str(number_of_files) + " files to process:"




cluster = False
if "tamsa2.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    cluster=True
    if not running_batch:
        number_of_cores = 1
        print "Can only run 1 job when running on " +  str(os.getenv("HOSTNAME"))

if not "SKTreeMaker" in cycle: 
    if ncore_def == 1:
        number_of_cores = 1
############################################################
### Correct user if ncores is > nfiles
############################################################
if number_of_cores > number_of_files:
    number_of_cores = number_of_files
print "Splitting job into " + str(number_of_cores) + " subjobs"






############################################################
### set number of files per job
############################################################
nfilesperjobs= 0
for i in range(1,number_of_files):
    if not i%number_of_cores:
        nfilesperjobs+=1

if number_of_cores == 1:
    nfilesperjobs = number_of_files

if nfilesperjobs == 0:
    nfilesperjobs=1
    
files_torun = (nfilesperjobs*number_of_cores)
remainder = number_of_files - (nfilesperjobs*number_of_cores)

if DEBUG == "True":
    print "Each job will process  " + str(nfilesperjobs) + "/" + str(nfilesperjobs+1) + " files"


###################################################
## counters
###################################################
nfiles=0
count=1
total_nsamples=0
filesprocessed=0
nfiles_file=0
n_remainder_files=0
check_array = []

###################################################
# Setup work area on var tmp
###################################################

    
workspace = workoutput_mounted+"/CAT_SKTreeOutput/"+ getpass.getuser() +"/"

if not (os.path.exists(workspace)):
        os.system("mkdir " + workspace)
out_end=sample


output=workspace + sample + "_" + now() + "_" + os.getenv("HOSTNAME")  + "/"
output_tag=2
while os.path.exists(output):
    output=workspace + sample + "_" + now() + "_" + os.getenv("HOSTNAME")  + "v" +str(output_tag) +"/"
    output_tag=output_tag+1

mergeoutputdir = mergetmpwork + sample + "_" + now() + "_" + os.getenv("HOSTNAME")  + "/"
merge_tags=2
while os.path.exists(mergeoutputdir):
    mergeoutputdir=mergetmpwork + sample + "_" + now() + "_" + os.getenv("HOSTNAME")  +"v" +str(merge_tags)+"/"
    merge_tags=merge_tags+1


if not mc:
    output=workspace + new_channel+ "_"+ sample + "_" + now() + "_" + os.getenv("HOSTNAME") + "/" 
    mergeoutputdir = mergetmpwork + new_channel+ "_"+ sample + "_" + now() + "_" + os.getenv("HOSTNAME") + "/"
    merge_tags=2
    while os.path.exists(mergeoutputdir):
        mergeoutputdir = mergetmpwork + new_channel+ "_"+ sample + "_" + now() + "_" + os.getenv("HOSTNAME") + "v" +str(merge_tags)+"/"


outputdir= output+ "output/"
outputdir_tmp= output+ "output_tmp/"
if not (os.path.exists(output)):
    os.system("mkdir " + output)
    print "Making tmp working directory to run Job  : " + output
if not (os.path.exists(mergeoutputdir)):
    os.system("mkdir " + mergeoutputdir)
    print "Making tmp merge dir " + mergeoutputdir


if(os.path.exists(outputdir)):
    number_of_outputfiles = sum(1 for item in os.listdir(outputdir) if isfile(join(outputdir, item)))
    if  not number_of_outputfiles ==0:
       os.system("rm " + outputdir + "/*.root")
       print "Emptying output directory as this should be empty for new job"
              
if not (os.path.exists(outputdir)):
    os.system("mkdir " + outputdir)
    os.system("mkdir " + outputdir_tmp)

###################################################
## Make subjob directories
###################################################
printedworkdir =  output + "Job_[" + str(1) + "-" + str(number_of_cores) + "]/"
for i in range(1,number_of_cores+1):
    workdir =  output + "Job_" + str(i) + "/"    
    if not (os.path.exists(workdir)):
            os.system("mkdir " + workdir)
            if i==1:
                if DEBUG == "True":
                    print "making sub work directories " + printedworkdir


####################################################
## Creat separate input lists/macros for each subjob
####################################################

### read inputlist.txt which contains all input files
fr = open(local_sub_dir + '/inputlist.txt', 'r')

outsamplename = sample
if runnp == "True":
    outsamplename = "nonprompt_" + outsamplename
    print "sample --> " + outsamplename
if runcf == "True":
    outsamplename = "chargeflip_" + outsamplename
    print "sample --> " + outsamplename
if not mc:
    outsamplename = outsamplename +  "_" + new_channel
    if useCATv742ntuples == "True":
        outsamplename = outsamplename + "_cat_" + output_catversion

else:
    if useCATv742ntuples == "True":
                outsamplename = outsamplename + "_cat_"+ output_catversion

if  "SKTreeMaker" in cycle:            
    outsamplename = outsamplename +  os.getenv("CATTAG")

if tmp_filename != "":
    outsamplename = outsamplename +"_"+tmp_filename

### specify the location of the macro for the subjob     
printedrunscript = output+ "Job_[1-" + str(number_of_cores)  + "]/runJob_[1-" + str(number_of_cores)  + "].C"

for line in fr:
    if nfiles < files_torun:
        if nfiles == 0 :
            runscript = output+ "Job_" + str(count) + "/runJob_" + str(count) + ".C"
            filelist = output+ "Job_" + str(count) + "/" + sample + "_%s" % (count) + ".txt"
            fwrite = open(filelist, 'w')
            configfile=open(runscript,'w')
            configfile.write(makeConfigFile(loglevel, outsamplename, filelist, tree, cycle, count, outputdir_tmp, outputdir, number_of_events_per_job, logstep, skipev, datatype, original_channel, data_lumi, totalev, xsec, tar_lumi, eff_lumi, useskinput, runevent, list_of_extra_lib, runnp,runcf,runtau, skflag,tagger,useskim)) #job, input, sample, ver, output
            configfile.close()
            if DEBUG == "True":
                print "Making file : " + printedrunscript
            fwrite.write(line)
            filesprocessed+=1
            nfiles_file+=1            
            nfiles+=1
            if files_torun == 1:
                fwrite.close()
            continue

        #End of file
        if not nfiles % nfilesperjobs:
            if not nfiles == number_of_files :
                # set counters
                nfiles_file=0
                count+=1        
                # close files
                fwrite.close()
                ### Make next set of scripts
                runscript = output+ "Job_" + str(count) + "/runJob_" + str(count) + ".C"
                filelist = output+ "Job_" + str(count) + "/" + sample + "_%s" % (count) + ".txt"
                fwrite = open(filelist, 'w')
                configfile=open(runscript,'w')
                configfile.write(makeConfigFile(loglevel,outsamplename, filelist, tree, cycle, count, outputdir_tmp,outputdir, number_of_events_per_job, logstep, skipev, datatype , original_channel, data_lumi, totalev, xsec, tar_lumi, eff_lumi, useskinput, runevent,list_of_extra_lib, runnp, runcf,runtau, skflag,tagger,useskim))
                configfile.close()
                fwrite.write(line)
                filesprocessed+=1
                nfiles_file+=1
            else:
                fwrite.write(line)
                filesprocessed+=1
                nfiles_file+=1
                if DEBUG == "True":
                    print "File " + filelist + " contains " + str(nfiles_file) + " files"            
                
        else:
            fwrite.write(line)
            filesprocessed+=1
            nfiles_file+=1

        if nfiles == number_of_files-1 :
            fwrite.close()

    else:
        n_remainder_files+=1
        filelist = output+ "Job_" + str(n_remainder_files) + "/" + sample + "_%s" % (n_remainder_files) + ".txt"
        #runscript = output+ "Job_" + str(count) + "/runJob_" + str(count) + ".C"
        fwrite = open(filelist, 'a')
        fwrite.write(line)
        #configfile=open(runscript,'w')
        #configfile.write(makeConfigFile(loglevel,sample, filelist, tree, cycle, count, outputdir_tmp,outputdir, number_of_events_per_job, logstep, skipev, datatype , original_channel, data_lumi, totalev, xsec, tar_lumi, eff_lumi, useskinput, runevent,list_of_extra_lib, runnp, runcf, skflag))
        #configfile.close()
        filesprocessed+=1
        fwrite.close()        
    nfiles+=1        
fr.close()

#################################################################### 
### Check Final input files have no duplicates
#################################################################### 
no_duplicate=False
for check in range(1, number_of_cores+1):
    filelist = output+ "Job_" + str(check) + "/" + sample + "_%s" % (check) + ".txt"
    fcheck = open(filelist, 'r')
    nsamples=0
    for line in fcheck:
        nsamples+=1
        total_nsamples+=1
        no_duplicate= True
        for s in check_array:
            if s == line :
                print "DUPLICATE file : " + s
                no_duplicate=False
                sys.exit()
        check_array.append(line)
    if DEBUG == "True":
        print  "File " + filelist + " contains " + str(nsamples) + " files"
    fcheck.close()
if DEBUG == "True":
    print "Total Number of input files = " + str(total_nsamples)     

if no_duplicate:
    #print "Checking for duplicates: "
    #print "Checking for duplicates:...... "
    if DEBUG == "True":
        print "Checking for duplicates: NONE found"
else:
     #print "Checking for duplicates: "
     #print "Checking for duplicates:...... "
     print "Checking for duplicates: Duplicate files found. Check script "
            
if DEBUG == "True":
    print "Total number of files processed = " + str(filesprocessed)

###################################################
### Run each .C file in background
###################################################
import thread,time
start_time = time.time()
start_running_time=0.
job_time=0
last_job_time=0
wait_sub = 1
if number_of_cores < 10:
    wait_sub = 5


if DEBUG == "True":
    print "Running CATAnalyzer jobs for: " + getpass.getuser()

array_batchjobs = []

for i in range(1,number_of_cores+1):

    batchscript =  output+ "Job_" + str(i) + "/runJob_" + str(i) + ".sh"
    batchfile=open(batchscript,'w')
    batchfile.write(make_batch_script(output+ "Job_" + str(i) , outsamplename+ "_Job_" + str(i),str(os.getenv("LQANALYZER_DIR")),"runJob_" + str(i) + ".C",cluster))
    batchfile.close()

    
    script = output+ "Job_" + str(i) + "/runJob_" + str(i) + ".C"
    log = output+ "Job_" + str(i) + "/runJob_" + str(i) +".log"
    logbatch="Job_" + str(i) + "/"+outsamplename+ "_Job_" + str(i)+".o[batchID]"  

    runcommand = "nohup root.exe -l -q -b " +  script + "&>" + log + "&"
    if running_batch:
        if queue == "":
            runcommand = "qsub -V " + batchscript   + "&>" + log 
        else:
            runcommand = "qsub -V  " + queue_command +" " + batchscript   + "&>" + log
    jobID=0
    first_jobid=0


    if singlejob:
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print "Running single job " + script 
        runcommand = "root.exe -l -q -b " +  script 
        os.system(runcommand)
    else:
        if i== number_of_cores:
            if running_batch:
                print "Running " + script + " . Log file --->  " + logbatch
            else:
                print "Running " + script + " . Log file --->  " + log

        elif i==1:
            if running_batch:
                print "Running " + script + " . Log file --->  " + logbatch
            else:
                print "Running " + script + " . Log file --->  " + log

        elif i==2:
            print "......"
        os.system(runcommand)
        if running_batch:
            for line in open(log, 'r'):
                entries = line.split()
                if len(entries) > 0:
                    if jobID ==0:
                        first_jobid= entries[2]
                    jobID=entries[2]
                    array_batchjobs.append(jobID)
                    k_batchscript =  output+ "JobKill.sh"
                    if i==1:
                        k_batchfile=open(k_batchscript,'w')
                        k_batchfile.write("qdel  " +str(jobID) +';')
                    else:
                        k_batchfile=open(k_batchscript,'a')
                        k_batchfile.write("qdel  " +str(jobID)+';')
                    k_batchfile.close()

if running_batch: 
    
    if not os.path.exists(workoutput_mounted+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/"):
        os.system("mkdir " + workoutput_mounted+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/")
    

    print "@@@@@@@@@@@@@@@@@@@@@@@@@"
    path_jobids = workoutput_mounted+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + original_sample + "jobid.txt"

    file_jobids = open(path_jobids,"w")
    for ijob in array_batchjobs:
        file_jobids.write(ijob + "\n")
        print "Job ["+str(ijob)+"] added to list......."
    file_jobids.close()    
    print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    print "In case user wants to kill job do : source " + output+ "JobKill.sh"
    print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"




##########################################################
## wait and do merging (also remove old log file/rootfiles
##########################################################
check_log= os.getenv("LQANALYZER_LOG_PATH") + "/" + outsamplename + "/"
if number_of_cores > 1:
    if (os.path.exists(check_log)):
        if sum(1 for item in os.listdir(check_log) if isfile(join(check_log, item))) > 0:
            os.system("rm " + os.getenv("LQANALYZER_LOG_PATH") + "/" + outsamplename + "/*.o*")

if DEBUG == "True":
    print "Waiting for all jobs to finish before Merging."

ncomplete_files=0
JobSuccess=False
JobOutput=True
CompletedJobs=[]
doMerge=False

print "Checking Job status:"

if singlejob:    
    check_outfile = outputdir + outsamplename +  "_1.root"
    if (os.path.exists(check_outfile)):
        JobSuccess=True;
else:        
    sys.stdout.write('\r'+ '0% of events processed...' )
    sys.stdout.flush()

low_cpu=0
ncycle=0
file_iterator=0
files_done= []
clear_line='                                                                          '
failed_macro=""
failed_log=""

JobCrash=False
checkJob=True

check_cluster="_"
n_cms1=0
n_cms2=0
n_cms3=0
n_cms4=0
n_cms5=0
n_cms6=0

it_counter=0

while not JobSuccess:
    it_counter+=1
    sub_counter=0
    if running_batch == False:
        os.system("ps ux &> " + local_sub_dir + "/log")
    else: 
        os.system("qstat -u '*' > " +  local_sub_dir + "/log")
        
    filename = local_sub_dir +'/log'
    running = False
    
    if running_batch == False:

        for line in open(filename, 'r'):
            if "root.exe" in line:
                running = True            
                splitline  = line.split()
                if splitline[2] < 0.1:
                    low_cpu+=1
        if low_cpu > 3:
            running = False

    else:
        for line in open(filename, 'r'):
            for job_id in array_batchjobs:
                if job_id in line:
                    running = True
                    entries = line.split()
                    ### Check the job is not held/suspended
                    if len(entries) ==  9:  
                        idcheck="_"+job_id+"_"
                        if not idcheck in check_cluster:
                            if "cms-0-1" in entries[7]:
                                n_cms1+=1
                            if "cms-0-2" in entries[7]:
                                n_cms2+=1
                            if "cms-0-3" in entries[7]:
                                n_cms3+=1
                            if "cms-0-4" in entries[7]:
                                n_cms4+=1
                            if "cms-0-5" in entries[7]:
                                n_cms5+=1
                            if "cms-0-6" in entries[7]:
                                n_cms6+=1

                        check_cluster=check_cluster+job_id+"_"

                        if entries[4] == "h":
                            print "Job " + str(job_id) + " is in held state: killing all jobs"
                            os.system("source " + output+ "JobKill.sh")
                            running = False
                        if entries[4] == "s":    
                            print "Job " + str(job_id) + " is in suspended state: killing all jobs"
                            os.system("source " + output+ "JobKill.sh")
                            running = False

    if not running:
        check_outfile = outputdir + outsamplename +  "_1.root"
        if not (os.path.exists(check_outfile)):
            JobSuccess=True
            JobOutput=False
            
    os.system("rm  " + local_sub_dir + "/log")
    
    for i in range(1,number_of_cores+1):
        skipcheck=False
        for check in CompletedJobs:
            if i== check: skipcheck=True
        if not skipcheck:
            if checkJob:
                if sum(1 for item in os.listdir(outputdir) if isfile(join(outputdir, item))) > 0:
                    if not start_running_time == 0.:
                        job_time = time.time()
                        checkJob=False
            check_outfile = outputdir + outsamplename +  "_" +  str(i) + ".root"   
            if (os.path.exists(check_outfile)):
                CompletedJobs.append(i)
                ncomplete_files+=1
                files_done.append("Job [" + str(i) + "] completed. Output ="  + check_outfile)
            
    if ncomplete_files== number_of_cores :
        last_job_time = time.time()
        sys.stdout.write('\r' + clear_line)
        sys.stdout.flush()
        sys.stdout.write('\r'+ '100% of events processed. \n' )
        sys.stdout.flush()
        print "Job finished"
        doMerge=True
        if ncycle == 0:
            print "Job ran in less than 10 seconds. Assumed bug:"
            if number_of_cores == 1:
                JobOutput=True
            else:
                JobOutput=False
        JobSuccess=True
                                    
    else:
        nevents_total=0.
        nevent_processed=0.

        sub_counter+=1


        if running_batch:
            ### print jobs running/in queue .... once all running print % completeion
            for i in range(1,number_of_cores+1):
                check_outfile = output + "/Job" +  "_" +  str(i) + "/" + outsamplename + "_Job_"+ str(i) +".o"+array_batchjobs[i-1]
                    
                if start_running_time == 0.:    
                    if i == 1:
                        start_running_time = time.time()    
                while not os.path.exists(check_outfile):
                    if i == 1:
                        start_running_time = time.time()

                    sys.stdout.write('\r' + clear_line)
                    sys.stdout.flush()
                    sys.stdout.write('\r'+ 'Current jobs running : [' + str(i-1) + '/' + str(number_of_cores) + ']... '+ str(number_of_cores-i+1) + ' in queue' )
                    sys.stdout.flush()
                    if checkJob:
                        if sum(1 for item in os.listdir(outputdir) if isfile(join(outputdir, item))) > 0:
                            job_time = time.time()
                            checkJob=False
                    time.sleep(1.)
            if ncycle == 0:
                sys.stdout.write('\r' + clear_line)
                sys.stdout.flush()
                sys.stdout.write('\r'+ 'Current jobs running : [' + str(number_of_cores) + '/' + str(number_of_cores) + ']... ')
                sys.stdout.flush()
                time.sleep(1.)
            #### check job is running. Halted or suspended and if not running is output file missing?
            os.system("qstat -u " + getpass.getuser()  + " > " +  local_sub_dir + "/log")
            filename = local_sub_dir +'/log'


            sub_counter+=1

            for i in range(1,number_of_cores+1):
                if not JobOutput:
                    break
                job_id_c=array_batchjobs[i-1]
                job_finished=True
                for line in open(filename, 'r'):
                    if job_id_c in line:
                        job_finished=False
                        entries = line.split()
                        if len(entries) ==  9:
                            if entries[4] == "h":
                                print "Job " + str(job_id_c) + " is in held state: killing all jobs"
                                os.system("source " + output+ "JobKill.sh")
                                number_of_cores=0
                                JobSuccess=True
                                JobOutput=False
                            if entries[4] == "s":
                                print "Job " + str(job_id_c) + " is in suspended state: killing all jobs"
                                os.system("source " + output+ "JobKill.sh")
                                number_of_cores=0
                                JobSuccess=True
                                JobOutput=False
                                
                if job_finished:
                    ### job id not in qstat output. Check if rootfile is missing. If so kill job
                    check_outfile = outputdir + outsamplename +  "_" + str(i)+".root"
                    if not (os.path.exists(check_outfile)):
                        failed_macro= output+ "Job_" + str(i) + "/runJob_" + str(i) + ".C"
                        failed_log= outsamplename+ "_Job_" + str(i)+".o"+str(array_batchjobs[i-1])
                        JobSuccess=True
                        JobOutput=False
                        print "Job " + str(job_id_c) + " is not running or in queue. Output " + str(check_outfile)+ " is missing."
                        print "Most likely a crash occurred.  So killing all jobs." 
                        os.system("source " + output+ "JobKill.sh")

                        check_error_outfile = output + "/Job" +  "_" +  str(i) + "/"+ outsamplename+ "_Job_"+ str(i) +".e"+array_batchjobs[i-1] 
                        print "Error file for job ["+str(job_id_c)+"] shows:"
                        for line in open(check_error_outfile, 'r'):
                            print line
                        print "Run in non-batch mode by setting njobs = 1 in submit.sh file and retry" 
                        JobCrash=True
                        
            os.system("rm  " + local_sub_dir + "/log")
            
        ##### Run over log to get % completion                    
        for i in range(1,number_of_cores+1):
            check_outfile = output + "/Job" +  "_" +  str(i) + "/runJob_"+ str(i) +".log"
            if running_batch == True:
                check_outfile = output + "/Job" +  "_" +  str(i) + "/" + outsamplename + "_Job_"+ str(i) +".o"+array_batchjobs[i-1]

            os.system('tail -100 ' + check_outfile + ' > ' + local_sub_dir + '/outlog.txt')
            nevent_processed_i=0.
            nevents_total_i=0.
            for line in open(local_sub_dir + '/outlog.txt', 'r'):
                if "Processing entry" in line:
                    if "LQCycleController" not in line:
                        entries = line.split()
                        if len(entries)> 6:                        
                            num = entries[7]
                            s = num.replace("/", " ")
                            event_split = s.split()
                            if len(event_split) < 2:
                                print "Error [2002] " + event_split 
                                print line
                                os.system("cp " + local_sub_dir + '/outlog.txt ~/error_log_'+str(array_batchjobs[i-1]))
                            nevent_processed_i = float(event_split[0])
                            nevents_total_i= float(event_split[1])
            nevent_processed+=nevent_processed_i                
            nevents_total+=nevents_total_i

        ## calculate time left to run jobs    
        end_time = time.time()
        remaining = (end_time - start_time)
        if nevent_processed !=0:
            remaining*= (1./nevent_processed)
            remaining*= ( nevents_total  - nevent_processed)
            b = str(round(float(100.*nevent_processed/nevents_total), 2)) + "% of events processed. Estimated cpu time remaining =  " + str(round(remaining,2)) + "[s]" 
            sys.stdout.write('\r' + clear_line)
            sys.stdout.flush()
            sys.stdout.write('\r'+b)            
            sys.stdout.flush()
        else:
            mess="Job Initialising 0% of events processed."
            if JobOutput:
                sys.stdout.write('\r' + clear_line)
                sys.stdout.flush()
                sys.stdout.write('\r'+mess)
                sys.stdout.flush()
                time.sleep(1.)
        if ncomplete_files > file_iterator:
            #print str(ncomplete_files) + "/" + str(number_of_cores) + " jobs completed.  " #Wait " + str(timeWait) + " second..."
            #print ""
            file_iterator=ncomplete_files
        time.sleep(1.)
        ncycle+=1
        sub_counter+=1


if not JobOutput:

    if not running_batch:
        failed_macro= output+ "Job_1/runJob_1.C"
        failed_log= outsamplename+ "_Job_" + str(0)+".log"


   
    print ""
    print "Job Failed...."
    if not os.path.exists(os.getenv("LQANALYZER_LOG_PATH")):
        os.system("mkdir " + os.getenv("LQANALYZER_LOG_PATH"))
    if not os.path.exists(os.getenv("LQANALYZER_LOG_PATH")+ "/" + outsamplename):
        os.system("mkdir " + os.getenv("LQANALYZER_LOG_PATH")+ "/" + outsamplename)
        
    if not number_of_cores == 1:
        os.system("mv "+ output + "/*/*.o* " + os.getenv("LQANALYZER_LOG_PATH") + "/" + outsamplename)    
        os.system("mv "+ output + "/*/*.e* " + os.getenv("LQANALYZER_LOG_PATH") + "/" + outsamplename)    
    print "###########################################################################################################"
    print "Check crash by running root -q -b " + failed_macro 
    print "after typing: export cat_version=$CATVERSION"
    print 'export CATAnalyzerPeriod="B"'
    print "Logfile of failed job is can be found at " + os.getenv("LQANALYZER_LOG_PATH") + "/" + outsamplename   + failed_log 
    print "###########################################################################################################"
    JobCrash=True
    #os.system("rm -r " + output)    
    #os.system("rm -r " + local_sub_dir)    
    #os.system("rm -r " + timestamp_dir)
    
else:

    for line in files_done:
        if DEBUG == "True":
            print line

    SKTreeOutput_pre = workoutput_mounted+"/CatNtuples/" + sample_catversion
    if not os.path.exists(SKTreeOutput_pre):
        os.system("mkdir " + SKTreeOutput_pre)

    SKTreeOutput_pre2 = workoutput_mounted+"/CatNtuples/" + sample_catversion + "/SKTrees/"
    if not os.path.exists(SKTreeOutput_pre2):
        os.system("mkdir " + SKTreeOutput_pre2)
                    
    SKTreeOutput = workoutput_mounted+"/CatNtuples/" + sample_catversion + "/SKTrees/"        
    
    #do not merge the output when using tree maker code
    if cycle == "SKTreeMaker":
        if not os.path.exists(SKTreeOutput):
            os.system("mkdir  " + SKTreeOutput)
        doMerge=False
        if not mc:
            Finaloutputdir = SKTreeOutput + "Data/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
                
            if original_channel =="DoubleEG":
                Finaloutputdir += "DoubleEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleMuon":
                Finaloutputdir += "DoubleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="MuonEG":
                Finaloutputdir += "MuonEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="SingleMuon":
                Finaloutputdir += "SingleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="SingleElectron":
                Finaloutputdir += "SingleElectron/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="SinglePhoton":
                Finaloutputdir += "SinglePhoton/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            Finaloutputdir += "period" + original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
        else:
            Finaloutputdir = SKTreeOutput + "MC/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
                                                                                            
            Finaloutputdir +=  original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
                os.system("chmod 777 -R " +  Finaloutputdir)
    if cycle == "SKTreeMakerNoCut":
        doMerge=False
        if not os.path.exists(SKTreeOutput):
            os.system("mkdir " + SKTreeOutput)
        if not mc:
            Finaloutputdir = SKTreeOutput + "DataNoCut/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleElectron":
                Finaloutputdir += "DoubleElectron/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleMuon":
                Finaloutputdir += "DoubleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="MuonEG":
                Finaloutputdir += "MuonEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="SinglePhoton":
                Finaloutputdir += "SinglePhoton/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
                     
            Finaloutputdir += "period" + original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
        else:
            Finaloutputdir = SKTreeOutput + "MCNoCut/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            Finaloutputdir += original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
                os.system("chmod 777 -R " +  Finaloutputdir)
    if cycle == "SKTreeMakerDiLep":
        doMerge=False
        if not os.path.exists(SKTreeOutput):
            os.system("mkdir " + SKTreeOutput)
        if not mc:
            Finaloutputdir = SKTreeOutput + "DataDiLep/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleEG":
                Finaloutputdir += "DoubleEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleMuon":
                Finaloutputdir += "DoubleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
                    
            if original_channel =="SingleMuon":
                Finaloutputdir += "SingleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="SingleElectron":
                Finaloutputdir += "SingleElectron/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
                    
            if original_channel =="MuonEG":
                Finaloutputdir += "MuonEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="SinglePhoton":
                Finaloutputdir += "SinglePhoton/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            Finaloutputdir += "period" + original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
        else:
            Finaloutputdir = SKTreeOutput + "MCDiLep/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            Finaloutputdir +=  original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
                os.system("chmod 777 -R " +  Finaloutputdir)



    if cycle == "SKTreeMakerHNDiLep":
        doMerge=False
        if not os.path.exists(SKTreeOutput):
            os.system("mkdir " + SKTreeOutput)
        if not mc:
            Finaloutputdir = SKTreeOutput + "DataHNDiLep/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleEG":
                Finaloutputdir += "DoubleEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleMuon":
                Finaloutputdir += "DoubleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)

            if original_channel =="SingleMuon":
                Finaloutputdir += "SingleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="SingleElectron":
                Finaloutputdir += "SingleElectron/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="MuonEG":
                Finaloutputdir += "MuonEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)

            Finaloutputdir += "period" + original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)

        else:
            Finaloutputdir = SKTreeOutput + "MCHNDiLep/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            Finaloutputdir +=  original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
                os.system("chmod 777 -R " +  Finaloutputdir)
    if cycle == "SKTreeMakerFakeHN":
        doMerge=False
        if not os.path.exists(SKTreeOutput):
            os.system("mkdir " + SKTreeOutput)
        if not mc:
            Finaloutputdir = SKTreeOutput + "DataHNFake/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleEG":
                Finaloutputdir += "DoubleEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleMuon":
                Finaloutputdir += "DoubleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)

            if original_channel =="SingleMuon":
                Finaloutputdir += "SingleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="SingleElectron":
                Finaloutputdir += "SingleElectron/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)

            Finaloutputdir += "period" + original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
        else:
            Finaloutputdir = SKTreeOutput + "MCHNFake/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            Finaloutputdir +=  original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
                os.system("chmod 777 -R " +  Finaloutputdir)
                

    if cycle == "SKTreeMakerTriLep":
        doMerge=False
        if not os.path.exists(SKTreeOutput):
            os.system("mkdir " + SKTreeOutput)
        if not mc:
            Finaloutputdir = SKTreeOutput + "DataTriLep/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleEG":
                Finaloutputdir += "DoubleEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="MuonEG":
                Finaloutputdir += "MuonEG/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            if original_channel =="DoubleMuon":
                Finaloutputdir += "DoubleMuon/"
                if not os.path.exists(Finaloutputdir):
                    os.system("mkdir " + Finaloutputdir)
            Finaloutputdir += "period" + original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
                                 
        else:
            Finaloutputdir = SKTreeOutput + "MCTriLep/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
            Finaloutputdir +=  original_sample + "/"
            if not os.path.exists(Finaloutputdir):
                os.system("mkdir " + Finaloutputdir)
                os.system("chmod 777 -R " +  Finaloutputdir)
                
    if not os.path.exists(Finaloutputdir):
        os.system("mkdir " + Finaloutputdir)


    outfile = cycle + "_" + filechannel + outsamplename + ".root"
    if doMerge:
        print "doing merge"
        if not mc:
            outfile = cycle + "_" + outsamplename + ".root"
        if os.path.exists(Finaloutputdir + outfile):
            os.system("rm  "  +  Finaloutputdir   + outfile)
        os.system("mv " +outputdir + "*.root" + " " + mergeoutputdir)
        os.system("hadd " + mergeoutputdir +  outfile  + " "+ mergeoutputdir + "*.root")
        
        if mc:
            if not os.path.exists( "/data2/CAT_SKTreeOutput/"+os.getenv("USER")+"/Histdir" + tagger ):
                os.system("mkdir " +  "/data2/CAT_SKTreeOutput/"+os.getenv("USER")+"/Histdir" + tagger)
            os.system("source "+os.getenv("LQANALYZER_DIR")+"/scripts/Counter.sh " + mergeoutputdir +  outfile + " > /data2/CAT_SKTreeOutput/"+os.getenv("USER")+"/Histdir" + tagger + "/"+original_sample+"Hist.txt"   )
            os.system("source "+os.getenv("LQANALYZER_DIR")+"/scripts/CutFlow.sh " + mergeoutputdir +  outfile + " > /data2/CAT_SKTreeOutput/"+os.getenv("USER")+"/Histdir" + tagger + "/"+original_sample+"CutFlow.txt"   )  

        if os.getenv("USER") == "jalmond":

            transout=Finaloutputdir.replace("/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/","/afs/cern.ch/work/j/jalmond/CAT/")
            
            catpath=os.getenv("LQANALYZER_DIR")+"/bin/catconfig"
            readcatpath=open(catpath,"r")
            lxmachine=""
            for rline in readcatpath:
                if "localcpu" in rline:
                    srline = rline.split()
                    lxmachine=srline[2]
            readcatpath.close()
            if not "OPT/" in Finaloutputdir:
                os.system("ssh  jalmond@"+lxmachine+".cern.ch  mkdir " + transout )
                if "OPT" in skflag:
                    transout = transout+"/OPT/"

            os.system("ssh  jalmond@"+lxmachine+".cern.ch  mkdir " + transout )
            os.system("scp -r "+mergeoutputdir +  outfile + " jalmond@"+lxmachine+".cern.ch:"+transout)
            

        os.system("mv " + mergeoutputdir +  outfile + " "  + Finaloutputdir)
        os.system("ls -lh " + Finaloutputdir +  outfile + " > " + path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/filesize_" + original_sample+ tagger +".txt")
        f = ROOT.TFile(Finaloutputdir +  outfile)
        t = f.Get("CycleInfo/CycleVirtualMemoryUsage")
        t2 = f.Get("CycleInfo/CyclePhysicalMemoryUsage")
        memoryusage_v=(t.GetBinContent(8)/ number_of_cores)
        memoryusage_p=(t2.GetBinContent(8)/ number_of_cores)
        f.Close()
        print "Merged output :" + Finaloutputdir + outfile
    else:
        print "not doing merge"
        
        if not mc:
            outfile = cycle + "_" + outsamplename + ".root"
        if number_of_cores == 1 and setnumber_of_cores:
            os.system("mv " + outputdir + outsamplename + "_1.root " + Finaloutputdir + outfile )
            
            os.system("ls -lh " + Finaloutputdir +   outfile + " > " + path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/filesize" + tagger+".txt")
            f = ROOT.TFile(Finaloutputdir +  outfile)
            t = f.Get("CycleInfo/CycleVirtualMemoryUsage")
            t2 = f.Get("CycleInfo/CyclePhysicalMemoryUsage")
            memoryusage_v=(t.GetBinContent(8)/ number_of_cores)
            memoryusage_p=(t2.GetBinContent(8)/ number_of_cores)
            f.Close()

        else:
            os.system("rm " + Finaloutputdir + "/*.root")
            os.system("mv " + outputdir + "*.root " + Finaloutputdir )
            os.system("chmod -R 777 " + Finaloutputdir )
            os.system("ls -lh " + Finaloutputdir +  outsamplename + "_1.root > " + path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/filesize_" + original_sample+ tagger+".txt")
            f = ROOT.TFile(Finaloutputdir +  outsamplename + "_1.root")
            t = f.Get("CycleInfo/CycleVirtualMemoryUsage")
            t2 = f.Get("CycleInfo/CyclePhysicalMemoryUsage")
            memoryusage_v=(t.GetBinContent(8)/ number_of_cores)
            memoryusage_p=(t2.GetBinContent(8)/ number_of_cores)
            f.Close()

        if DEBUG == "True":
            print "Non merged output :" +Finaloutputdir
    
    if remove_workspace == "True":
        if not os.path.exists(os.getenv("LQANALYZER_LOG_PATH")):
            os.system("mkdir " + os.getenv("LQANALYZER_LOG_PATH"))
            
        if not os.path.exists(os.getenv("LQANALYZER_LOG_PATH")+ "/" + outsamplename):
            os.system("mkdir " + os.getenv("LQANALYZER_LOG_PATH")+ "/" + outsamplename)
        
        if not number_of_cores == 1:    
            os.system("mv "+ output + "/*/*.o* " + os.getenv("LQANALYZER_LOG_PATH") + "/" + outsamplename)
        os.system("rm -r " + output)
        os.system("rm -r " + local_sub_dir)
        os.system("rm -r " + timestamp_dir)
        print "Log files are sent to  --> "  + os.getenv("LQANALYZER_LOG_PATH")+ "/" + outsamplename
        if doMerge:
            print "All sampless finished: OutFile:"  + cycle + "_" + filechannel + outsamplename + ".root -->" + Finaloutputdir
        else:
            if number_of_cores == 1:
                print "All sampless finished: OutFiles "+ outsamplename + "_1.root -->" + Finaloutputdir + outfile
            else:
                print "All sampless finished: OutFiles "+ outsamplename + "*.root -->" + Finaloutputdir
            
    else:
        print "TMP directory " + output + "is not removed. "
        os.system("rm -r " + local_sub_dir)
        
if os.path.exists(local_sub_dir):
    os.system("rm -r " + local_sub_dir)
    os.system("rm -r " + timestamp_dir)
end_time = time.time()
total_time=end_time- start_time
print "Using " + str(number_of_cores) + " cores: Job time = " + str(total_time) +  " s"
print ""

statfile=path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/" + str(tagger) + "/statlog_"+ original_sample + tagger +".txt"
statfile_time=path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/" + str(tagger) + "/statlog_timetmp_"+original_sample+tagger +".txt"
statfile_time_complete=path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/" + str(tagger) + "/statlog_time_"+original_sample+tagger +".txt"

if not os.path.exists(path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() ):
    os.system("mkdir " + path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser())
if not os.path.exists(path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/" + str(tagger)):
    print "mkdir " + path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser()+ "/" + str(tagger)
    os.system("mkdir " + path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser()+ "/" + str(tagger))

statwrite = open(statfile, 'r')
statwrite_time = open(statfile_time, 'w')
for line in statwrite:
    statwrite_time.write(line)
statwrite.close()

statwrite_time.write("time " + str(total_time) + " \n")
statwrite_time.write("job_time  " + str(job_time-start_running_time)  + " \n")
statwrite_time.write("last_job_time  " + str(last_job_time-start_running_time)  + " \n")
print "time " + str(total_time)
print "job_time  " + str(job_time)
print "start_running_time " + str(start_running_time)
print "last_job_time " + str(last_job_time)

pathfilesize=path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/filesize_"+original_sample + tagger +".txt"

if os.path.exists(pathfilesize):
    readfilesize = open(pathfilesize, "r")
    for line in readfilesize:
        splitline = line.split()
        if len(splitline) == 9:
            statwrite_time.write("outputfile_size " + str(splitline[4])  + " \n")
    readfilesize.close()
    os.system("rm " + pathfilesize)
else:
    statwrite_time.write("outputfile_size 0.0 \n")

statwrite_time.write("Njobs " +  str(number_of_cores) + " \n") 
statwrite_time.write("NFiles " + str(number_of_files)  + " \n")

statwrite_time.write("ClusterIDs " + str(n_cms1) + ":" + str(n_cms2) + ":"+ str(n_cms3) + ":"+ str(n_cms4) + ":"+ str(n_cms5) + ":"+ str(n_cms6) + " \n")
statwrite_time.write("memoryusage_v " + str(memoryusage_v/1000)  + "MB \n") 
statwrite_time.write("memoryusage_p " + str(memoryusage_p/1000)  + "MB \n") 
if JobCrash:
    statwrite_time.write("Success= False \n")
    if not os.path.exists(workoutput_mounted+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/"):
        os.system("mkdir " + workoutput_mounted+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/")
    os.system("mkdir " + workoutput_mounted+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + original_sample+"_crash")
    crash_log= workoutput_mounted+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + original_sample+"_crash/crashlog.txt"
    writecrashlog = open(crash_log,"w")
    writecrashlog.write("###########################################################################################################")
    writecrashlog.write("Check crash by running root -q -b " + failed_macro)
    writecrashlog.write("Logfile of failed job is can be found at " + os.getenv("LQANALYZER_LOG_PATH") + "/" + outsamplename   + failed_log)
    writecrashlog.write("###########################################################################################################")
    writecrashlog.close()
else:
    statwrite_time.write("Success= True \n")

statwrite_time.close()
print "mv " + statfile_time + " " + statfile_time_complete

os.system("mv " + statfile_time + " " + statfile_time_complete)

GeneralStatFile = os.getenv("LQANALYZER_DIR")+ "/python/StatFile.py"

number_of_jobs_for_statfile=number_of_cores

if number_of_events_per_job > 0:
    number_of_jobs_for_statfile = 1
if skipev > 0:
    number_of_jobs_for_statfile = 1
if JobCrash:
    number_of_jobs_for_statfile = 1


print "python " + GeneralStatFile + " -x " + tagger + " -s " + original_sample + " -n " + str(number_of_jobs_for_statfile)

os.system("python " + GeneralStatFile + " -x " + tagger + " -s " + original_sample + " -n " + str(number_of_jobs_for_statfile))

#set_logfile=path_jobpre +"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/" + str(tagger)+ "/statlog_"+original_sample+ tagger + ".txt"
#if os.path.exists(set_logfile):
#    os.system("rm " + set_logfile)
