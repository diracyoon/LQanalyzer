import os

def MakeDirectory(dirpath):
    predir=dirpath
    predir=predir.replace(os.getenv("yeartag"),"")
    if not os.path.exists(predir):
        os.system("mkdir " + predir)
    if not os.path.exists(dirpath):
        os.system("mkdir " + dirpath)

binDir = os.getenv("LQANALYZER_DIR")+ "/bin/"
fr = open(binDir + 'Branch.txt', 'r')
sline=0
for line in fr:
    if sline ==0:
        print "Local code uses CATAnalyser version : " +  line
        print "Look up code at https://github.com/jalmond/LQanalyzer/tree/" + line
    if sline == 1:
        print "..."
    sline=sline+1

print "Running on tag : " +  os.getenv("CATTAG")    
   

if os.path.exists("LQRun/Macros/"):
    os.system("rm -r LQRun/Macros/")
if os.path.exists("LQRun/job_output/"):
    print "Cleaning up directory that failed to be removed by git merge"
    os.system("rm -r LQRun/job_output/")
if os.path.exists("LQRun/runJob_1.C"):
    print "Cleaning up file that failed to be removed by git merge"
    os.system("rm LQRun/runJob_1.C")
if os.path.exists("LQCycle/"):
    print "Cleaning up directory that failed to be removed by git merge"
    os.system("rm -r LQCycle/")

tag_dir  = os.getenv("LQANALYZER_LIB_PATH")+ "/" + os.getenv("CATTAG");
yeartag= str(os.getenv("yeartag"))

localfiledir = os.getenv("LQANALYZER_FILE_DIR")
snufiledir = os.getenv("FILEDIR")
snulumifiledir = os.getenv("LQANALYZER_DIR")+ "/data/Luminosity/"+yeartag
snufakefiledir = os.getenv("LQANALYZER_DIR")+ "/data/Fake/"+yeartag
snutriggerfiledir = os.getenv("LQANALYZER_DIR")+ "/data/Trigger/"+yeartag
snupileupfiledir= os.getenv("LQANALYZER_DIR")+ "/data/Pileup/"+yeartag
snuidfiledir= os.getenv("LQANALYZER_DIR")+ "/data/ID/"+yeartag
snubtagfiledir = os.getenv("LQANALYZER_DIR")+ "/data/BTag/"+yeartag
rochdir=os.getenv("LQANALYZER_DIR")+ "/data/rochester/"+yeartag

txtfiledir = os.getenv("LQANALYZER_DIR")+ "/LQRun/txt/"
old_lib_slc5=os.getenv("LQANALYZER_DIR")+ "/LQLib/slc5/"
old_lib_slc6=os.getenv("LQANALYZER_DIR")+ "/LQLib/slc6/"
old_lib_machine_1=os.getenv("LQANALYZER_DIR")+ "/LQLib/slc5_cms2/"
old_lib_machine_2=os.getenv("LQANALYZER_DIR")+ "/LQLib/slc5_cms3/"
old_lib_machine_3=os.getenv("LQANALYZER_DIR")+ "/LQLib/slc5_cms4/"
old_lib_machine_4=os.getenv("LQANALYZER_DIR")+ "/LQLib/slc6_cms5/"
old_lib_machine_5=os.getenv("LQANALYZER_DIR")+ "/LQLib/slc6_cms6/"
old_lib_machine_6=os.getenv("LQANALYZER_DIR")+ "/LQLib/slc6_cms1/"

if not os.path.exists(tag_dir):

    libpath=os.getenv("LQANALYZER_LIB_PATH")
    if os.path.exists(libpath):
        os.system("rm -r " + libpath)
        os.system("mkdir  " +libpath)
    os.system("mkdir " + tag_dir)

    print "Copying all latest rootfiles for use in analysis"

    if not os.path.exists(os.getenv("LQANALYZER_DIR")+ "/data/Luminosity/80X/") or not os.path.exists(os.getenv("LQANALYZER_DIR")+ "/data/Luminosity/76X/"):
        os.system("rm -r " + os.getenv("LQANALYZER_DIR")+ "/data/")
    MakeDirectory((os.getenv("LQANALYZER_DIR")+ "/data/"))
    MakeDirectory(snulumifiledir)
    os.system("cp " + localfiledir + "/Luminosity/*"+str(os.getenv("CATVERSION"))+".txt " + snulumifiledir)
    MakeDirectory(snufakefiledir)
    os.system("cp " + localfiledir + "/Fake/*.root " + snufakefiledir)
    MakeDirectory(snutriggerfiledir)
    os.system("cp " + localfiledir + "/Trigger/*.root " + snutriggerfiledir)
    os.system("cp " + localfiledir + "/Trigger/*.txt " + snutriggerfiledir)
    MakeDirectory(snupileupfiledir)
    os.system("cp " + localfiledir + "/Pileup/*.root "+ snupileupfiledir)
    MakeDirectory(snuidfiledir)
    os.system("cp " + localfiledir + "/ID/*.root " + snuidfiledir)
    if os.path.exists(snubtagfiledir):
        os.system("rm -r " + snubtagfiledir)
    MakeDirectory(snubtagfiledir)
    os.system("cp " + localfiledir + "/BTag/*.csv " + snubtagfiledir)
    MakeDirectory(rochdir)
    os.system("cp -r " + localfiledir + "/rochester/rcdata.2016.v3/ " + rochdir) 

    if os.path.exists(snufiledir+"/cMVAv2.csv"):
        os.system("rm  "+snufiledir+"/*.csv")
    if os.path.exists(snufiledir +"/triggers_catversion2016_802.txt") or os.path.exists(snufiledir +"/lumi_catversion2016_802.txt"):
        os.system("rm " + snufiledir+"/*.txt")
    if os.path.exists(snufiledir +"/Luminosity/triggers_catversion2016_802.txt"):
        os.system("rm " +snufiledir +"/Luminosity/*2016*")

    if os.path.exists(snufiledir):
        os.system("rm -r " + snufiledir)


    logdir =  os.getenv("LQANALYZER_LOG_8TeV_PATH")
    if os.path.exists(logdir):
        os.system("rm -r "+logdir)

    old_out=os.getenv("LQANALYZER_DIR")+"/data/output/CAT/"

    
    mount_name="/data2"
    if "cmscluster.snu.ac.kr" in str(os.getenv("HOSTNAME")):
        mount_name="/data4"

    new_out=mount_name+"/CAT_SKTreeOutput/JobOutPut/"+os.getenv("USER")
    print "cleaning up home directory"
    if not os.path.exists(new_out):
        os.system("mkdir " + new_out)

    new_out=mount_name+"/CAT_SKTreeOutput/JobOutPut/"+os.getenv("USER")+"/LQanalyzer/"
    if not os.path.exists(new_out):
        os.system("mkdir " + new_out)
        new_out=mount_name+"/CAT_SKTreeOutput/JobOutPut/"+os.getenv("USER")+"/LQanalyzer/data/"
        os.system("mkdir " + new_out)
        new_out=mount_name+"/CAT_SKTreeOutput/JobOutPut/"+os.getenv("USER")+"/LQanalyzer/data/output/"
        os.system("mkdir " + new_out)
        new_out=mount_name+"/CAT_SKTreeOutput/JobOutPut/"+os.getenv("USER")+"/LQanalyzer/data/output/CAT/"
        os.system("mkdir " + new_out)
        os.system("mv "+ old_out + "/* " + new_out)
        print "Moving output to " + new_out
        if os.path.exists(old_out):
            os.system("rm -r " + os.getenv("LQANALYZER_DIR")+"/data/output/")
            
        
    if os.path.exists(old_lib_slc5):
        os.system("rm -r " + old_lib_slc5)
    if os.path.exists(old_lib_slc6):
        os.system("rm -r " + old_lib_slc6)
    
    if os.path.exists(old_lib_machine_1):
        os.system("rm -r " + old_lib_machine_1)
    if os.path.exists(old_lib_machine_2):
        os.system("rm -r " + old_lib_machine_2)
    if os.path.exists(old_lib_machine_3):
        os.system("rm -r " + old_lib_machine_3)
    if os.path.exists(old_lib_machine_4):
        os.system("rm -r " + old_lib_machine_4)
    if os.path.exists(old_lib_machine_5):
        os.system("rm -r " + old_lib_machine_5)
    if os.path.exists(old_lib_machine_6):
        os.system("rm -r " + old_lib_machine_6)

    print "using branch for first time: All codes are being recompiled"
    os.system("source bin/Make/make_clean_newbranch.sh")
    

fakelib = os.getenv("LQANALYZER_LIB_PATH") + "/libHNCommonLeptonFakes.so"

if not os.path.exists(fakelib):
    os.system("source bin/Make/make_fake_lib.sh")



