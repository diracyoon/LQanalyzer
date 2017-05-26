import os,getpass
from CleanUp import *

path_jobpre="/data1/"
if "tamsa2.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    path_jobpre="/data2/"

LQANALYZER_DIR= str(os.getenv("LQANALYZER_DIR"))
LQANALYZER_LOG= str(os.getenv("LQANALYZER_LOG_PATH"))

if not LQANALYZER_DIR == "None" :
	datadir = LQANALYZER_DIR + "/data/"	
	if not (os.path.exists(datadir)):
		print "This is the first time running LQAnalyzer in this location"
		print "Making data directory in $LQANALYZER_DIR"
		os.system("mkdir " + datadir)
        
	outfiledir= LQANALYZER_DIR +"/data/output/"
	lumifiledir= LQANALYZER_DIR +"/data/Luminosity/"+ str(os.getenv("yeartag"))
	if not os.path.exists(LQANALYZER_DIR +"/data/Luminosity/"):
		os.system("mkdir " +LQANALYZER_DIR +"/data/Luminosity/")
	if not os.path.exists(lumifiledir):
		os.system("mkdir " + lumifiledir)
	btagfiledir = LQANALYZER_DIR +"/data/BTag/"+ str(os.getenv("yeartag"))
	if not os.path.exists(LQANALYZER_DIR +"/data/BTag/"):
		os.system("mkdir " +LQANALYZER_DIR +"/data/BTag/")
	if not os.path.exists(btagfiledir):
		os.system("mkdir " + btagfiledir)

	if not (os.path.exists(outfiledir)):
		os.system("mkdir " + outfiledir)
		print "Making data/output directory in $LQANALYZER_DIR"


	EightTeVdataOne=path_jobpre+"" + getpass.getuser() + "/LQ_SKTreeOutput/"
	EightTeVdataTwo="/data2/" + getpass.getuser() + "/LQ_SKTreeOutput/"
	 
	if os.path.exists(os.getenv("LQANALYZER_DIR")+ "/nohup.out"):
		os.system("rm " +os.getenv("LQANALYZER_DIR")+ "/nohup.out")


        CleanUpLogs(path_jobpre+"LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/"+ getpass.getuser()+ "/")
        CleanUpJobLogs(LQANALYZER_LOG)
	if os.getenv("HOSTNAME") == "cms.snu.ac.kr":
            CleanUpLogs(path_jobpre+"CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs("/data2/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs("/data7/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs("/data8/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs(os.getenv("LQANALYZER_BATCHLIB_PATH"))
            CleanUpLogs(EightTeVdataOne)
            CleanUpLogs(EightTeVdataTwo)
        else:
            CleanUpLogs("/data4/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
        localfiledir = os.getenv("LQANALYZER_FILE_DIR")
	datasetfiledir = os.getenv("LQANALYZER_DATASETFILE_DIR")
	txtfiledir = os.getenv("LQANALYZER_DIR")+ "/LQRun/txt/"
	cltxtfiledir = os.getenv("LQANALYZER_DIR")+ "/LQRun/txt/Cluster/"
	seldir =os.getenv("LQANALYZER_DIR")+  "/CATConfig/SelectionConfig/"
	os.system("cp " + localfiledir + "/Luminosity/triggers_catversion_"+str(os.getenv("CATVERSION"))+"* "  + lumifiledir)
	os.system("cp " + localfiledir + "/Luminosity/lumi_catversion_"+str(os.getenv("CATVERSION"))+".txt "  + lumifiledir)
	os.system("cp " + datasetfiledir + "/list_all_mc_"+str(os.getenv("CATVERSION"))+".sh " + txtfiledir)
        if os.getenv("HOSTNAME") == "cms.snu.ac.kr":
            if not os.getenv("USER") == "jalmond":
                os.system("cp " + localfiledir + "/Selection/*.sel " + seldir)
        else:
            os.system("rm " + seldir  +"/electrons.sel") 
            while not os.path.exists(seldir  +"/electrons.sel"):
                os.system("scp 147.47.242.42:/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/Selection/*.sel  " + seldir)
                
	#os.system("cp " + localfiledir + "/*.csv " + btagfiledir)
	#os.system("source " +  os.getenv("LQANALYZER_DIR") + "/bin/IncludePrivateSamples.sh")
else:
	print "Area is not setup. Cannot make directories needed for analysis"

