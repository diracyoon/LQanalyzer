import os

def CheckPathInFile(channel,sample,mc):
    CopySKTrees(channel, sample, mc , "False")
    

    

def CopySKTrees(channel,sample,mc,docopy):
    
    if not docopy == "False":
        docopy = "True"
    
    import os
    filename = 'txt/datasets_cms21.txt'
    tmpxsec=0
    inputdir=""
    if not mc:
        for line in open(filename, 'r'):
            if not line.startswith("#"):
                entries = line.split()
                if len(entries)==3:
                    if channel ==entries[0] and sample == entries[1]:
                        inputdir = entries[2]
                        print inputdir 
    else:
        for line in open(filename, 'r'):
            if not line.startswith("#"):
                entries = line.split()
                if len(entries)==3:
                    tmpxsec=entries[1]
                    if sample == entries[0]:
                        inputdir = entries[2]

    sample_dir = os.getenv("LQANALYZER_DIR")+ "/data/input/"
    if not os.path.exists(sample_dir):                
        os.system("mkdir " + sample_dir)
    if not mc:
        sample_dir+= "data/"
        if not os.path.exists(sample_dir):
            os.system("mkdir " + sample_dir)
        sample_dir += channel+sample
        if not os.path.exists(sample_dir):
                os.system("mkdir " + sample_dir)
    else: 
        sample_dir += "mc/"    
        if not os.path.exists(sample_dir):  
            os.system("mkdir " + sample_dir)
        sample_dir += sample
        if not os.path.exists(sample_dir):
            os.system("mkdir " + sample_dir)     

    if docopy == "True":
        user = raw_input("Please enter username of cms21 account")    
        run = "scp " + str(user) + "@147.47.242.42:" + inputdir + "/*.root " + sample_dir
        print run
        os.system(run)
    
    infile=False
    if not mc:    
        for line in open('txt/datasets_mac.txt', 'r'):
                if not line.startswith("#"):
                    entries = line.split()
                    if len(entries)==3:
                        if sample == entries[1]:    
                            if channel in entries[0]:
                                infile=True
    else:
        for line in open('txt/datasets_mac.txt', 'r'):
            if not line.startswith("#"):
                entries = line.split()
                if len(entries)==3:
                    if sample == entries[0]:
                       infile=True 

    if not infile:                   
        f = open('txt/datasets_mac.txt', 'a') 
        try:
            if not mc:
                towrite = channel + " " + sample + " " + sample_dir
            else:
                towrite = sample + " " + tmpxsec + " " + sample_dir
            f.write(towrite)
        finally:    
            f.closed()  


def make_batch_script(workdir, jname, lqdir, macroname, cluster):

    config='#!/bin/sh'
    config+='\n'
    config+='#$ -S /bin/bash \n'
    config+='#$ -N ' + jname + '  \n'
    config+='#$ -wd ' + workdir + ' \n'
    config+='#$ -o ' + workdir + ' \n'
    config+='#$ -e ' + workdir + ' \n'

    config+='echo "Job started at " `date` \n'
    if cluster:
        config+='cd /share/apps/root_v5_34_32/root/ \n'
    else:
        config+='cd /share/apps/root_v5-34-32/root/ \n'

    config+='. bin/thisroot.sh \n'
    config+='cd ' + lqdir + '\n'

    config+='source setup.sh ' + os.getenv("CATVERSION") + ' \n'
    config+='echo "PWD= "$PWD \n'
    config+='cd ' + workdir + ' \n'
    config+='root -l -q -b ' + macroname +' \n'
    config+='echo "Ran macro 2" \n'
    config+=''

    return config


def makeConfigFile(log,sample, input, tree, cycle, ver, output_tmp, output, nevents, outstep, skipev, datatype, channel, period, totalmcevents, xsec, tar_lumi, eff_lumi, useSKinput, runevent, libraries, runnp, runcf, runtau, skflag):

    config='{\n'
    config+='    gEnv->SetValue("TFile.AsyncPrefetching", 1);\n'

    config+='   //### Load Libraries\n'
    config+='   gSystem->Load("libSKTree.so");\n'
    config+='   gSystem->Load("libHist.so");\n'
    config+='   gSystem->Load("libAnalysisCore.so");\n'
    config+='   gSystem->Load("libNtuplecore.so");\n'
    config+='   gSystem->Load("libSelection.so");\n'
    config+='   gSystem->Load("libPlotting.so");\n'
    config+='   gSystem->Load("libHNCommonLeptonFakes.so");\n'
    config+='   gSystem->Load("librochcor2016");\n'
    config+='   gSystem->Load("libBTagSFUtil.so");\n'
    for lib in libraries:
        config+='   gSystem->Load("' + lib + ' + .so");\n'
        
    config+='   gSystem->Load("libAnalyzerTools.so");\n'
    config+='   gSystem->Load("libLQAnalysis.so");\n'
    config+='   gSystem->Load("libSKTreeMaker.so");\n'
    config+='   gSystem->Load("libValidation.so");\n'
    
    config+='   gSystem->Load("libPyROOT.so");\n'
    config+='   \n'
    config+='   TString filename = "' + input + '";\n'
    config+='   LQController analysis;\n'
    if useSKinput == "True":        
        config+='   analysis.SetLQInput(false);\n'
    elif useSKinput == "true":
        config+='   analysis.SetLQInput(false);\n'    
    else :
        config+='   analysis.SetLQInput(true);\n'
    config+='   analysis.SetJobName("' + sample + "_" + cycle+'");\n'
    config+='   analysis.SetInputList(TString(filename));\n'
    config+='   analysis.SetTreeName("'+ tree +'");\n'
    config+='   analysis.SetCycleName("' + cycle + '");\n'
    config+='   analysis.SetLogLevel("'+ log +'");\n'
    if not datatype =="mc" :
        config+='   analysis.SetChannel("'+ channel +'");\n'
    
    if not eff_lumi == -1.:
        config+='   analysis.SetEffectiveLuminosity(' + str(eff_lumi)+');\n'
    if not runevent ==-1:
        config+='   analysis.RunEvent(' + str(runevent)+');\n'
    if not nevents ==-1:
        config+='   analysis.SetNEventsToProcess(' + str(nevents) + ');\n'
    if not skipev ==-1:
        config+='   analysis.SkipEvents(' + str(skipev) + ');\n'
    if not outstep ==-1 :
        config+='   analysis.SetOutPutStep(' + str(outstep) + ');\n'
    if not period  =="":
        config+='   analysis.SetDataPeriod("' + period + '");\n'
    if not tar_lumi ==-1.:
        config+='   analysis.SetTargetLuminosity('+ str(tar_lumi) +');\n'
    if not datatype =="":
        config+='   analysis.SetDataType("' + datatype + '");\n'
    if not totalmcevents == -1 :
        config+='   analysis.SetTotalMCEvents(' + str(totalmcevents) +');\n'
    if not xsec == -1.:
        config+='   analysis.SetMCCrossSection(' + str(xsec) +');\n'
    config+='   analysis.SetUserFlag("'+skflag+'");\n'    
    config+='   analysis.RunNonPrompt("' +runnp+'");\n'
    config+='   analysis.RunTauMode("' +runtau+'");\n'
    config+='   analysis.RunChargeFlip("' +runcf+'");\n'
    config+='   analysis.SetName("' + sample + '",'+ str(ver) +',"'+ output_tmp +'");\n'                        
    config+='   analysis.Initialize();\n'
    config+='   analysis.ExecuteCycle();\n'

    config+='   gSystem->Exec("mv ' + output_tmp + sample + '_' + str(ver) + '.root ' + output + '");\n'
    config+='  \n}'

    return config

#Import date
def now():
    from datetime import datetime
    return str(datetime.now().month) + str(datetime.now().day)+ str(datetime.now().hour)+str(datetime.now().minute)
