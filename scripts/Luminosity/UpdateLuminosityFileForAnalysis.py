import os,sys, filecmp

tmpcatversion=str(os.getenv("CATVERSION"))

newsamplelist=os.getenv("LQANALYZER_DATASETFILE_DIR") +"/datasets_snu_CAT_mc_"+tmpcatversion+"new.txt"
if(os.path.exists(newsamplelist)):
    os.system("chmod 777 " + newsamplelist)
    os.system("rm " + newsamplelist)

from EmailNewEffLumiList import *

def CheckForDuplicates(printDuplicates):

    samplelist=os.getenv("LQANALYZER_DATASETFILE_DIR") +"/datasets_snu_CAT_mc_"+catversion+".txt"
    os.system("chmod 777 " + samplelist)
    copy_samplelist=[]
    file_samplelist = open(samplelist,"r")
    for line in file_samplelist:
        copy_samplelist.append(line)
    file_samplelist.close()
    
    listsamples=[]
    duplicates=[]
    xsec_duplicates=[]
    rd_samplelist=os.getenv("LQANALYZER_DATASETFILE_DIR") +"/datasets_snu_CAT_mc_"+catversion+"_removeduplicates.txt"
    file_rd_samplelist= open(rd_samplelist,"w")
    for xline in copy_samplelist:
        split_current_line=xline.split()
        if len(split_current_line) == 6:
            for x in listsamples:
                if x == split_current_line[0]:
                    duplicates.append(x)
                    xsec_duplicates.append(split_current_line[3])
                    if printDuplicates==1:
                        print "Sample " + x + " is duplicated..... Fixing --> " + rd_samplelist
                    if printDuplicates==2:
                        print "Sample " + x + " is duplicated but xsec is different..... Please Fix by hand"

            listsamples.append(split_current_line[0])
        
    copied_duplicates=[]
    copied_xsecduplicates=[]
    for cl in copy_samplelist:
        split_current_line=cl.split()
        AddToFile=True
        csamples=""
        xseccsamples=""
        if len(split_current_line) == 6:
            csamples=split_current_line[0]
            xseccsamples=split_current_line[3]
            for x in duplicates:
                if x in cl:
                    if float(xseccsamples) == 1.:
                        AddToFile=False
            for c in range(0,len(copied_duplicates)):
                if copied_duplicates[c] == csamples and float(copied_xsecduplicates[c]) == float(xseccsamples):
                    AddToFile=False
        if AddToFile:
            file_rd_samplelist.write(cl)
            if csamples:
                copied_duplicates.append(csamples)
                copied_xsecduplicates.append(xseccsamples)
            
    file_rd_samplelist.close()
    
    if printDuplicates==1:
        print "Replacing " + samplelist + " with " + rd_samplelist
        os.system("rm " + samplelist)
        os.system("mv " + rd_samplelist + " " + samplelist)
    else:
        os.system("chmod 777 " + rd_samplelist)
        os.system("rm " + rd_samplelist)

def UpdateLumiFile(modlistpath, catversion,NewSampleList,whichfile):

    ### xseclist should contain lines that are updated in xsec
    ### samplelist should contain lines for new samples
    samplelist=os.getenv("LQANALYZER_DATASETFILE_DIR") +"/datasets_snu"+whichfile+"_CAT_mc_"+catversion+".txt"
    os.system("chmod 777 "  + samplelist)
    newsamplelist=os.getenv("LQANALYZER_DATASETFILE_DIR") +"/datasets_snu"+whichfile+"_CAT_mc_"+catversion+"new.txt"                                                               
    #samplelist="/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/datasets_snu_CAT_mc_"+catversion+".txt"
    #newsamplelist="/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/datasets_snu_CAT_mc_"+catversion+"tmp.txt"
    
    ### Make a copy of the original dataset list
    copy_samplelist=[]
    file_samplelist = open(samplelist,"r")
    for line in file_samplelist:
        copy_samplelist.append(line)
    file_samplelist.close()

    if True:
        file_samplelist = open(newsamplelist,"w")
        ### Make a new tmp list
        for xline in copy_samplelist:
            split_current_line=xline.split()
            replace_line=False
            if len(split_current_line) == 6:
                modlist=open(modlistpath,"r")
                for line in modlist:
                    split_modline=line.split()
                    if len(split_modline) == 6:
                        if split_modline[0] == split_current_line[0]:
                            replace_line=True
                            print "replace line : " + xline
                modlist.close()
            if not replace_line:
                file_samplelist.write(xline)
            else:
                modlist=open(modlistpath,"r")
                for line in modlist:
                    split_modline=line.split()
                    if len(split_modline) == 6:
                        if split_modline[0] == split_current_line[0]:
                            file_samplelist.write(line)
                            print "with " + line
                modlist.close()
        file_samplelist.close()    

   ### Make a copy of the original dataset list                                                                                                                                                                                                                                
    copy_samplelist=[]
    file_samplelist = open(samplelist,"r")
    for line in file_samplelist:
        copy_samplelist.append(line)
    file_samplelist.close()


    if len(NewSampleList)>0 :
        file_samplelist = open(newsamplelist,"w")
        for xline in copy_samplelist:
            isNewSample=False
            if "#### CATTuples" in xline:
                file_samplelist.write(xline)
                modlist=open(modlistpath,"r")
                for line in modlist:
                    isNewSample=False

                    split_modline=line.split()
                    if len(split_modline) == 6:
                        for xn in NewSampleList:
                            if (split_modline[0] + " ")  in xn or split_modline[1] in xn:
                                isNewSample=True
                    if not "SKTree" in line:
                        if isNewSample:
                            file_samplelist.write(line)
                modlist.close()
            elif "#### Single" in xline:
                file_samplelist.write(xline)
                modlist=open(modlistpath,"r")
                for line in modlist:
                    isNewSample=False
                    split_modline=line.split()
                    if len(split_modline) == 6:
                        for xn in NewSampleList:
                            xnsplit = xn.split()
                            if len(xnsplit) == 4:
                                if (xnsplit[0] + " ")  in line or xnsplit[1] in line:
                                    isNewSample=True

                    if "SKTree" in line:
                        if not "MCDiLep" in line:
                            if isNewSample:
                                file_samplelist.write(line)
                modlist.close()
            elif "#### Di" in xline:
                file_samplelist.write(xline)
                modlist=open(modlistpath,"r")
                for line in modlist:
                    isNewSample=False
                    split_modline=line.split()
                    if len(split_modline) == 6:
                        for xn in NewSampleList:
                            xnsplit = xn.split()
                            if len(xnsplit) == 4:
                                if (xnsplit[0] + "_dilep")  in line or xnsplit[1] in line:
                                    isNewSample=True
                                                                        
                                    
                                    
                    if "SKTree" in line:
                        if "MCDiLep" in line:
                            if isNewSample:
                                file_samplelist.write(line)
                modlist.close()
            elif "#### Tri" in xline:
                file_samplelist.write(xline)
                modlist=open(modlistpath,"r")
                for line in modlist:
                    isNewSample=False
                    split_modline=line.split()
                    if len(split_modline) == 6:
                        for xn in NewSampleList:
                            xnsplit = xn.split()
                            if len(xnsplit) == 4:
                                if (xnsplit[0] + "_trilep")  in line or xnsplit[1] in line:
                                    isNewSample=True



                    if "SKTree" in line:
                        if "MCTriLep" in line:
                            if isNewSample:
                                file_samplelist.write(line)
                modlist.close()

            elif "#### NoCut" in xline:
                file_samplelist.write(xline)
                modlist=open(modlistpath,"r")
                for line in modlist:
                    isNewSample=False
                    split_modline=line.split()
                    if len(split_modline) == 6:
                        for xn in NewSampleList:
                            xnsplit = xn.split()
                            if len(xnsplit) == 4:
                                if (xnsplit[0] + "_nocut")  in line or xnsplit[1] in line:
                                    isNewSample=True



                    if "SKTree" in line:
                        if "MCNoCut" in line:
                            if isNewSample:
                                file_samplelist.write(line)
                modlist.close()


            else:
                file_samplelist.write(xline)

        file_samplelist.close()


def CheckFileFormat(filepath):

    if not os.path.exists(filepath):
        print "File " + filepath + " does not exist"
        sys.exit()

    file_fulllist = open(filepath,"r")
    for line in file_fulllist:
        if not "END" in line:
            splitline = line.split()
            if not len(splitline) == 4:
                print "Error reading file " + filepath + " too many parts of line :"
                print line
                sys.exit()

    file_fulllist.close()

### ExtractListFromDatasetFile Makes a list of samples located in /data1/LQAnalyzer_rootfiles_for_analysis/DataSetLists//dataset-$CATVERSION

import ExtractListFromDatasetFile



catversion=str(os.getenv("CATVERSION"))

path_full_sample_list=os.getenv("LQANALYZER_DATASET_DIR")+"/cattuplist_"+catversion+".txt"
path_full_sample_list_user=os.getenv("LQANALYZER_DATASET_DIR")+"/"+ os.getenv("USER")  +"/cattuplist_"+catversion+ os.getenv("USER")+".txt"


CheckFileFormat(path_full_sample_list_user)



change_in_xsec=False
new_sample=False
new_catversion=False
newxsec_list=[]
newsample_list= []


if not os.path.exists(os.getenv("LQANALYZER_DATASETFILE_DIR") +"/datasets_snu_CAT_mc_"+catversion+".txt"):
    os.system("rm " + path_full_sample_list)

if os.path.exists(path_full_sample_list):

    check_file = open(path_full_sample_list,"r")
    correct_file_format=False
    for line in check_file:
        if "DY" in line:
            correct_file_format=True
    check_file.close()
    if not correct_file_format:
        os.system("rm " + path_full_sample_list)


if os.path.exists(path_full_sample_list):
    CheckFileFormat(path_full_sample_list)
    #### If list just made is identical to old list no update is done and code exits
    if filecmp.cmp(path_full_sample_list,path_full_sample_list_user):
        print "List of files in /data1/LQAnalyzer_rootfiles_for_analysis/DataSetLists/datasets_"+catversion + " is unchanged. No update needed"
        os.system("rm " + path_full_sample_list_user)
        CheckForDuplicates(1)
        sys.exit()
    else:
        
        
        #### check all samples that are in previous sample list are in new sample list
        #### the only reason they will not be is if a dataset file was removed or wrongly modified
        file_fulllist = open(path_full_sample_list,"r")
        removed_samples=[]
        for fullline in file_fulllist:
            line_exists=False
            s_fullline = fullline.split()
            s_fullline_1=""
            if len(s_fullline) > 0:
                s_fullline_1 = s_fullline[0]
            #### sline1 now = alias name in master list    
            file_user = open(path_full_sample_list_user,"r")
            for userline in file_user:
                s_userline = userline.split()
                s_userline_1=""
                if len(s_userline) > 0:
                    s_userline_1 = s_userline[0]
                #### scline2 now = alias name in new list 
                if s_fullline_1 == s_userline_1:
                    line_exists=True
                ### sample in master list exists in new list
            file_user.close()                    
            if not line_exists:
                removed_samples.append(fullline)    
                #### sample in master list does not  exist in new list 
        file_fulllist.close()
        
        if len(removed_samples) > 0:
            print "Some samples have been removed from previous sample list. These are:"
            for xline in removed_samples:
                print xline
            print "The code is exiting. Please add the samples back."
            sys.exit()


        ### check xsecs in new list compared to master list 
        file_fulllist = open(path_full_sample_list,"r")

        xsec=0.
        xsec_2=0.
        sample1=""
        sample2=""
        ### loop over maste list
        for line in file_fulllist:
            splitline=line.split()
            if len(splitline) == 1:
                ### we want to consider only lines with split size of 4 and remove size 1
                continue
            if len(splitline) == 4:
                xsec=float(splitline[2])
                sample1=splitline[0]
                ### set alias and xsec from master list
            file_user = open(path_full_sample_list_user,"r")
            for cline in file_user:
                splitcline=cline.split()
                ### we want to consider only lines with split size of 4 and remove size 1 
                if len(splitcline) == 1:
                    continue
                if len(splitcline) == 4:
                    xsec_2=float(splitcline[2])
                    sample2=splitcline[0]
                ### check if alias' are the same and then compare xsec's    
                if sample1 == sample2:
                    if not xsec_2 == xsec:
                        new_xsec=True
                        newxsec_list.append(cline)
            file_user.close()
        file_fulllist.close()
        

        ### check for new samples in new list compared to master list
        file_user = open(path_full_sample_list_user,"r")
        sample1=""
        sample2=""
        for line in file_user:
            line_exists=False
            splitline=line.split()
            if len(splitline) == 1:
                continue
            if len(splitline) == 4:
                sample1=splitline[0]

            file_cat = open(path_full_sample_list,"r")
            for cline in file_cat:
                splitcline=cline.split()
                if len(splitcline) == 1:
                    continue
                if len(splitcline) == 4:
                    sample2=splitcline[0]
                if sample1 == sample2:
                    line_exists=True
            file_cat.close()

            file_mc = open(os.getenv("LQANALYZER_DATASETFILE_DIR") +"/datasets_snu_CAT_mc_" + catversion + ".txt","r") 
            inMCfile=False
            for cline in file_mc:
                splitcline=cline.split()
                if len(splitcline) > 1:
                    if sample1 == splitcline[0]:
                        inMCfile=True
            file_mc.close()

            if not inMCfile:
                line_exists=False

            if not line_exists:
                new_sample=True
                newsample_list.append(line)
        file_user.close()

        #### check also datasetmc file and add missed samples to newsample list
        

        
        if not os.path.exists(os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/log"):
            os.system("mkdir "+ os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/log")
            
        path_newfile=os.getenv("LQANALYZER_DATASET_DIR")+"/"+os.getenv("USER")+"/newfile.txt"
        #### rename variables later
        path_newfile2=os.getenv("LQANALYZER_DATASET_DIR")+"/"+os.getenv("USER")+"/newfile2.txt"
        path_newfile3=os.getenv("LQANALYZER_DATASET_DIR")+"/"+os.getenv("USER")+"/newfile3.txt"
        file_newfile = open(path_newfile,"w")
        file_newfile2 = open(path_newfile2,"w")
        file_newfile3 = open(path_newfile3,"w")
        print "Lines to be updated:"
        for x in newxsec_list:
            print "xsec " + x
            newx=x.split()
            file_newfile.write(x)
            file_newfile2.write(x)
        for x in newsample_list:
            print x
            file_newfile.write(x)
            file_newfile3.write(x)
        file_newfile.close()
        file_newfile2.close()
        file_newfile3.close()

        if os.path.exists(os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_CAT_mc_" + catversion + "new.txt"):
            os.system("rm " + os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_CAT_mc_" + catversion + "new.txt")
        print "\n"
        os.system("source " + os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/runGetEffLumi.sh " + path_newfile + " new ")
        print "source " + os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/runGetEffLumi.sh " + path_newfile + " new "
        print "\n"
        
        file_exists=False
        if   os.path.exists(os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_nonsig_CAT_mc_" + catversion + "new.txt"):
            file_exists=True
        if   os.path.exists(os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_sig_CAT_mc_" + catversion + "new.txt"):
            file_exists=True


        if not  file_exists:
            print  os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/runGetEffLumi.sh was meant to produce file "+ os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_CAT_mc_" + catversion + "new.txt"
            print "This file does not exists: exiting...."
            sys.exit()
        else:     
            ### update lumifile:
            isnewsample= len(newsample_list) > 0
            
            UpdateLumiFile(os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_nonsig_CAT_mc_" + catversion + "new.txt", catversion, newsample_list,"_nonsig")
            UpdateLumiFile(os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_sig_CAT_mc_" + catversion + "new.txt", catversion, newsample_list,"_sig")
            UpdateLumiFile(os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_CAT_mc_" + catversion + "new.txt", catversion, newsample_list,"")

            os.system("rm " + os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_nonsig_CAT_mc_" + catversion + "new.txt")
            os.system("rm " + os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_sig_CAT_mc_" + catversion + "new.txt")
            os.system("rm " + os.getenv("LQANALYZER_DIR")+"/scripts/Luminosity/datasets_snu_CAT_mc_" + catversion + "new.txt")
            samplelist=os.getenv("LQANALYZER_DATASETFILE_DIR") +"/datasets_snu_CAT_mc_"+catversion+".txt"
            os.system("chmod 777 " + samplelist)
            newsamplelist=os.getenv("LQANALYZER_DATASETFILE_DIR") +"/datasets_snu_CAT_mc_"+catversion+"new.txt"
            print "Is the following list of differences correct:"
            print "\n"
            print "diff " + samplelist + " " + newsamplelist
            os.system("diff " + samplelist + " " + newsamplelist)
            print "\n"
            if not os.getenv("USER") =="jalmond":
                input = raw_input("If Yes : Type Y and Enter. (not typing Y will not update the file: ")
                if input == "Y":
                    print "replacing " + samplelist + " with  " + newsamplelist
                    os.chmod(newsamplelist, 0777)
                    os.system("cp " + newsamplelist + " " + samplelist)
                    os.system("chmod  -R 777 " + os.getenv("LQANALYZER_DATASETFILE_DIR") )
                    #os.chmod(samplelist, 0777)
                    print "replacing " + path_full_sample_list + " with " + path_full_sample_list_user
                    os.system("cp " + path_full_sample_list_user + " " + path_full_sample_list)
                    os.system("chmod 777 " + path_full_sample_list)
                    os.system("rm " + path_full_sample_list_user)
                else:
                    print "You ignored changes. The sample list will not be updated"
                    os.system("chmod 777 " + newsamplelist)
                    os.system("rm " + newsamplelist)
                    sys.exit()
                    
            else:
                print "replacing " + samplelist + " with  " + newsamplelist
                os.chmod(newsamplelist, 0777)
                print  "removing " + samplelist
                os.system("rm " + samplelist)
                print "replacing " + samplelist + " with " + newsamplelist
                os.system("cp " + newsamplelist + " " + samplelist)
                os.chmod(samplelist, 0777)
                print "replacing " + path_full_sample_list + " with " + path_full_sample_list_user
                os.system("cp " + path_full_sample_list_user + " " + path_full_sample_list)
                os.system("chmod 777 " + path_full_sample_list)
                os.system("rm " + path_full_sample_list_user)
                    
            os.system("chmod 777 " + samplelist)
            os.system("chmod 777 " + newsamplelist)
            os.system("rm " + newsamplelist)   
        
        print "Running runInputListMaker.sh: Note this may take several minutes..."    
        os.system("source " + os.getenv("LQANALYZER_DIR")+"/scripts/runInputListMaker.sh")
        ##### now check file has no duplicates
        CheckForDuplicates(1)
        CheckForDuplicates(2)

        list_new=[]
        if len(newsample_list) > 0:
            file_newlist = open(path_newfile3,"r")
            for line in file_newlist:
                sline = line.split()
                if len(sline) == 4:
                    list_new.append(sline[0])
            file_newlist.close()                                                                                                      
            os.system("cp " + os.getenv("LQANALYZER_DIR")+"/LQRun/txt/list_user_mc.sh " + os.getenv("LQANALYZER_DIR")+"/LQRun/txt/list_user_mctmp.sh")
            file_userlist = open(os.getenv("LQANALYZER_DIR")+"/LQRun/txt/list_user_mc.sh","a")
            addstring = "declare -a new_list=("
            runSKTreemaker=False
            for l in list_new:
                if os.path.exists("/data2/CatNtuples/"+str(os.getenv("CATVERSION"))+"/SKTrees/MC/"+str(l)):
                    print "Not remaking sktree as this already exists (name must have been deleted from list by  hand to recalculate lumi"
                else:
                    print "Sample " + str(l) + " is new. Making SKtree"
                    addstring+="'"+l+"' "
                    runSKTreemaker=True
            addstring+=")\n"
            print "Adding list "
            print str(addstring)
            file_userlist.write(addstring)
            file_userlist.close()
            if runSKTreemaker:
                os.system("bash " + os.getenv("LQANALYZER_DIR")+"/bin/submitSKTree.sh -M True -a SKTreeMaker -list new_list -c " + catversion + " -m ' first time sample is made in current catversion'")
                
            os.system("mv " +  os.getenv("LQANALYZER_DIR")+"/LQRun/txt/list_user_mctmp.sh " + os.getenv("LQANALYZER_DIR")+"/LQRun/txt/list_user_mc.sh")

        os.system("ls -l " + samplelist + " > testperm")
        checkperm  = open("testperm","r")
        permok=False
        for line in checkperm:
           if "-rwxrwxrwx" in line:
               permok=True
        checkperm.close()
        os.system("rm testperm")
        if not permok:
            print "####"*100
            print "----"*100
            print "Error in permission file: please type :  chmod 777 " + samplelist + " in terminal"
            print "Check permission is -rwxrwxrwx using ll -rth " + samplelist
            print "If permission is not changed please make this known"
            print "----"*100
            print "####"*100
            
        if os.getenv("USER") == "jalmond":
            sys.exit()
        if len(newxsec_list) > 0:
            EmailNewXsecList(catversion,path_newfile2)
        if len(newsample_list) > 0:
            EmailNewSampleList(catversion,path_newfile3)                
else:

    
    ### if sample list does not exist then this is first  time it is run with new catversion so cp new list to main list

    #path_full_sample_list=os.getenv("LQANALYZER_DATASET_DIR")+"/cattuplist_"+catversion+".txt"
    #path_full_sample_list_user=os.getenv("LQANALYZER_DATASET_DIR")+"/"+ os.getenv("USER")  +"/cattuplist_"+catversion+ os.getenv("USER")+".txt"
    
    
    print "NEW SAMPLES"
    
    os.system("cp " + path_full_sample_list_user + " " + path_full_sample_list)
    os.system("chmod 777 " + path_full_sample_list)
    new_catversion=True
    
    lqdir = str(os.getenv("LQANALYZER_DIR"))
    if not os.path.exists(lqdir+"/scripts/Luminosity/log"):
        os.system("mkdir "+ lqdir+"/scripts/Luminosity/log")
 
    print "source " + lqdir+"/scripts/Luminosity/runGetEffLumi.sh " + os.getenv("LQANALYZER_DATASET_DIR")+"/cattuplist_"+str(os.getenv('CATVERSION'))+".txt"
    
    #### runGetEffLumi.sh creates dataset input list for analysis
    ### This takes the list of datasets created above
    ### It collects samples at SNU and counts events to calculate effective luminosity of each sample
    os.system("source " + lqdir+"/scripts/Luminosity/runGetEffLumi.sh " + os.getenv("LQANALYZER_DATASET_DIR")+"/cattuplist_"+str(os.getenv('CATVERSION'))+".txt")

    print "Finished GetEff"
    if os.path.exists(lqdir+"/scripts/Luminosity/log"):
        os.system("rm -r "+lqdir+"/scripts/Luminosity/log")

    print "Running scripts/runInputListMaker.sh"
    if os.path.exists(lqdir+"/scripts/Luminosity/inputlist_efflumi.txt"):
        os.system("rm " + lqdir+"/scripts/Luminosity/inputlist_efflumi.txt")
    
    os.system("source " + os.getenv("LQANALYZER_DIR")+"/scripts/runInputListMaker.sh")

    os.system('bash ' + os.getenv('LQANALYZER_DIR')+'/bin/submitSKTree.sh -M True -a  SKTreeMaker -list all_mc  -c '+catversion+' -m "First set of cuts with '+catversion+'cattuples"')
    os.system('bash  ' + os.getenv('LQANALYZER_DIR')+'/bin/submitSKTree.sh -M True -a  SKTreeMakerDiLep -list all_mc  -c '+catversion+'  -m "First set of cuts with '+catversion+' cattuples"')
    os.system('bash  ' + os.getenv('LQANALYZER_DIR')+'/bin/submitSKTree.sh -M True -a  SKTreeMakerTriLep -list all_mc  -c '+catversion+'  -m "First set of cuts with '+catversion+' cattuples"')

    EmailNewList(catversion)    
