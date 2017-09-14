import os,sys

def  MakeCutFile(cutlist, jobname):

    if not os.path.exists(str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname +"/dat/"):
        os.system("mkdir " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname +"/dat/")

    file_cuts=open(str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname +"/dat/cutfile.txt","w")
    file_cuts.write(cutlist+"\n")
    file_cuts.write("END")
    file_cuts.close()

    return



def  MakeConfFile(inputdir,jobname,datastream,analyzer,tag,tag2,catversion,periodtag,filedate,cut,leglist, plottag, isblind):

    if not os.path.exists(str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname +"/Config/"):
        os.system("mkdir " +  str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname +"/Config/")
    file_conf = open(str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname +"/Config/configfile.txt" ,"w"  )
    file_conf.write("# Specify_location_of_data/mc\n")
    file_conf.write("mcpath          "+inputdir+"/\n")
    file_conf.write("datapath          "+inputdir+"/"+analyzer+"_data_"+datastream+"_cat_"+catversion+".root\n")
    file_conf.write("datadrivenpath          "+inputdir+"/\n")
    file_conf.write("# Specify_the_list_of_hists_to_plot \n")
    file_conf.write("plottingpath    dat/histfile.txt \n")
    file_conf.write("cutpath         dat/cutfile.txt\n")
    file_conf.write("# file_prefix\n")
    file_conf.write("prefix          " + analyzer + tag+"\n")
    file_conf.write("postfix           " + tag2+"_cat_"+catversion+".root\n")
    file_conf.write("# NAME_THE_OUTPUT_DIR\n")
    file_conf.write("outputdir  HN13TeV_"+analyzer+"_"+catversion+"_"+periodtag+"_"+filedate+"_"+cut+"\n")
    file_conf.write("histdir    CAT2016_"+catversion+plottag+"\n")
    file_conf.write("# CONFIGURE_HISTOGRAMS\n")
    # CONFIGURE_HISTOGRAMS                                                                                                                                                  
    if isblind == "true":
        file_conf.write("showdata   false \n")
    else:
        file_conf.write("showdata   true \n")
    file_conf.write("usenp      false \n")
    file_conf.write("ylog       true \n")
    file_conf.write("# List_the_samples_to_be_includedv\n")
    for i in leglist:
         file_conf.write("samples    " + i+"\n")
    file_conf.write("END")
    return


from optparse import OptionParser
parser = OptionParser()
parser.add_option("-x", "--x", dest="x", default="TEST",help="tag")
parser.add_option("-b", "--b", dest="b", default="false",help="tag")
parser.add_option("-i", "--i", dest="i", default="123",help="tag")
parser.add_option("-d", "--d", dest="d", default="",help="tag")
parser.add_option("-D", "--D", dest="D", default="",help="tag")
parser.add_option("-s", "--s",  dest="s", default="",help="tag")
parser.add_option("-a", "--a",  dest="a", default="",help="tag")
parser.add_option("-S", "--S",  dest="S", default="",help="tag")
parser.add_option("-c", "--c",  dest="c", default="",help="tag")
parser.add_option("-p", "--p",  dest="p", default="",help="tag")
parser.add_option("-C", "--C",  dest="C", default="",help="tag")
parser.add_option("-M", "--M",  dest="M", default="",help="tag")
parser.add_option("-t", "--t",  dest="t", default="",help="tag")


(options, args) = parser.parse_args()
jobdir = options.x
inputdir = options.d
inputfile = options.i
stream = options.s
analyzer = options.a
skim = options.S
catversion=str(os.getenv("CATVERSION"))
period=options.p
filedate=options.D
cutlist=options.C
plottag=options.c
configinputfile=options.M
isblind=options.b
siginputfile=options.t

nodata=False
if isblind =="true":
    nodata=True


cap_file = open("caption.txt","r")
text_caption=""
for ca in cap_file:
    text_caption=ca
    break

histlist=[]
xtitlelist=[] 

input_configfile = open(configinputfile,"r")
for line in input_configfile:
    if "## h" in line:
        sline = line.split()
        histlist.append(sline[1])
        xtitle=""
        for r in range(5,len(sline)):
            xtitle+=sline[r]+" "
        xtitle=xtitle[:-1]
        xtitlelist.append( xtitle)            
input_configfile.close()

list_of_legends=[]
all_list_of_legends=[]
all_list_of_uncer=[]
list_of_legends_alias=[]
list_of_colors=[]


input_config = open(inputfile,"r")
for line in input_config:
    sline = line.split()
    if "#" in sline[0]:
        continue
    if "alias" in line:
        continue
    if "END" in line:
        break
    legend=""
    ni=0
    legalias=""
    for i in sline:
        if "|" in i:
            ni+=1
        elif ni == 1:
            legend=legend+i+" "
        elif ni == 2:
            legalias=i
    legend=legend[:-1]
    l_in_list=Truel_in_list=False
    for l in list_of_legends:
        if l == legend:
            l_in_list=True
    if not  l_in_list:
        list_of_legends.append(legend)
        list_of_legends_alias.append(legalias)
input_config.close()


empty_list=[]
empty_list_alias=[]

iterator_leg=-1
for i in list_of_legends:
    iterator_leg=iterator_leg+1
    list_of_samples=[]
    list_of_uncert=[]
    input_config = open(inputfile,"r")
    for line in input_config:
        if (" " + i) in line:
            sline = line.split()
            if "#" in sline[0]:
                continue

            color=""
            uncer=""
            ni=0
            for j in sline:
                if "|" in j:
                    ni+=1
                elif ni == 3:
                    color=j
                elif ni == 4:
                    uncer=j
            sample_exists=False        
            if(len(sline) ) > 0:
                tag="_"
                tag2=""
                if  skim == "SKTree_LeptonSkim":
                    tag="_SK"
                    tag2=""
                elif skim == "SKTree_DiLepSkim":
                    tag="_SK"
                    tag2="_dilep"
                elif skim == "SKTree_TriLepSkim":
                    tag="_SK"
                    tag2="_trilep"


                prefix= analyzer + tag
                postfix= tag2+"_cat_"+catversion+".root"

                if os.path.exists(inputdir+"/"+prefix+sline[0]+postfix):
                    sample_exists=True
                    list_of_samples.append(sline[0])
                else:
                    print "Missing " + inputdir+"/"+prefix+sline[0]+postfix

            if(len(sline) ) > 8:
                if sample_exists:
                    list_of_uncert.append(uncer)
    input_config.close()
    if len(list_of_samples) > 0:
        all_list_of_legends.append(list_of_samples)
        all_list_of_uncer.append(list_of_uncert)
        list_of_colors.append(color)
    else:
        empty_list.append(i)
        empty_list_alias.append(list_of_legends_alias[iterator_leg])
        
for i in empty_list:
    list_of_legends.remove(i)
for i in empty_list_alias:
    list_of_legends_alias.remove(i)

    
sigscales=[]
#sig1scale=0.
#sig2scale=0.
nsig=[]
#sig1=False
#sig2=False
sigpaths=[]
#sig1path=""
#sig2path=""
siglegsig=[]
#legsig1=""
#legsig2=""
ns=-1
drawsig=False
print "siginputfile = " + siginputfile

siginputfile="/home/jalmond/HeavyNeutrino/13TeV/LQAnalyzer_cat/LQanalyzer/Macros/CatPlotter/PlotConfig/"+siginputfile
siginput_configtmp = open(siginputfile,"r")

if siginputfile!="":
    drawsig=True
    siginput_config = open(siginputfile,"r")
    for line in siginput_config:
        if "END" in line:
            break
        if "###" in line :
            continue
        legsig=""
        ns=ns+1
        nsig.append(True)
        sline = line.split()
        nss=0
        for s in sline:
            if "|" == s:
                nss=nss+1
                continue
            if nss==0:
                sigpaths.append(s)
            if nss==1:
                sigscales.append(s)
            if nss==2:
                if s != "|":
                    legsig = legsig + s
                    legsig = legsig + " "
                    
        print legsig
        siglegsig.append(legsig)


plot_comfig_dir = str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/PlotConfig/"

path_macroC_file = str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir+"/MakeDataMCCompPlots.cpp" ;
path_macroH_file = str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir+"/MakeDataMCCompPlots.h" ;
path_setup = str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir+"/setenv.sh"
path_CMSconfig = str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir+"/CMS_lumi.h"


if not os.path.exists(str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir):
    os.system("mkdir " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir)
    os.system("mkdir " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir +"/Tables/")

os.system("cp " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/Code/setenv.sh " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir)
os.system("cp "+ str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/Code/MakeDataMCCompPlots.h "+ str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir)
os.system("cp "+ str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/Code/CMS_lumi.h  "+ str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir)
os.system("cp "+ str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/Code/Hist*  "+ str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir)
os.system("cp "+ str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/Code/make*  "+ str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir)

path_skeleton_macroC =  str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/Code//MakeDataMCCompPlots.cpp" ;
skeleton_macroC = open(path_skeleton_macroC,"r")
new_macroC = open(path_macroC_file,"w")
for line in skeleton_macroC:
    if "scp -r " in line:
        if os.getenv("USER") == "jalmond":
            new_macroC.write(line+"")
    elif "Open plots in" in line:
        if os.getenv("USER") == "jalmond":
            new_macroC.write(line+"")

    elif "vector<TString> legorder;" in line:
        new_macroC.write(line+"\n")
        for i in list_of_legends:
            new_macroC.write('legorder.push_back("' + i + '");\n')

    elif "  //removedata" in line:
        if nodata:
            new_macroC.write('hdata= dynamic_cast<TH1*>((h_nominal)->Clone((string(h_nominal->GetName())+"data").c_str()));;\n')
            new_macroC.write('hdata->GetXaxis()->SetRangeUser(xmin, xmax);;\n')
            new_macroC.write('hdata->GetYaxis()->SetRangeUser(ymin, ymax);;\n')
            new_macroC.write('SetTitles(hdata, hname);\n')
            new_macroC.write('showdata=false;\n')
              
              
    elif "//// SET TITLES" in line:
        new_macroC.write(line+"\n")
        for x in range(0,len(xtitlelist)):
            new_macroC.write('if(name.find("' + histlist[x] + '")!=string::npos) xtitle='+xtitlelist[x] +';\n')


    elif "/// SIGNAL" in line:
        new_macroC.write(line+"\n")
        new_macroC.write('')
        
        new_macroC.write('vector<TH1*> hsig ;\n' )   
        new_macroC.write('float int_bkg = hup->Integral()/2.; \n')
        for x in range (0, len(nsig)):
            if nsig[x]:
                new_macroC.write('TFile* file_sig'+str(x)+' =  TFile::Open(("'+sigpaths[x]+'")); \n')
                new_macroC.write('TH1* hsig'+str(x)+' = dynamic_cast<TH1*> ((file_sig'+str(x)+'->Get(name.c_str()))->Clone()); \n')
                new_macroC.write('hsig'+str(x)+'->Rebin(rebin); \n')
                new_macroC.write('hsig'+str(x)+'->Scale('+sigscales[x]+'); \n')
                new_macroC.write('FixOverUnderFlows(hsig'+str(x)+', xmax); \n')
                new_macroC.write('ymax = GetMaximum(hsig'+str(x)+', hsig'+str(x)+', ylog, name, xmax, xmin); \n')

                new_macroC.write('hsig'+str(x)+'->SetLineColor('+str(x+2)+'); \n')

                new_macroC.write('hsig'+str(x)+'->SetLineWidth(3.); \n')
                new_macroC.write('hsig'+str(x)+'->GetXaxis()->SetRangeUser(xmin,xmax); \n')
                new_macroC.write('hsig'+str(x)+'->GetYaxis()->SetRangeUser(ymin,ymax); \n')
            else:
                new_macroC.write('TH1* hsig1;\n');

            new_macroC.write('hsig.push_back(hsig'+str(x)+');\n')
            
    elif "bool drawsig" in line:
        if drawsig:
             new_macroC.write('bool drawsig=true;\n')
        else:
            new_macroC.write('bool drawsig=false;\n')
    elif "/// Draw sig" in line:
        new_macroC.write(line+"\n")
        for x in range (0, len(nsig)):
            new_macroC.write('hsigs['+str(x)+']->Draw("hist9same"); \n')
            new_macroC.write('legend->AddEntry(hsigs['+str(x)+'], "'+siglegsig[x]+'","l");\n')
            
            
    elif "/// Draw sig(1)" in line:
        new_macroC.write(line+"\n")
        for x in range (0, len(nsig)):
            new_macroC.write('hsigs['+str(x)+']->Draw("hist9same"); \n')

    elif "/// Draw(2) sig"in line:
        new_macroC.write(line+"\n")
        for x in range (0, len(nsig)):
            new_macroC.write('hsigs['+str(x)+']->Draw("hist9same"); \n')

                         
    elif "caption=" in line:
        new_macroC.write('caption="' + text_caption + '\n')

    elif "vector<pair<TString,float> > list;" in line:
        new_macroC.write("")
        new_macroC.write(line+"\n")
        for i in range(0,len(list_of_legends)):
            new_macroC.write('if(sample.Contains("'+list_of_legends_alias[i]+'")){')
            for j in range(0,len(all_list_of_legends[i])):
                new_macroC.write('    list.push_back(make_pair("' +str(all_list_of_legends[i][j])+'",'+all_list_of_uncer[i][j]+'));\n')
            new_macroC.write('}')
    elif "for( unsigned int i = 0; i < listofsamples.size(); i++){"       in line :
        new_macroC.write(line+"")
        for x in range(0,len(list_of_legends_alias)):
            if not "nonprompt" in list_of_legends_alias[x]:
                if not "chargeflip" in list_of_legends_alias[x]:
                    new_macroC.write('   if(listofsamples.at(i) =="'+list_of_legends_alias[x]+'")samples.push_back(make_pair(make_pair('+list_of_legends_alias[x]+','+list_of_colors[x]+'),"'+list_of_legends[x]+'"));\n')
                else:
                    new_macroC.write('   if(listofsamples.at(i) =="'+list_of_legends_alias[x]+'")samples.push_back(make_pair(make_pair(cf,'+list_of_colors[x]+'),"'+list_of_legends[x]+'"));\n')
            else:
                new_macroC.write('   if(listofsamples.at(i) =="'+list_of_legends_alias[x]+'")samples.push_back(make_pair(make_pair(np,'+list_of_colors[x]+'),"'+list_of_legends[x]+'"));\n')
                
    elif " /// Setup list of samples: grouped into different processes" in line:
        new_macroC.write(line+"")
        
        for x in range(0,len(list_of_legends_alias)):
            if not "nonprompt" in list_of_legends_alias[x]:
                if not "chargeflip" in list_of_legends_alias[x]:
                    new_macroC.write('vector<pair<TString,float> >  '+list_of_legends_alias[x]+' = InitSample(" '+list_of_legends_alias[x]+'"); \n')    
                else:
                    new_macroC.write('vector<pair<TString,float> > cf; \n')      
                    new_macroC.write('cf.push_back(make_pair("chargeflip",0.12)); \n')      

            else:
                new_macroC.write('/// NP is nonprompt \n')
                new_macroC.write('vector<pair<TString,float> > np;\n')
                new_macroC.write('np.push_back(make_pair("'+list_of_legends_alias[x]+'",0.34));\n')



    else:
        new_macroC.write(line+"")

skeleton_macroC.close()
new_macroC.close()


os.system("source  " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/Code/run.sh " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobdir)

MakeCutFile(cutlist,jobdir)
tag="_"
tag2=""
if  skim == "SKTree_LeptonSkim":
    tag="_SK"
    tag2=""
elif skim == "SKTree_DiLepSkim":
    tag="_SK"
    tag2="_dilep"
elif skim == "SKTree_TriLepSkim":
    tag="_SK"
    tag2="_trilep"

MakeConfFile(inputdir,jobdir,stream,analyzer,tag,tag2,catversion,period,filedate,cutlist,list_of_legends_alias, plottag, isblind)

os.system("source  " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/Code/runjob.sh " + jobdir)


#os.system("rm -r " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/" + jobdir)
