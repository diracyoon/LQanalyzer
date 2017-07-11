import os,sys

def  MakeHistFile(histlist, binlist,xminlist, xmaxlist,jobname):

    if not os.path.exists(str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname ):
        os.system("mkdir " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname)

    if not os.path.exists(str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname +"/dat/"):
        os.system("mkdir " + str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname +"/dat/")
    file_cuts=open(str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/"+jobname +"/dat/histfile.txt","w")
    for i in range(0,len(histlist)):
        file_cuts.write(histlist[i] + " " + binlist[i] + " " + xminlist[i] + " " +xmaxlist[i]+"\n")
    file_cuts.write("END")
    file_cuts.close()

    return


if os.getenv("USER") == "jalmond":
    os.system("cat ~/.ssh/config > check_connection.txt")
    
    ch_connect = open("check_connection.txt",'r')
    
    cpath="/tmp/"
    for line in ch_connect:
        if "ControlPath" in line:
            if "~/ssh" in line:
                cpath="~/"
            elif "/tmp/" in line:
                cpath="/tmp/"
            else:
                print "Modify the cms21 connection since  ControlPath in ~/.ssh/cofig is set to something other than tmp or home dir"
            
    ch_connect.close()
    os.system("rm check_connection.txt")

    
    os.system("ls " + cpath + " > check_snu_connection.txt")
    snu_connect = open("check_snu_connection.txt",'r')
    connected_lxplus=False

    for line in snu_connect:
        if "ssh-"+os.getenv("USER")+"@lxplus" in line:
            connected_lxplus=True
            
    os.system("rm check_snu_connection.txt")    
    if connected_lxplus == False:    
        print "No connection to lxplus please make connection in screen and run script again"
        quit()



from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--i", dest="i", default="123",help="tag")
parser.add_option("-x", "--x", dest="x", default="TEST",help="tag")

(options, args) = parser.parse_args()
configinputfile = options.i
jobdir = options.x

if not os.path.exists(configinputfile):
    print configinputfile + " is missing"
    sys.exit()

inputfile=""
histlist=[]
binlist=[]
xminlist=[]
xmaxlist=[]
cutlist=[]
skim=""
periodtag=""
analyzer=""
stream=""
isblind="false"
plottag="Default"
outputlist=[]
input_configfile = open(configinputfile,"r")
for line in input_configfile:
    if "*********************" in line:
        continue
    if "cut " in line and "# " in line:
        sline = line.split()
        cut = sline[2]
        cutlist.append(cut)
    elif "blinded" in line and "# " in line:
        sline = line.split()
        if sline[2] == "true":
            isblind="true"
    elif "samples" in line and "# " in line:
        sline = line.split()
        inputfile = str(os.getenv("LQANALYZER_DIR")) + "/Macros/CatPlotter/PlotConfig/"+sline[2]
    elif "skim" in line and "# " in line:
        sline = line.split()
        skim  = sline[2]
        print "skim = " + skim
    elif "periodtag" in line:
        sline = line.split()
        periodtag  =sline[2]
    elif "analyzer" in line:
        sline = line.split()
        analyzer  =sline[2]
    elif "plottag" in line:
        sline = line.split()
        plottag  =sline[2]
    elif "stream" in line:
        sline = line.split()
        stream  =sline[2]
    elif "caption" in line:
        sline = line.replace("# caption","")
        sline = sline.replace("...",'.";')
        cap= open("caption.txt","w")
        cap.write(sline)
        cap.close()
    elif "## h" in line:
        sline = line.split()
        histlist.append(sline[1])
        binlist.append(sline[2])
        xminlist.append(sline[3])
        xmaxlist.append(sline[4])

    elif "#################" in line:
        if len(cutlist) == 0:
            continue
        inputdir="/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/"+ analyzer +"/"+periodtag
        MakeHistFile(histlist, binlist,xminlist, xmaxlist,jobdir)

        if not os.path.exists("~/CATAnalyzerPlots/"):
            os.system("mkdir ~/CATAnalyzerPlots/")
        outputlist.append("https://jalmond.web.cern.ch/jalmond/SNU/CATAnalyzerPlots/HN13TeV_" + analyzer+"_"+os.getenv("CATVERSION")+"_"+periodtag+"_"+cut+"/histograms/CAT2016_"+os.getenv("CATVERSION")+plottag+"/indexCMS.html")

    
        os.system('python  ' + os.getenv("LQANALYZER_DIR") + '/Macros/CatPlotter/setupplotter.py -i ' + inputfile + ' -d ' + inputdir + ' -x ' + jobdir + ' -s ' + stream + ' -a ' + analyzer + ' -S ' + skim + ' -p ' + periodtag + ' -C ' + cutlist[0] + ' -M ' + configinputfile + ' -c ' + plottag + ' -b ' + isblind)
        histlist=[]
        binlist=[]
        xminlist=[]
        xmaxlist=[]
        cutlist=[]

print "Plots are found in:"
for x in outputlist:
    print x
