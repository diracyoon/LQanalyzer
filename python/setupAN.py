import os,time

version = "Version1"
os.chdir("/data7/Users/jalmond/HNDileptonAN/2016_DiLep/myDir/")
os.system("bash setup.sh")
os.chdir("/data7/Users/jalmond/HNDileptonAN/2016_DiLep/myDir/notes/AN-17-999/trunk")

os.system("bash compile.sh")

if os.path.exists("index.html"):
    os.remove("index.html") 
file_html= open("index.html","w")
file_html.write('<embed src="AN-17-999_'+version+'.pdf" width="800px" height="2100px" />')
file_html.close()

date_tag=time.strftime("%d-%m-%Y")
if os.path.exists("index_"+date_tag+".html"):
    os.remove("index_"+date_tag+".html")
file_html= open("index_"+date_tag+".html","w")
file_html.write('<embed src="AN-17-999_'+version+'_'+date_tag+'.pdf" width="800px" height="2100px" />')
file_html.close()

catpath=os.getenv("LQANALYZER_DIR")+"/bin/catconfig"
readcatpath=open(catpath,"r")
lxmachine=""
for rline in readcatpath:
    if "localcpu" in rline:
        srline = rline.split()
        lxmachine=srline[2]
readcatpath.close()


os.system("ssh jalmond@" + lxmachine +".cern.ch mkdir /afs/cern.ch/user/j/jalmond/www/HNDiLepton2016_AN/"+version)
os.system("scp  index.html jalmond@" + lxmachine +".cern.ch:/afs/cern.ch/user/j/jalmond/www/HNDiLepton2016_AN/"+version)
os.system("scp  index_"+date_tag+".html jalmond@" + lxmachine +".cern.ch:/afs/cern.ch/user/j/jalmond/www/HNDiLepton2016_AN/"+version)
os.system("scp  /data7/Users/jalmond/HNDileptonAN/2016_DiLep/myDir/notes/tmp/AN-17-999_temp.pdf jalmond@" + lxmachine +".cern.ch:/afs/cern.ch/user/j/jalmond//www/HNDiLepton2016_AN/"+version + "/AN-17-999_"+version+".pdf")
os.system("scp  /data7/Users/jalmond/HNDileptonAN/2016_DiLep/myDir/notes/tmp/AN-17-999_temp.pdf jalmond@" + lxmachine +".cern.ch:/afs/cern.ch/user/j/jalmond//www/HNDiLepton2016_AN/"+version + "/AN-17-999_"+version+"_"+date_tag+".pdf")
print "https://jalmond.web.cern.ch/jalmond/HNDiLepton2016_AN/"+version + "/index.html"
print "https://jalmond.web.cern.ch/jalmond/HNDiLepton2016_AN/"+version + "/index_"+date_tag+".html"

if not os.path.exists("/data7/Users/jalmond/HNDileptonAN/2016_DiLep/BackUp/"+date_tag):
    os.mkdir("/data7/Users/jalmond/HNDileptonAN/2016_DiLep/BackUp/"+date_tag)
    os.system("scp -r /data7/Users/jalmond/HNDileptonAN/2016_DiLep/myDir/notes/AN-17-999/trunk /data7/Users/jalmond/HNDileptonAN/2016_DiLep/BackUp/"+date_tag)



