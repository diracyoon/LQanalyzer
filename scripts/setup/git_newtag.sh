if [[ $LQANALYZER_DIR == "" ]];
    then
    source $LQANALYZER_DIR/setup.sh
fi

########## Tag index
itag=".23"
tagname=$CATVERSION$itag


diff  $LQANALYZER_DIR/setup.sh $LQANALYZER_DIR/scripts/setup/tag_setup.sh >> SetupCheck.txt

setup_is_different="False"
while read line
do
    setup_is_different="True"
done  < SetupCheck.txt

rm SetupCheck.txt
if [[ $setup_is_different == "True" ]];
then
    echo "setup.sh is changed. Make changes to scripts/setup/tag_setup.sh"
    return
fi

notnew_tag=False
while read line
do
    if [[ $line  == *$tagname* ]]; then
	notnew_tag=True
    fi	
done  < /data1/LQAnalyzer_rootfiles_for_analysis/CATTag/LatestTag80X.txt

if [[ $notnew_tag == "False" ]];then
    echo "$tagname (HEAD)" >> LatestTag80X.txt
    
    while read line
    do
	if [[ $line != *"HEAD"* ]];then
	    echo "$line" >> LatestTag80X.txt
	else
	    suffix="(HEAD)"
	    sline=${line%$suffix}
	    echo "$sline" >> LatestTag80X.txt
	fi
    done  < /data1/LQAnalyzer_rootfiles_for_analysis/CATTag/LatestTag80X.txt
    mv LatestTag80X.txt /data1/LQAnalyzer_rootfiles_for_analysis/CATTag/LatestTag80X.txt
else
    echo "Not adding a new tag"
fi

rm $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "export CATVERSION="$CATVERSION >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "### If there is a small bug/new code then new subtag is made"  >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "export tag_numerator='"$itag"'"  >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "if [[ \$1 == '"branch"' ]];"  >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "    then" >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "    export CATTAG=" >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "else" >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "    export CATTAG=$CATVERSION"$itag >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "fi" >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "" >> $LQANALYZER_DIR/scripts/setup/SetBrachAndTag.sh

sendemail=false

tag_exists=false
git tag>> gitcheck.txt
while read line
  do
  if [[ $line  == *$tagname* ]];
      then
      tag_exists=true
  fi
done < gitcheck.txt

echo $tag_exists
rm gitcheck.txt


cp $LQANALYZER_DIR/scripts/setup/tag_setup.sh $LQANALYZER_DIR/setup.sh
git commit -a

#"New Tag: "$tagname  

if [[ $1 == "" ]];
    then
    echo "Making new tag "$tagname
fi



if [[ $tag_exists == "true" ]];
    then
    echo "Deleting branch "$tagname
    git tag -d $tagname
    git push origin :refs/tags/$tagname
    
fi

if [[ $1 == "-e" ]];
    then
    sendemail=true;
fi


git tag $tagname
git push --tags

if [[ $sendemail == "true" ]]; then
    python $LQANALYZER_DIR/python/NewTagEmail.py -t $tagname
fi