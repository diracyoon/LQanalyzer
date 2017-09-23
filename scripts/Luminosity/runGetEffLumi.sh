lumidir="/data2/LQ_SKTreeOutput/Lumi/"${USER}
if [[ ! -d "${lumidir}" ]]; then
    mkdir $lumidir
fi
cd $LQANALYZER_DIR/scripts/Luminosity/
root -b -l <<EOF
.L GetEffectiveLuminosity.C+g
GetEffectiveLuminosity("$1",true, "$2","$3");
GetEffectiveLuminosity("$1",false,"$2","$3");
.q;
EOF
#bash cleanup.sh
cd -
#python EmailNewEffLumiList.py