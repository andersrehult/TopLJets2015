#!/bin/bash

#parse input arguments
ERA=2017
while getopts "o:y:" opt; do
    case "$opt" in
        o) WHAT=$OPTARG
            ;;
        y) ERA=$OPTARG
            ;;
    esac
done


#check an operation has been given
if [ -z "$WHAT" ]; then
    echo "steerAnalysis.sh -o <TEST/SEL/MERGE/PLOT> [ -y 2016/2017/2017lowpu/2018 ] ";
    echo "   TEST          - test locally the code on a single file";
    echo "   SEL           - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "   MERGE         - merge output"
    echo "   PLOT          - make plots"
    exit 1; 
fi

#configuration parameters
queue=tomorrow
samples=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/summer2019/samples_${ERA}.json
outdir=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/summer2019/analysis_${ERA}
if [[ ${ERA} = "2016" ]]; then
    githash=0c522df
    eosdir=/store/cmst3/group/top/RunIIReReco/2016/${githash}
    dataeosdir=${eosdir}
    lumi=35882
    lumiUnc=0.025
    testtag=Data13TeV_2016B_SingleMuon
    testfile=${eosdir}/${testtag}/Chunk_0_ext0.root
elif [[ ${ERA} = "2017" ]]; then
    githash=ab05162
    eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
    dataeosdir=${eosdir}
    lumi=41367
    lumiUnc=0.025
    testtag=Data13TeV_2017B_SingleMuon
    testfile=${eosdir}/${testtag}/Chunk_0_ext0.root
elif [[ ${ERA} = "2017lowpu" ]]; then
    ERA=2017
    githash=ab05162
    eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
    dataeosdir=/store/cmst3/group/top/RunIIReReco/2017/newproton_calib/
    lumi=220
    lumiUnc=0.025
    testtag=Data13TeV_2017H_SingleMuon_v2
    testfile=${dataeosdir}/${testtag}/MiniEvents_101.root
fi

#run the operation required
RED='\e[31m'
NC='\e[0m'
case $WHAT in
    
    TEST )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
            -i ${testfile} --tag ${testtag} \
            -o testsel_${ERA}.root --genWeights genweights_${githash}.root \
            --njobs 1 -q local --debug \
            --era era${ERA} -m RunTopSummer2019;
        ;;
    
    SEL )
        baseOpt="--genWeights genweights_${githash}.root"
        baseOpt="${baseOpt} -o ${outdir} -q ${queue} --era era${ERA} -m RunTopSummer2019"
        baseOpt="${baseOpt} --only ${samples}";
	python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} -i ${eosdir};
        if [ "${eosdir}" != "{dataeosdir}" ]; then
	    python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} -i ${dataeosdir} --farmappendix data;
        fi
	;;

    CHECKSELINTEG )
        python scripts/checkLocalAnalysisInteg.py ../../../FARM${githash}/
        ;;

    MERGE )
	mergeOutputs.py ${outdir} True;
	;;

    PLOT )
	commonOpts="-i ${outdir} -l ${lumi} --mcUnc ${lumiUnc}"
	python scripts/plotter.py ${commonOpts} -j ${samples}    -O ${outdir}/plots --saveLog;
	;;

esac
