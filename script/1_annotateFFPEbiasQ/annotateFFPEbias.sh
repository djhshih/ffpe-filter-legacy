#!/bin/sh

source /broad/software/scripts/useuse
use Python-2.7 
use Java-1.7

while getopts "i:b:c:r:f:t:d:o:x:s:j:J:?" Option
do
    case $Option in
        i    ) ID=$OPTARG;;
        b    ) BAM=$OPTARG;;
        c    ) CTX=$OPTARG;;
        r    ) REF=$OPTARG;;
        f    ) FROM=$OPTARG;;
        t    ) TO=$OPTARG;;
        d    ) DBSNP=$OPTARG;;
        o    ) OUT=$OPTARG;;
        x    ) EXT=$OPTARG;;
        s    ) MAX=$OPTARG;;
        j    ) JAR=$OPTARG;;
        J    ) JAVA=$OPTARG;;
        ?    ) echo "Invalid option: -$OPTARG" >&2
               exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               exit 0;;
    esac
done


# Have script fail if error is created in subtask.
set -e

echo "ID:       	${ID}"
echo "Bam:         	${BAM}" 
echo "Reference:	${REF}" 
echo "dbSNP:    	${DBSNP}" 
echo "Context:     	${CTX}" 
echo "FromBase:     	${FROM}" 
echo "ToBase:     	${TO}" 
echo "extension:    	${EXT}" 
echo "max sites:    	${MAX}" 
echo "Output area:	${OUT}" 


oopt=" -o $OUT "
if [[ -z $OUT ]]; then
   echo "no output full path"
   oopt=""
   OUT="."
fi   

copt=" -c $CTX "
if [[ -z $CTX ]]; then
   echo "no context option, use .CG  "
   copt=" -c .CG "
fi   

xopt=" -x ${EXT}.txt "
if [[ -z $EXT ]]; then
   echo "no extension "
   EXT="metrics"
   xopt=" -x metrics.txt "
fi   


if [[ -z $MAX ]]; then
   echo "no maximum sites specified use 2147483647 "
   MAX = "2147483647"
fi   


# parse bam path to construct oxog_metric path 
bamdir=$(dirname $BAM)
bamfile=$(basename $BAM)
stub="${bamfile%.*}"
ffpeQ=${stub}.${EXT}
    
    
Dir=`dirname $0`


if [ ! -f $bamdir/$ffpeQ ]; then

    echo "generate metrics file from bam"
        
    ${JAVA} -Xmx3600M -jar ${JAR} INPUT=${BAM} OUTPUT=${OUT}/${ffpeQ} REF_BASE=${FROM} ART_BASE=${TO} REFERENCE_SEQUENCE=${REF} DB_SNP=${DBSNP} MINIMUM_QUALITY_SCORE=20 MINIMUM_MAPPING_QUALITY=30 MINIMUM_INSERT_SIZE=60 MAXIMUM_INSERT_SIZE=600 USE_OQ=true CONTEXT_SIZE=1 STOP_AFTER=${MAX} VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
    if [[ $? -ne 0 ]] ; then
        exit 1
    fi
else
	echo "link to existing oxog_metrics file with bam"
	ln -s  $bamdir/${ffpeQ} ${OUT}/.

fi

# oxog_metric file in output area
FFPEQ=${OUT}/${ffpeQ}

ls -latr ${FFPEQ}



echo ""
echo "annotateFFPE  command line: "
echo "python $Dir/annotateFFPEbias.py -i $ID -m $FFPEQ  $oopt $copt $xopt"

python $Dir/annotateFFPEbias.py -i $ID -m $FFPEQ  $oopt $copt  $xopt

echo ""

echo "original report count"

cut -f1-5,11-12 $FFPEQ

echo "FFPE Q: ${CTX}"
cat ${OUT}/${ID}.${EXT}.txt
echo " "
echo "done"
