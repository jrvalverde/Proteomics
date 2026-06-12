db=${1}

USE_CD_HIT=${USE_CD_HIT:-NO}
USE_USEARCH=${USE_USEARCH:-YES}

if [ "$USE_CD_HIT" = "YES" ] ; then
    bin/remove_redundance_cdhit.sh $db
    #mv ${db}.1c ${db}.1
elif [ "$USE_USEARCH" = "YES" ] ; then
    bin/remove_redundance_usearch.sh $db
    #mv ${db}.1u ${db}.1
fi
exit
