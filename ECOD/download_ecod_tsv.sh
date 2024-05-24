
if [ -z "$1" ] 
  then 
    echo "Downloading ecod ver. 291..."
    curl prodata.swmed.edu/ecod/distributions/ecod.develop291.domains.txt > ecod_domains291.tsv

  else
    echo "Downloading ecod ver. $1..."
    curl prodata.swmed.edu/ecod/distributions/ecod.develop"$1".domains.txt > ecod_domains"$1".tsv
fi
