#!/bin/bash -l

if [ $# -lt 3 ]; then
    exit 1
fi

file=$1
chain=$2
out=$3

## file=Effect/ctg_ad.bed.gz
## chain=hg19ToHg38.over.chain.gz
## out=Effect_HG38/ctg_ad.bed.gz

mkdir -p $(dirname $out)
coord=$(dirname $out)/$(basename $out .gz)_coord.gz
lifted=$(dirname $out)/$(basename $out .gz)_lifted

mkdir -p /broad/hptmp/compbio/data/sort/$out/
temp=$(mktemp -d "/broad/hptmp/compbio/data/sort/$out/temp.XXXXXXXX")

zcat $file | tail -n+2 | awk -F'\t' '{ print $1 "\t" int($2) "\t" int($3) "\t" ($1 ":" int($2) "-" int($3)) }' | sort -T $temp -k1,1 -k2,2n | gzip > $coord || exit 1

./bin/liftOver ${coord} ${chain} ${lifted} ${lifted}_remove || exit 1

cat $file | gzip -d | awk -vLIFTED=${lifted} -F'\t' 'BEGIN {
  while((getline line < LIFTED) > 0) {
    split(line,arr,"\t");
    valid[arr[4]] = arr[1] "\t" arr[2] "\t" arr[3];
  }
}
NR == 1 { 
  print $0
}
NR > 1 && (($1 ":" int($2) "-" int($3)) in valid){
  printf "%s", valid[($1 ":" int($2) "-" int($3))];
  for(j=4; j<=NF; ++j) printf "\t" $j;
  printf "\n";
}' | sort -T $temp -k1,1 -k2,2n | bgzip -c > $out || exit 1

[ -d $temp ] && rm -r $temp
[ -f $lifted ] && rm ${lifted} ${lifted}_remove
[ -f $coord ] && rm $coord
