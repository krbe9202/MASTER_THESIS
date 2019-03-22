#!/bin/bash -l 

 # module load bioinfo-tools
 # module load BEDTools/2.27.1

peakPath=$1


# Run HOMER for motif analysis and peak annotation 
for peakFile in $peakPath/*.txt; do 
	motifPath="./homer/motifs/$(basename $peakFile | grep -oP '(?<=_).*.?(?=\.)')"
	annoPath="./homer/annotatePeaks/annotatePeaks_$(basename $peakFile | grep -oP '(?<=_).*.?(?=\.)')"

	findMotifsGenome.pl $peakFile mm9 $motifPath
	annotatePeaks.pl $peakFile mm9 > $annoPath.tsv
done  

for file in ./homer/annotatePeaks/*.tsv; do 
 
	awk -F '\t' -v OFS='\t' 'NR==FNR{a[$2FS$3FS$4]=$1;next}$2FS$3FS$4 in a{print $0,a[$2FS$3FS$4]}' $file ./homer/annotatePeaks/*.tsv > "common_"$(basename "$file")".xls"
	
done  

