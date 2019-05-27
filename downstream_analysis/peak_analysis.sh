#!/bin/bash -l 

 # module load bioinfo-tools
 # module load BEDTools/2.27.1

peakPath=$1

# Run HOMER for motif analysis and peak annotation 
for peakFile in $peakPath/*.bed; do 
	motifPath="./homer/motifs/$(basename $peakFile | grep -oP '(?<=_).*.?(?=\.)')"
	annoPath="./homer/annotatePeaks/annotatePeaks_$(basename $peakFile | grep -oP '(?<=_).*.?(?=\.)')"

	findMotifsGenome.pl $peakFile mm9 -len 12,16,18,21 -size 50 -N 20000 -mis 4 $motifPath
	annotatePeaks.pl $peakFile mm9 -size 1000 -go ./G4_go > $annoPath.tsv
done  

for file in ./homer/annotatePeaks/*.tsv; do 
 
	awk -F '\t' -v OFS='\t' 'NR==FNR{a[$2FS$3FS$4]=$1;next}$2FS$3FS$4 in a{print $0,a[$2FS$3FS$4]}' $file ./homer/annotatePeaks/*.tsv > "common_"$(basename "$file")".xls"
	
done  

# Overlap with GQP and G4Seq
bedFile="G4Seq.mm9.bed"
for f in ./bed_files/*.bam_peaks.bed; do
	echo "Creating venn diagram for G4 prediction file and $(basename "$file")"
	intervene venn -i ./bed_files/$bedFile $f -o ./venn_diagrams 
done 

# Create BED file with scores for annotation 
awk -F '\t' '$8 ~ "TSS" { print } ' annotatePeaks_G4_06_ALL.tsv | cut -f 2-5 > 06_TSS.bed

awk -F '\t' -v OFS='\t' '$8 == "Annotation"{$7 = $7 FS "Score"; print ; } $8 ~ "TSS"{ $7 = $7 FS "1";  print ; } $8 ~ "Intergenic"{ $7 = $7 FS "2"; print ; } $8 ~ "intron"{ $7 = $7 FS "3"; print ; } $8 ~ "3'"'"' UTR"{ $7 = $7 FS "4"; print ; } $8 ~ "non.coding"{ $7 = $7 FS "5"; print ; } $8 ~ "5'"'"' UTR"{ $7 = $7 FS "6"; print ; } $8 ~ /exon \(NM_.*\)/{ $7 = $7 FS "7"; print ; } $8 == "NA"{$7 = $7 FS "0"; print ; } $8 ~ "TTS"{$7 = $7 FS "TTS"; print } ' annotatePeaks_G4_06_ALL.tsv | cut -f 2-5,7 > 06_ALL.bed

