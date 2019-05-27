 #!/bin/bash

# module load bio-info tools
# module load MACS/2.1.2

# Input path to BAM folder and output directory for peak files 
bamPath=$1
MACS2_out_dir="./peaks" # Specify output directory

set -- "$bamPath"/*.bam

# Handle the case where no files were found matching our glob
[[ -e $1 || -L $1 ]] || { echo "No bam files found!" >&2; exit 1; }

# IPtarget 
IPtarget="G4"
inputName="IN"

# If same input file for all files please specifiy glob as X in "IN_X.bam"
glob="ALL.mm9"

# Peak calling parameters
gen_size="1.87e9" #mm9
pvalue="1e-5"

i=0
pairs=()

for bamFile in $bamPath/*.bam; do 
	type=$(echo $(basename "$bamFile") | cut -d _ -f 1)
	name=$(echo $(basename "$bamFile") | grep -oP '(?<=_).*.?(?=\.)') 

	if [ ${#pairs[@]} -eq 0 ]; then 		  
		pairs[i++]=$(find $bamPath/*.bam -type f \( -name "$IPtarget"_"$name.bam" -o -name "$inputName"_"$name"."bam" \) | sort -u | tr '\n' ' ' ):$(find $bamPath/*.bam -type f \( -name "$IPtarget"_"$name.bam" -o -name "$inputName"_"$name"."bam" \) | wc -l )
	else
    	if ! [[ "${pairs[@]:0}" == *"$bamFile"* ]]; then 
    		pairs[i++]=$(find $bamPath/*.bam -type f \( -name "$IPtarget"_"$name.bam" -o -name "$inputName"_"$name"."bam" \) | sort -u | tr '\n' ' ' ):$(find $bamPath/*.bam -type f \( -name "$IPtarget"_"$name.bam" -o -name "$inputName"_"$name"."bam" \) | wc -l ) 
    	fi 
	fi

	echo "---- Generating bigWig tracks for $(basename "$bamFile")"

	# Generate bigWig tracks
	# bamCoverage -b $bamFile -o ./bw/$type"_"$name.bw 
done  

for key in ${!pairs[@]}; do 
	pairNr=$(echo ${pairs[$key]} | cut -d : -f 2)
	pairFiles=$(echo ${pairs[$key]} | cut -d : -f 1)
		if [ $pairNr == 2 ]; then 
			for file in ${pairFiles[@]}; do
				if [ $(echo $(basename "$file") | cut -d _ -f 1) == $IPtarget ]; then
					chip="$file"
				else 
					input="$file"
					
					echo "---- Performing peak calling on $(basename "$chip") as treatment and $(basename "$input") as control"
					# Peak calling 
					macs2 callpeak -t $chip -c $input -g $gen_size -p $pvalue -n $(basename "$chip") --outdir $MACS2_out_dir 
					chip=""
					input=""
				fi 
			done
 		elif [ $pairNr == 1 ]; then
 				if [ $(echo $(basename "$pairFiles") | cut -d _ -f 1) == $IPtarget ]; then
					chip="$pairFiles"
					echo "---- Performing peak calling on $(basename "$chip") without control"
					# Peak calling 
					# macs2 callpeak -t $chip -g $gen_size -p $pvalue -n $(basename "$pairFiles") --outdir $MACS2_out_dir
				else 
					echo "---- Cannot perform peak calling with only input file $(basename "$pairFiles")"
				fi
 		fi
done 

# Generate bed files from MACS2 peak files

for peakFile in $MACS2_out_dir/*.narrowPeak; do 
	
	echo "---- Generating bed file and sorted bed file for $(basename "$peakFile")"
	
	cut -f 1-3 $peakFile > ./bed_files/$(basename "$peakFile" | grep -oP '(?<=_).*.?(?=\.)').bed
	# Filter and sort based on q-value 
	awk '$9>2' $peakFile | sort -t $'\t' -k9,9rn $peakFile | cut -f 1-3  >  ./bed_files/$(basename "$peakFile" | grep -oP '(?<=_).*.?(?=\.)')_sorted.bed
done  
bedtools multiIntersectBed -i ./bed_files/*_peaks.bed > consensus_peaks.bed 
bedtools multiIntersectBed -i  ./bed_files/*_sorted.bed > consensus_peaks_sorted.bed
# Assessment of reproducibility  
