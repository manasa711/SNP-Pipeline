#!/bin/bash

# Function for doing your getopts
get_input () {

	while getopts 'a:b:r:eo:f:zvih' OPTION;
	do
		case "$OPTION" in
			a) reads1=$OPTARG;; #takes in the file: pair1
			b) reads2=$OPTARG;; #takes in the file: pair2
			r) ref=$OPTARG;; #for the reference genome
			e) realign=1;; # re-alignment
			o) output=$OPTARG;; #output file name
			f) millsFile=$OPTARG;;  #mills file location
			z) gunzip=1 ;; #used for output file to be gunzipped
			v) v=1;; #verbose mode
			i) index=1;; #providing a variable for indexing the output BAM file
			h) h=1 ; echo "Display help" ;;
			*) echo "Please assign the defined arguments. Select -h for help"
				exit 1;;
		esac
	done

	if (( h ))
	then
		echo " The bash script takes in the following arguments:
			-a : Reads input for fasta file: pair 1
			-b : Reads input for fasta file: pair 2
			-r : Input for the refernce genome file
			-e : for performing read re-alignment
			-o : for specifying the output VCF file name
			-f : Specifies the millsFile location
			-z : Output VCF file will be gunzipped
			-v : Verbose mode; let's you know what the script is doing
			-i : index the output BAM file
			-h : prints usage information

			The bash script needs to be run in the following format while also specifying the arguments given above : ./snp_pipeline.bash "
		exit
	fi

}

# Function for checking for presence of input files, reference genome,
# and the output VCF file
#
# Input: File locations (string)
# Output: True, if checks pass; False, if checks fail (bool)

check_files () {

	if (( v ))
	then
		echo "Checking for the presence of input files, reference genome and output VCF file"
	fi

	if [ -f "$reads1" ]  #checks if the first FASTQ file exists
	then
		echo "$reads1 exists"
	else
		echo "$reads1 is missing"
		exit 1
	fi

	if [ -f "$reads2" ] #checks if the second FASTQ file exists
	then
		echo "$reads2 exists"
	else
		echo "$reads2 is missing"
		exit 1
	fi

	if [ -f "$ref" ] #checks if the reference genome exists
	then
		echo "$ref exists"
	else
		echo "$ref is missing"
		exit 1
	fi

	if [ -f "$millsFile" ] #checks if the indels millsFile exists
	then
		echo "$millsFile exists"
	else
		echo "$millsFile is missing"
		exit 1
	fi

	if [ -f "$output" ] #checks for the presence of the output file
	then
		echo "$output exists"
		echo "If you want to overwrite the existing output file, please enter '1' or if you wish to exit the bash script, enter '0'"
		read -r answer
		if [ "$answer" == '0' ]
		then
			exit
		fi
	else
		echo "$output will be created at the end of the program"

	fi

}

	# Preparing the tempary directory

prepare_temp () {

	if (( v ))
	then
		echo "Preparing the tempary directory and indexing the reference genome for mapping"
	fi

	if [ -f ${ref}.bwt ]
	then
		echo "Indexed file exists."
	else

		bwa index "$ref" #indexing the reference genome for mapping
	fi
}

# Function for the mapping step of the SNP-calling pipeline
# Input: File locations (string), Verbose flag (bool)
# Output: File locations (string)

mapping () {

	if (( v ))
	then
		echo "The mapping step of the SNP-calling pipeline is running"
	fi

	if [ -f lane.sam ]
	then
		echo "Exists."
	else

		bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' "$ref" "$reads1" "$reads2" > lane.sam #mapping the reads
	fi

	if [ -f lane_fixmate.bam ]
	then
		echo "E."
	else
		samtools fixmate -O bam lane.sam lane_fixmate.bam #cleaning read pair information and flags
	fi

	if [ -f lane_sorted.bam ]
	then
		echo "E."
	else
		samtools sort -O bam -o	lane_sorted.bam -T lane_temp lane_fixmate.bam #sorting them according to their coordinates
	fi

	echo "Mapping"
}

# Function for improving the number of miscalls
# Input: File locations (string)
# Output: File locations (string)

improvement () {

	if (( v ))
	then
		echo "Preparing the fasta file to be used as reference by creating a dictionary and a fasta index file"
	fi

	#Preparing the fasta file to be used reference

	if [ -f chr17.dict ]
	then
		echo "e."
	else
		samtools dict "$ref" -o chr17.dict  #creating the fasta sequence dictionary file
	fi

	if [ -f $ref.fai ]
	then
		echo "e."
	else
		samtools faidx "$ref" "$ref".fai #creating the fasta index file
	fi

	if (( index ))
        then
                if (( v ))
                then
                        echo "Indexing the BAM file"
                fi

                        samtools index lane_sorted.bam lane_sorted.bam.bai
        fi


	if (( realign ))
	then
	#realigning raw gapped alignment with Broad's GATK Realigner to reduce the number of miscalls of INDELs
		if (( v ))
		then
			echo "Realigning raw gapped alignment"
		fi

		java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R "$ref" -I lane_sorted.bam -o lane.intervals --known "$millsFile" --log_to_file "SNP_Pipeline_1.log"
	  java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R "$ref" -I lane_sorted.bam -targetIntervals lane.intervals -known "$millsFile" -o lane_realigned.bam --log_to_file "SNP_Pipeline_2.log"
	fi

	#Using the mark duplicates tool to compile all of the reads from each librbary into one BAM file and also marks PCR and optical duplicates

	#java -Xmx2g -jar MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT INPUT=lane_1.bam INPUT=lane_2.bam INUT=lane_3.bam OUTPUT=library.bam

	#merge step

	#samtools merge sample.bam library1.bam library2.bam library3.bam
	#samtools index sample.bam

	echo "Improvement"
}

# Function to call variants
# Input: File locations (string)
# Ouput: None

call_variants () {

	if (( v ))
	then
		echo "The function to call variants is running"
	fi

	if (( gunzip ))
	then
		if (( realign ))
		then
			bcftools mpileup -Ou -f "$ref" lane_realigned.bam | bcftools call -vmO z -o "$output".vcf.gz
		else
			bcftools mpileup -Ou -f "$ref" lane_sorted.bam | bcftools call -vmO z -o "$output".vcf.gz
		fi
	else
		if (( realign ))
		then
			bcftools mpileup -Ou -f "$ref" lane_realigned.bam | bcftools call -vmO v -o "$output".vcf
		else
			bcftools mpileup -Ou -f "$ref" lane_sorted.bam | bcftools call -vmO v -o "$output".vcf
		fi
	fi
	#using mpileup to produce a BCF file that contains all the locations of the genome inorder to convert the BAM file into genomic positions

	tabix -p vcf "$output".vcf.gz #indexing the VCF

	bcftools filter -O z -o study_filtered..vcf.gz -s LOWQUAL -i'%QUAL>10' "$output".vcf.g #filtering the data

	echo "Calling Variants"
}


main() {

	get_input "$@"
	check_files
	prepare_temp
	mapping
	improvement
	call_variants
}

# Calling the main function
main "$@"



bats_test (){
    command -v bats
}
