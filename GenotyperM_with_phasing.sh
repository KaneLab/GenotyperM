#This program turns a vcf table into a fasta file given a reference sequence
#Inputs are
#	$1: multi fasta containing all references unwrapped
#	$2: vcf output of an mpileup for one or multiple organisms versus a reference
#	$3: relative or absolute filepath to the sorted.bam files.  At present, the files must be in the format <sampleID>.sorted.bam

#get array of reference genes and length of array
refs=(`grep ">" $1 | tr -d ">"`)
refLen=${#refs[@]}
echo "Building genotypes for $refLen sequences, ${refs[*]}."

#Filter the .vcf
echo "Formatting input vcf."
grep -v "#" $2 | awk '$6 == 999 && $4 != "N" && $5 !~ ","' | grep -v "INDEL" > workingVcf
echo "Formatting completed successfully."

#Get the header info and format it to work with RefSeq2
sampleArray=(`grep -v "##" $2 | head -n 1 | sed 's/ /\t/g' | sed 's/.sorted.bam//g' | cut -f10-`)
sampleLen=${#sampleArray[@]}
grep -v "##" $2 | head -n 1 | sed 's/ /\t/g' | cut -f 10- > header
sed -i 's/.sorted.bam//g' header
bash transpose header > sampleIDs; rm header
echo "Each sequence has $sampleLen individuals."


#Run bam parser for each sorted BAM in directory.
mkdir tmp_input
perl run_bam_parser.pl sampleIDs $3 workingVcf tmp_input
cat tmp_input/*.count.txt > count.txt
cat tmp_input/*.jump1.txt > jump1.txt
cat tmp_input/*.jump2.txt > jump2.txt
cat tmp_input/*.read.txt > read.txt
grep -v "#" $2 | awk '$6 == 999 && $4 != "N" && $5 !~ ","' | grep -v "INDEL" | cut -f 2,4,5 > sites.txt


#Run hapseq2 phasing program
hapseq2 --readCounts count.txt --polymorphicSites sites.txt --readHap jump1.txt --seqReadFile read.txt --seed 1 --readForwOpt 1 --mcmcHap 2 --mcmcScale 05 -o HapSeqOutput --seqError 0.01 --phase --geno --quality -r 20
rm count.txt; rm jump1.txt; rm jump2.txt; rm read.txt


#####Now reformat the output from hapseq2 to input into GenotyperM#####

echo "New stuff 1"
#Introduce delimitations into haplotypes so that we can transpose the rows into columns
awk '{print $3}' HapSeqOutput | sed 's/[a-z]/& /g' > HaplotypesSpaced
bash transpose HaplotypesSpaced > HaplotypesSpaced_T
rm HaplotypesSpaced

#Tranpose name list into columns
#awk '{print $2}' HapSeqOutput > nameList
#bash transpose nameList > nameList_T
#rm nameList

echo "New stuff 2"
#Double each sample identifier, since each sample takes up 2 columns (1 column per chromosome) and transpose
sed 'h;:a;s/[^\n]\+/&/2;t;G;ba' sampleIDs > sampleIDs_D
bash transpose sampleIDs_D > sampleIDs_D_T
#rm sampleIDs_D

echo "New stuff 3"
#Concatenate name list and haplotype info
cat sampleIDs_D_T HaplotypesSpaced_T > preGenotyperFormat
#rm sampleIDs_D_T; rm HaplotypesSpaced_T

#Get positions list
grep -v "##" $2 | cut -f 1,2 > posList

echo "New stuff 4"
#Get list of positions of variant sites from original .vcf file
cat posList preGenotyperFormat > preGenotyperFormat_S
#rm preGenotyperFormat; rm posList

echo "Last new stuff"
#Finally, we get the correct format for GenotyperM
paste preGenotyperFormat_S HaplotypesSpaced_T > workingVcf
#rm preGenotyperFormat_S; rm HaplotypesSpaced_T

#reference loop
for(( r=0; r<$refLen; r++ )); do

        #store current gene in working name
        cRef=${refs[$r]}

        #make new multifastas to write genotypes to
        > ${cRef}_L.fa
	> ${cRef}_R.fa

        #get reference in the form to genotype (supports single fasta reference for now)
        grep -A1 $cRef $1 | tail -n1 > tempRef.fa


        #only take rows that include that referece header from the fasta input
        grep $cRef workingVcf > preSnp

                #get positions of snps and snps into array and get lengths
                posArray=(`awk '{print $2}' preSnp`)
                posLength=${#posArray[@]}

		#get list of alternate nucleotides to replace the reference with
		altArray=(`awk '{print $3}' preSnp`)
		
		#make a genotyped sequence for LEFT and RIGHT columns of each sample
		for((s=0; s<$sampleLen; s++));do
			echo "Current iteration is $s out of $sampleLen"
			#get templates for overwriting with unique genotypes
		        typeL=`cat tempRef.fa`
		        typeR=`cat tempRef.fa`
	

			#get current sample
			currSample=${sampleArray[$s]}
			echo "currSample is $currSample"
			#echo "currsamp is $currSample"
			#find column to work on
			currColumnL=$((2*$s + 3))
			currColumnR=$((2*$s + 4))
			#echo "$currColumnL and $currColumnR"
                	snpArrayL=(`awk -v currColumn=$currColumnL '{print $currColumn}' preSnp`) #get LEFT snp column in vcf for current sample.
			snpArrayR=(`awk -v currColumn=$currColumnR '{print $currColumn}' preSnp`) #get RIGHT snp column in vcf for current sample.
			#echo "${snpArrayL[*]} and ${snpArrayR[*]}"
	                #iteratively sub reference nucleotide for current snp at the current position
                	for((i=0;i<$posLength;i++));do
                        	currentPos=${posArray[$i]}
	                        currentSNPL=${snpArrayL[$i]}
				currentSNPR=${snpArrayR[$i]}
				#echo "$currentSNPL $currentSNPR"
				#replace 
	                        	typeL=`echo $typeL | sed "s/./$currentSNPL/$currentPos"`
                                        typeR=`echo $typeR | sed "s/./$currentSNPR/$currentPos"`

	                        #end loop on SNP positions
	        	        done

        	        #Create new file for the genotyped individual, use $3 command line arg (prefix) for header and file name
	                #vcfName=`echo $cVcf | awk -F "\." '{print $1}'`
	                echo ">$currSample" >> ${cRef}_L.fa; echo ">$currSample" >> ${cRef}_R.fa;
	                echo $typeL >> ${cRef}_L.fa
			echo $typeR >> ${cRef}_R.fa
		done

                rm tempRef.fa
		rm preSnp
done
