#!/bin/bash

#this script thanks to a design file, will give for each sample of each condition, the counts with kallisto.
#the results will be in the output directory that you have set
#each sample has its own sub-directory, and in each, the counts are in *_abundance.tsv


############ input variables #############
#################################

#the design should be something like this (tab-separated) :
# sample_name		condition
# A	cond1
# B	cond1
# O	cond2
# P	cond2
# ...
#design="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/kallisto_transcripts_counts_sample.tsv"

#reads path
#reads_path="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/"

#genome
#genome_fasta="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_riboseq_DIS3_diagenode_merged_runs_15_11_2021_and_16_12_2021/Dominika_ribotricer_analysis_gencode_scallop_eRNA_selected_length_from_unique_25_34nt_extendedStartCodon/ORFs_in_fasta_format/gencode26_scallop_eRNA_ORFs_extendedStartCodon.fa"

#output files
#counter_outputs="/media/marcgabriel/saylar17/Dominika_GTEX_kallisto_counts_ribotricer_gencode26_DIS3_scallop_HCT116_eRNAs_conservative_refined_candidate_orfs/"

#not used anymore (unless pseudomapping is "on")
chrom_size="/home/marcgabriel/Documents/gencode32/chrom_size.tsv"
input_gff="/home/marcgabriel/Documents/gencode32/gencode.v32.annotation_sorted.gff3"

#e.g : pre_R_1_, post_PD_1_
#separator_after_file_name="_"
separator_after_file_name="_"

#separator_after_file_name="."

#number of process (samples) to run in parallel
process_number_limit=8

#threads used by kallisto for one sample
kallisto_threads=4

#kallisto program
#kallisto="/home/marcgabriel/new_dekupl/dekupl-run/bin/kallisto"
kallisto="/home/marcgabriel/Downloads/kallisto/kallisto"

gffread="/home/marcgabriel/Desktop/gffread-0.9.12.Linux_x86_64/gffread"

do_bam="no"

#kmer size for kallisto index (default=31)
#kmer_size=15
kmer_size=31

#library type ; if it's unstranded, just put "", otherwise choose between --rf-stranded & --fr-stranded
#--rf-stranded = mate2 is in the same direction as the RNA
#library_type="--fr-stranded"
library_type=""
#library_type="--rf-stranded"

#clean the existing bam file, if you think it's not in a standard way ("0" pos fort start or end, no string in the CIGAR == Ving gives eerros becaause of that)
clean_existing_bam="no"

check_file="no"
#check_file="yes"

#check_reads="yes"
check_reads="no"

#fragment length (for single-end)
fragment_length=100

sd_of_frag_len=80

###############  process the inputs ###################
#######################################################

#standardize the output dir name
counter_outputs="${counter_outputs}/"
counter_outputs=$(echo $counter_outputs |sed 's/\/\//\//g')
if [ ! -d $counter_outputs ]; then mkdir -p $counter_outputs; fi 

#directory of scripts for each sample (to be laucnhed in parallel)
sample_scripts_dir="${counter_outputs}/sample_scripts/"
if [ -d $sample_scripts_dir ]; then rm -rf $sample_scripts_dir; fi 
mkdir $sample_scripts_dir


######### dev ###########

#if no kallisto gtf, construct it
myGtf=${counter_outputs}kallisto_gtf.gtf

if [[ ! -f ${counter_outputs}kallisto_gtf.gtf ]];then

        #take only exons to facilitate the downstream processing
		grep -P "\texon\t" $input_gff >${counter_outputs}raw_gtf.gtf

		#gff to gtf
		$gffread ${counter_outputs}raw_gtf.gtf --force-exons -T -F -o ${counter_outputs}raw_gtf.tmp && mv ${counter_outputs}raw_gtf.tmp ${counter_outputs}raw_gtf.gtf

		cat ${counter_outputs}raw_gtf.gtf| sed 's/[;]/\t/g'|cut -f1-11 |sort -k1,1 -k10,10 -k11,11 -k7,7 | bedtools groupby -g 1,10,11,7 -c 4,5 -o min,max|awk -F"\t" 'OFS="\t"{print $1,".","gene",$5,$6,".",$4,".",$2";"" "$3";"}'|sed 's/ gene_id/gene_id/g' >${counter_outputs}raw_gtf.tmp && mv ${counter_outputs}raw_gtf.tmp ${counter_outputs}raw_gtf_genes.gtf 

		cat ${counter_outputs}raw_gtf.gtf| sed 's/[;]/\t/g' |cut -f1-11|sort -k1,1 -k10,10 -k11,11 -k9,9 -k7,7 | bedtools groupby -g 1,10,11,9,7 -c 4,5 -o min,max|awk -F"\t" 'OFS="\t"{print $1,".","transcript",$6,$7,".",$5,".",$2";"" "$3";"" "$4";"}' |sed 's/ gene_id/gene_id/g'>${counter_outputs}raw_gtf.tmp && mv ${counter_outputs}raw_gtf.tmp ${counter_outputs}raw_gtf_transcripts.gtf 

		cat ${counter_outputs}raw_gtf_genes.gtf ${counter_outputs}raw_gtf_transcripts.gtf ${counter_outputs}raw_gtf.gtf |sort -k1,1 -k4,4n -k4,5nr >${counter_outputs}raw_gtf.tmp && mv ${counter_outputs}raw_gtf.tmp ${counter_outputs}raw_gtf.gtf


		awk '( $3 ~ /gene/ )' ${counter_outputs}raw_gtf.gtf > $myGtf
		awk '( $3 ~ /transcript/ )' ${counter_outputs}raw_gtf.gtf >> $myGtf
		awk '( $3 ~ /exon/ && $7 ~ /+/ )' ${counter_outputs}raw_gtf.gtf | sort -k1,1 -k4,4n >> $myGtf
		awk '( $3 ~ /exon/ && $7 ~ /-/ )' ${counter_outputs}raw_gtf.gtf | sort -k1,1 -k4,4nr >> $myGtf


fi



##########"end of dev #######"

#look for all conds ; keep the order in the file
all_conds=($(grep -v "^#" $design|cut -f2 |awk '!seen[$0]++'))

echo -e "all conds are : ${all_conds[*]}\n"

#array to store all files
all_files=()

#array with all new names
all_new_names=()

#set -x
#for each sample of each cond, assign a name based on the condition (for conds normal & tumoral with 2 samples each : normal_1, normal_2...tumoral_1, tumoral_2)
for one_cond in $(seq 0 $((${#all_conds[*]}-1)));do

	#for a given cond, take all its files (samples)
	#files_one_cond=($(grep -P "\t${all_conds[$one_cond]}$" $design |cut -f1|sort -n))
	files_one_cond=($(grep -v "^#" $design |grep -P "\t${all_conds[$one_cond]}$" |cut -f1))
 
 
	#number of the 1st rep of the condition
	cond_start=1
 
	#number of the last rep of the condition
	cond_end=${#files_one_cond[*]}
 
	#store all the new names for this condition
	cond_new_names=()
 
    #loop across the samples of the cond, the for each of them concatenate the name of the condition & the number
	for ((i=$cond_start; i<=$cond_end; i++)); do cond_new_names+=(${all_conds[$one_cond]}_${i}) ; done
 
	all_files+=(${files_one_cond[*]})
 
	all_new_names+=(${cond_new_names[*]})

done

#echo -e "all files are : ${all_files[*]}"

#echo -e "new names are : ${all_new_names[*]}"

echo -e "new design : "

paste <(echo -e ${all_files[*]}|sed 's/ /\n/g') <(echo -e ${all_new_names[*]} |sed 's/ /\n/g')

paste <(echo -e ${all_files[*]}|sed 's/ /\n/g') <(echo -e ${all_new_names[*]} |sed 's/ /\n/g') >${counter_outputs}new_design.tsv



echo -e "sample\tnb_reads_R1" >${counter_outputs}reads_summary.tsv



index_name=$(basename $genome_fasta|sed 's/\.fa//g'|sed 's/\.fasta//g')

#if the kallisto index of the genome isn't built, do it
index="${counter_outputs}${index_name}_kallisto_index.idx"

#build the index file if it's not there
if [ ! -f $index ];then 

   echo -e "\nwe are going to make the kallisto index with $genome_fasta\n..."

  $kallisto index -i $index --kmer-size=$kmer_size $genome_fasta

fi


#summary of the count files
>${counter_outputs}summary.tsv


##############  run a loop across the samples in order to process them #################
########################################################################################

all_results=()


#create a script for each sample, that will run Kallisto & process the output file
for one_sample in $(seq 0 $((${#all_files[*]}-1)));do

          #for a given sample of a given cond, take its new name
		  one_sample_name=$(echo "${all_new_names[$one_sample]}")
		  
		  #if [[ "$one_sample_name" != "Colon_49" ]];then
		  
		  #continue
		  
		  #fi
		  
		  library_type=""
		  
		  if [[ "$one_sample_name" =~ "LNCaP_FV" ]];then
		  
		    library_type="--rf-stranded"
		    
		    echo -e "we're going to change the lib type for LNcap Rocco...\n"
		  
		  fi
		  
		  #path to the sample
		  counter_outputs_OneSample="${counter_outputs}${one_sample_name}/"
		  
		  if [ ! -d $counter_outputs_OneSample ]; then mkdir $counter_outputs_OneSample; fi 
		  
		   echo -e "#!/bin/bash\n\n" >${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			#if the count file isn't there, use kallisto 
		   #echo -e "\nif [ ! -f ${counter_outputs_OneSample}abundance.tsv ] || [[ \$(wc -l ${counter_outputs_OneSample}abundance.tsv|awk '{print \$1}') -lt 1 ]];then\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
		   
		   echo -e "\nif [ ! -f ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv ] || [[ \$(head ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv|wc -l|awk '{print \$1}') -lt 1 ]];then\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
		   
		   
		   
		   #set -x
		  		 
			  #fastq_list_OneSample=($(find $reads_path -name "${all_files[$one_sample]}${separator_after_file_name}*" |grep ".*\.fastq"|sort -n)) || { echo "no fastq files 1 !!" 1>&2; exit; }
			  
			  fastq_list_OneSample=($(find $reads_path -name "${all_files[$one_sample]}${separator_after_file_name}*"|grep -v -i "report" |grep ".*\.fastq"|sort -n)) 
			  #|| { echo "no fastq files 1 !!" 1>&2; exit; }
			  
			  if [[ ${#fastq_list_OneSample[*]} -eq 0 ]];then
			  
			     echo -e "\n-> no fastq files with sep \"${separator_after_file_name}\" for $one_sample_name ($one_sample) , we're going to try with \".\"..."
			  
			     fastq_list_OneSample=($(find $reads_path -name "${all_files[$one_sample]}.*"|grep -v -i "report" |grep ".*\.fastq"|sort -n))|| { echo "no fastq files 2 !!" 1>&2; exit; }
			  
			  
			  fi
			  
			  #set +x
			  
			  echo -e "echo -e \"fastq files for condition $one_sample_name : \n ${fastq_list_OneSample[*]} \n--------------\n\"" >>${sample_scripts_dir}subscript_${one_sample_name}.sh

			  
			  #if [[ ${#fastq_list_OneSample[*]} -eq 2 ]];then
			  
			 
			  
				  #fastq_list_OneSample_R1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "1\.fastq"|sort -u)
				  
				  
				  #set -x
				  #fastq_list_OneSample_R1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "1_.*\.fastq"|sort -u|head -n1)
				  
				  fastq_list_OneSample_R1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "_R1.*\.fastq"|sort -u|head -n1)
				  
				  if [[ ! -f "$fastq_list_OneSample_R1" ]];then
				  
				    fastq_list_OneSample_R1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "1\.fastq"|sort -u|head -n1)
				  
				  fi
				  
				   
				 
				  #fastq_list_OneSample_R2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "2\.fastq"|sort -u)
				  
				  #fastq_list_OneSample_R2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "2_.*\.fastq"|sort -u|head -n1)
				  
				  fastq_list_OneSample_R2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "_R2.*\.fastq"|sort -u|head -n1)
				  
				  
				  if [[ ! -f "$fastq_list_OneSample_R2" ]];then
				  
				    fastq_list_OneSample_R2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "2\.fastq"|sort -u|head -n1)
				  
				  fi
				  #set +x
				  
				  if [[ "$check_file" == "yes" ]];then
				  
						if [[ $(gzip -t $fastq_list_OneSample_R1 2>&1|grep "not in") ]] || [[ $(gzip -t $fastq_list_OneSample_R1 2>&1|grep "unexpected end of") ]] ;then echo -e "-> issue with $fastq_list_OneSample_R1 !\n";fi
				  
				  fi
				  
				  if [[ "$check_reads" == "yes" ]];then
				  
						nb_reads=$(zcat $fastq_list_OneSample_R1|wc -l |awk '{OFMT="%f";print $1/4}' )
						
						echo -e "${one_sample_name}\t${nb_reads}" >>${counter_outputs}reads_summary.tsv
				  
				  fi
				  
	      
			      if [[ "$do_bam" == "yes" ]];then
			      
				  #run the quantif
			      echo -e "${kallisto} quant -i $index $library_type -t $kallisto_threads --chromosomes=$chrom_size --genomebam --gtf=$myGtf -o ${counter_outputs_OneSample} ${fastq_list_OneSample_R1} ${fastq_list_OneSample_R2}|| { echo \"kallisto failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			      
			      #remove reads without start & end positions, or CIGAR string
			      echo -e "mv ${counter_outputs_OneSample}/pseudoalignments.bam ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam && samtools view -@ $kallisto_threads -h ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam|awk 'OFS=\"\\\t\"{if((\$4!=0 && \$8 !=0 && \$6!=\"*\")||\$1~/^@/){print}}'|samtools view -@ $kallisto_threads -Sbh - >${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam && rm ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam && samtools index ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam || { echo \"kallisto bam failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh			      
			      
			      else
			      
			      
				  #run the quantif
			      echo -e "${kallisto} quant -i $index $library_type -t $kallisto_threads -o ${counter_outputs_OneSample} ${fastq_list_OneSample_R1} ${fastq_list_OneSample_R2}|| { echo \"kallisto failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			      
			      			      
			      
			      
			      
			      fi
				  
			  #elif [[ ${#fastq_list_OneSample[*]} -eq 1 ]];then
			  
			  
				  #if [[ "$check_file" == "yes" ]];then
				  
						#if [[ $(gzip -t ${fastq_list_OneSample[*]} 2>&1|grep "not in") ]];then echo -e "-> issue with $fastq_list_OneSample_R1 !\n";fi
				  
				  #fi			     
			  

			       
			       #if [[ "$do_bam" == "yes" ]];then
			       
			       #echo -e "${kallisto} quant -i $index $library_type -t $kallisto_threads --chromosomes=$chrom_size --genomebam --gtf=$myGtf --single --fragment-length=$fragment_length --sd=$sd_of_frag_len -o ${counter_outputs_OneSample} ${fastq_list_OneSample[*]}|| { echo \"kallisto failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			       
			       #echo -e "mv ${counter_outputs_OneSample}/pseudoalignments.bam ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam && samtools view -@ $kallisto_threads -h ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam|awk 'OFS=\"\\\t\"{if((\$4!=0 && \$8 !=0 && \$6!=\"*\")||\$1~/^@/){print}}'|samtools view -@ $kallisto_threads -Sbh - >${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam && rm ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam && samtools index ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam || { echo \"kallisto bam failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			       
			       
			       #else
			       
			       
			       #echo -e "${kallisto} quant -i $index $library_type -t $kallisto_threads --single --fragment-length=$fragment_length --sd=$sd_of_frag_len -o ${counter_outputs_OneSample} ${fastq_list_OneSample[*]}|| { echo \"kallisto failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh			       
			       
			       
			       #fi
			   
			  
			  #elif [[ ${#fastq_list_OneSample[*]} -eq 0 ]];then
			  
			     #echo -e "no fastq file for ${all_files[$one_sample]} $one_sample_name"
			     
			     ##exit 1
			     
			     #echo -e "\nfi\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			     
			     ##it could be usefull if you just want to process some files in the list
			     #continue
			  
			  
			  #elif [[ ${#fastq_list_OneSample[*]} -gt 2 ]];then
			  
			    #echo -e "too many fastq files for ${all_files[$one_sample]} / $one_sample_name : ${fastq_list_OneSample[*]}"
			    
			    #exit 1
			  
			  #fi
			  
			  echo -e "cp ${counter_outputs_OneSample}abundance.tsv ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  echo -e "\nfi\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh

			  
			  
			  
			  #process the count file in order to have columns like this : ID	counts
			  #echo -e "awk 'OFS=\"\\\t\"{if(NR>1){print \$1,\$4}}' ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv >${counter_outputs_OneSample}${one_sample_name}_abundance.tmp && mv ${counter_outputs_OneSample}${one_sample_name}_abundance.tmp ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  
			  echo -e "if [[ \$(head -n1 ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv|awk '{print NF}') -gt 2 ]];then\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  			  
				echo -e "awk 'OFS=\"\\\t\"{if(NR>1){print \$1,\$4}}' ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv >${counter_outputs_OneSample}${one_sample_name}_abundance.tmp && mv ${counter_outputs_OneSample}${one_sample_name}_abundance.tmp ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  echo -e "fi\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  all_results+=(${counter_outputs_OneSample}${one_sample_name}_abundance.tsv)
			  
			  
				
			  #give rights to the script	
			  chmod 755 ${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  echo -e "${counter_outputs_OneSample}${one_sample_name}_abundance.tsv\t${one_sample_name}" >>${counter_outputs}summary.tsv
			  
			  #ensure that the bam is well formed
			  if [[ "$do_bam" == "yes" ]];then
			  
			    #force to clean the existing bam file
			    if [[ "$clean_existing_bam" == "yes" ]];then
			  
			  			       echo -e "samtools view -@ $kallisto_threads -h ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam|awk 'OFS=\"\\\t\"{if((\$4!=0 && \$8 !=0 && \$6!=\"*\")||\$1~/^@/){print}}'|samtools view -@ $kallisto_threads -Sbh - >${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.test.bam && mv ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.test.bam ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam && samtools index ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam || { echo \"kallisto bam failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  			       
			    fi
			  
			  fi
					   
done



#run the subscripts in parallel

#here, if you are on a cluster, and you want to use a qsub command, just use a loop on the content of ${sample_scripts_dir}, and give the result to qsub
find ${sample_scripts_dir} -name "*subscript_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "executing bash sorting and/or counting subsritpts failure !" 1>&2; exit; }


all_IDs=($(cut -f2 ${counter_outputs}new_design.tsv))

echo -e "all samples are : ${all_IDs[*]}"





header=$(echo -e "${all_IDs[*]}"|sed 's/ /\t/g')

#final combined counts
echo -e "feature\t${header}" >${counter_outputs}header.tsv


#initialize the pos of the count files
#a=1

#header_table=(feature)

#for i in ${all_IDs[*]};do


    #one_file=$(find ${counter_outputs} -name "${i}_abundance.tsv"|grep -E "\/${i}_abundance.tsv")
    
    #echo -e "sample is : ${i}, file is $one_file"
    
    #if [[ -f $one_file ]];then


		#if [[ $a -eq 1 ]];then

			#cat $one_file >${counter_outputs}combined_counts.tsv
			
		#else
		
			#paste ${counter_outputs}combined_counts.tsv <(cut -f2 $one_file) >${counter_outputs}combined_counts.tmp && mv ${counter_outputs}combined_counts.tmp ${counter_outputs}combined_counts.tsv
		
		
		#fi
		
		#header_table+=(${i})
	  
	  
	  #a=$((a+1))
  
  
  #else
  
    #echo -e "${i}_abundance.tsv doesn't exist, next !\n"
  
  #fi

#done

#with too many files, you have to increase the limits
ulimit -Sn 10240 # The soft limit


#for a given line, from the 2nd column, if the column number has 0 after modulo 2, store the result (the count is at the columns with even numbers)
cat ${counter_outputs}header.tsv <(paste ${all_results[*]} | awk '{all=$1;for(i=2;i<=NF;i++){if(i%2==0){all=all"\t"$i}};print all}') >${counter_outputs}combined_counts.tsv



#cat ${counter_outputs}header.tsv ${counter_outputs}combined_counts.tsv >${counter_outputs}combined_counts.tmp && mv ${counter_outputs}combined_counts.tmp ${counter_outputs}combined_counts.tsv

#cat <(echo -e "${header_table[*]}"|sed 's/ /\t/g') ${counter_outputs}combined_counts.tsv >${counter_outputs}combined_counts.tmp && mv ${counter_outputs}combined_counts.tmp ${counter_outputs}combined_counts.tsv

exit


if [[ $(grep "\.gz" $genome_fasta|wc -l) -eq 0 ]];then


	cat <(echo -e "id\tlength") <(awk 'BEGIN{RS=">"}''{if(NR>1){sub("\n","\t"); gsub("\n",""); print $1"\t"length($2)}}' $genome_fasta) >${counter_outputs}matching_id_length.tsv

else

	cat <(echo -e "id\tlength") <(zcat $genome_fasta|awk 'BEGIN{RS=">"}''{if(NR>1){sub("\n","\t"); gsub("\n",""); print $1"\t"length($2)}}') >${counter_outputs}matching_id_length.tsv
	
fi

echo -e "\n== check file : ${counter_outputs}combined_counts.tsv ==\n"






