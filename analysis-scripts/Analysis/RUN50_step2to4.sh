#!/bin/bash

source activate longread_umi_HIV
#These should be fixed as standard
#set working directory
wrk=/data/thesisguest1/work/UMI

run_ID=RUN50
prj=$wrk/data_"$run_ID"

sample_id="$run_ID"_guppy_v6.5.7_sup_BC

raw_data_loc=$prj/raw_fastq_guppy_v6.5.7_sup


FW2=AAGTAGTGTGTGCC
RV2=GGAAAGTCCCCAGC

MIN_LENGTH=100
MAX_LENGTH=10000

IN_DIR="$prj"/demultiplex-participants_replicates/

THREADS=20

DEBUG="YES"

OUT_DIR="$prj"/demultiplex-participants_replicates/cluster


############# UMI SAMPLE DEMULTIPLEXING #############
#demultiplex UMI bins by replicate (defined by ONT barcode) and participant (defined by PCR barcode) 

mkdir $prj/demultiplex-participants_replicates

### Demultiplex the UMI bins present in one replicate (1 ONT fastq file) by participant (by PCR barcode)

for bc in {01..06}; do 
longread_umi demultiplex -r $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq -c $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/raconx3_medakax2/consensus_raconx3_medakax2.fa -n '1-8' -o $prj/demultiplex-participants_replicates/$sample_id"$bc" -t 20 -b $CONDA_PREFIX/longread_umi/scripts/barcodes_HIV.tsv -u $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/umi_binning/read_binning/umi_bin_map.txt
done

### Indicate the PCR replicate origins by different numbers for each participant

replicate=1

for bc in {01..06}; do  
  for pcr_barcode in {1..8}; do
    awk -v pat="$replicate" '/^>/ {$1=$1";rep0"pat""} 1' $prj/demultiplex-participants_replicates/$sample_id"$bc"/barcode"$pcr_barcode".fa > $prj/demultiplex-participants_replicates/barcode"$pcr_barcode".rep0"$replicate".fa
    awk -v pat="$replicate" '/^>/ {$1=$1";rep0"pat""} 1' $prj/demultiplex-participants_replicates/$sample_id"$bc"/undetermined.fa > $prj/demultiplex-participants_replicates/undetermined.rep0"$replicate".fa
  done
  replicate=$((replicate+1))
  find $prj/demultiplex-participants_replicates/ -type f -empty -delete
done

### Combine all UMI bins from different replicate into one file per participant sample

for pcr_barcode in {1..8}; do
  cat $prj/demultiplex-participants_replicates/barcode"$pcr_barcode".rep*.fa >  $prj/demultiplex-participants_replicates/"$run_ID".barcode"$pcr_barcode".replicates.fa
  cat $prj/demultiplex-participants_replicates/undetermined.rep*.fa >  $prj/demultiplex-participants_replicates/"$run_ID".undetermined.replicates.fa
  find $prj/demultiplex-participants_replicates/ -type f -empty -delete
done 

### list the number of UMI bins passing QC for each replicate per participant
today=$(date +"%Y_%m_%d")

 rm $prj/demultiplex-participants_replicates/umi_bin_replicate_stats_$today.txt

 echo -e "Run_ID\tBarcode\tID1\tID2\tID3\tID4\tID5\tID6\tID7\tID8\tundertermined" >> $prj/demultiplex-participants_replicates/umi_bin_replicate_stats_$today.txt
 
for bc in {01..06}; do 
    cd $prj/demultiplex-participants_replicates/$sample_id"$bc";
    ID1=$(seqkit stats barcode1.fa -T | awk 'NR!=1 {print $4}');
    ID2=$(seqkit stats barcode2.fa -T | awk 'NR!=1 {print $4}');
    ID3=$(seqkit stats barcode3.fa -T | awk 'NR!=1 {print $4}');
    ID4=$(seqkit stats barcode4.fa -T | awk 'NR!=1 {print $4}');
    ID5=$(seqkit stats barcode5.fa -T | awk 'NR!=1 {print $4}');
    ID6=$(seqkit stats barcode6.fa -T | awk 'NR!=1 {print $4}');
    ID7=$(seqkit stats barcode7.fa -T | awk 'NR!=1 {print $4}');
    ID8=$(seqkit stats barcode8.fa -T | awk 'NR!=1 {print $4}');
    undet=$(seqkit stats undetermined.fa -T | awk 'NR!=1 {print $4}');
    echo -e $run_ID'\t'BC"$bc"'\t'"$ID1"'\t'$ID2'\t'$ID3'\t'$ID4'\t'$ID5'\t'$ID6'\t'$ID7'\t'$ID8'\t'$undet >> $prj/demultiplex-participants_replicates/umi_bin_replicate_stats_$today.txt;  

done

****

############# CHECK FOR UMI CROSSOVER #############
# incorrect demultiplexing of ONT long-reads by guppy basecaller can lead to a small number of reads being incorrectly assigned to another ONT barcode, this code detects crossover UMI bins (identical UMI bins found in different ONT barcode (=replicates). Only retain UMI bin for the replicate supported with highest number of reads per bin.

mkdir $prj/demultiplex-participants_replicates/check_crossover

### screening for identicals UMIs in different replicates

replicate=1

for bc in {01..06}; do 
#Names of included umi bins
awk -v pat="$replicate" -F ';' '/^>/ {print substr($1,2,length($1)-5)";rep0"pat""} ' $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/raconx3_medakax2/consensus_raconx3_medakax2.fa > $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_$sample_id"$bc".list
#Obtain & relabel umi bins for each replicate
awk -v pat="$replicate" -F ';' '/^>/ {$0=$1";rep0"pat""} 1' $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/umi_binning/umi_ref/umi_ref.fa > $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_ref_renamed.fa
#extract the used umi bins
seqtk subseq  $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_ref_renamed.fa $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_$sample_id"$bc".list >> $prj/demultiplex-participants_replicates/check_crossover/umi.included.fa
#Names of included umi bins linked to real name
awk -v pat="$replicate" -F ';' '/^>/ {print substr($1,2,length($1)-5)";rep0"pat"\t"substr($1,2,length($1)-1)";"$2";rep0"pat} ' $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/raconx3_medakax2/consensus_raconx3_medakax2.fa >> $prj/demultiplex-participants_replicates/check_crossover/umi_$sample_id.binsizecol.list
replicate=$((replicate+1))
done

# rename umi ref names to real names
seqkit replace -p '(.+)$' -r '{kv}' -k $prj/demultiplex-participants_replicates/check_crossover/umi_$sample_id.binsizecol.list $prj/demultiplex-participants_replicates/check_crossover/umi.included.fa > $prj/demultiplex-participants_replicates/check_crossover/umi.included.renamed.fa


# Cluster UMI pairs and detect doubles == crosssover
usearch \
  -fastx_uniques \
   $prj/demultiplex-participants_replicates/check_crossover/umi.included.renamed.fa \
  -fastaout $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.fa \
  -minuniquesize 2 \
  -sizeout \
  -uc $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.txt \
  -strand both
 
# list of non-crossover bins
awk '$1=="C" && $3<2 {print $9}' $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.txt >$prj/demultiplex-participants_replicates/check_crossover/bins.noncrossover.ok.list

# export hits of identical UMI bins found in different replicates
awk '$1=="H" {print $9 "_" $10}' $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.txt > $prj/demultiplex-participants_replicates/check_crossover/UMI_cross.txt

### Filter out UMI bin crossovers
declare -A max_values
declare -A max_columns

# Read the file line by line
while IFS= read -r line; do
  # Split the line into individual entries
  entries=(${line//_/ })

  # Extract the column values
  column1=${entries[0]}
  column2=${entries[1]}

  # Extract the bins size from umi bins in column 1 and column 2
  number1=$(echo "$column1" | awk -F 'ubs=|;rep' '{print $2}')
  number2=$(echo "$column2" | awk -F 'ubs=|;rep' '{print $2}')

  number1=$(grep -oP 'ubs=\K\d+' <<< "$column1")
  number2=$(grep -oP 'ubs=\K\d+' <<< "$column2")

# Compare the bin size and update the maximum value and corresponding column name in the arrays
if [[ -z ${max_values[$column2]} || $number1 -gt ${max_values[$column2]} ]]; then
  if [[ $number2 -gt $number1 ]]; then
    max_values[$column2]=$number2
    max_columns[$column2]=$column2
  else
    max_values[$column2]=$number1
    max_columns[$column2]=$column1
  fi
fi
done < $prj/demultiplex-participants_replicates/check_crossover/UMI_cross.txt

# Print the UMI bins with highest coverage to the output file
for key in "${!max_values[@]}"; do
  echo ${max_columns[$key]} >> $prj/demultiplex-participants_replicates/check_crossover/bins.crossover.keep.list
done

### Remove crossover bins 

# list of final umi_bins to keep
cat $prj/demultiplex-participants_replicates/check_crossover/bins.crossover.keep.list $prj/demultiplex-participants_replicates/check_crossover/bins.noncrossover.ok.list > $prj/demultiplex-participants_replicates/check_crossover/final.list

seqtk subseq $prj/demultiplex-participants_replicates/check_crossover/umi.included.renamed.fa $prj/demultiplex-participants_replicates/check_crossover/final.list >> $prj/demultiplex-participants_replicates/check_crossover/umi.included.fa

for pcr_barcode in {1..8}; do
  seqtk subseq $prj/demultiplex-participants_replicates/"$run_ID".barcode"$pcr_barcode".replicates.fa $prj/demultiplex-participants_replicates/check_crossover/final.list >> $prj/demultiplex-participants_replicates/"$run_ID".barcode"$pcr_barcode".replicates.crossoverremoved.fa
done 
  seqtk subseq $prj/demultiplex-participants_replicates/"$run_ID".undetermined.replicates.fa $prj/demultiplex-participants_replicates/check_crossover/final.list >> $prj/demultiplex-participants_replicates/"$run_ID".undetermined.replicates.crossoverremoved.fa

### Tidy data

mv $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.* /$prj/demultiplex-participants_replicates
#rm -rf $prj/demultiplex-participants_replicates/barcode*.fa $prj/demultiplex-participants_replicates/check_crossover $prj/demultiplex-participants_replicates/undetermined*.fa $prj/demultiplex-participants_replicates/$sample_id*
#find $prj/demultiplex-participants_replicates/ -type f -empty -delete






for barcode in {1..8}; do
longread_umi clustering -c "$run_ID".barcode"$barcode".replicates.crossoverremoved.fa -F $FW2 -R $RV2 -m $MIN_LENGTH -M $MAX_LENGTH -o $IN_DIR/cluster_barcode"$barcode" -i $IN_DIR -t $THREADS -b
done 

for barcode in {1..8}; do
# create empty file 
echo Cluster_ID length length_max rep01 rep02 rep03 rep04 rep05 rep06 rep07 rep08 >> $IN_DIR/cluster_barcode"$barcode"/cluster_barcode"$barcode"_replicates.txt

cd $IN_DIR/cluster_barcode"$barcode"/clusters/
# for each cluster count number of binned virus per replicate
for cluster in $(ls Cluster*) 
do
count_rep01=0
count_rep02=0
count_rep03=0
count_rep04=0
count_rep05=0
count_rep06=0
count_rep07=0
count_rep08=0
#echo $cluster
	for bin in  $(grep -o -E "^>\w+.*" $cluster | tr -d ">")
	do count_rep=count_"$(echo $bin | grep -o -E "\rep+.*" )"
	((++$count_rep))
	done
umi_len=$(seqkit stat "$cluster" | awk 'NR>1 {print $7}')
umi_len_max=$(seqkit stat "$cluster" | awk 'NR>1 {print $8}')
echo $cluster $umi_len $umi_len_max $count_rep01 $count_rep02 $count_rep03 $count_rep04 $count_rep05 $count_rep06 $count_rep07 $count_rep08 >> $IN_DIR/cluster_barcode"$barcode"/cluster_barcode"$barcode"_replicates.txt
done
done



for barcode in {1..8}; do

# create empty file 

cd $IN_DIR/cluster_barcode"$barcode"/clusters/
# for each cluster count number of binned virus per replicate
for cluster in $(ls Duo* | grep -v '^Duo.*duo\.trimmed\.fa$') 
do
count_rep01=0
count_rep02=0
count_rep03=0
count_rep04=0
count_rep05=0
count_rep06=0
count_rep07=0
count_rep08=0
#echo $cluster
	for bin in  $(grep -o -E "^>\w+.*" $cluster | tr -d ">")
	do count_rep=count_"$(echo $bin | grep -o -E "\rep+.*" )"
	((++$count_rep))
	done
umi_len=$(seqkit stat "$cluster" | awk 'NR>1 {print $7}')
umi_len_max=$(seqkit stat "$cluster" | awk 'NR>1 {print $8}')

echo $cluster $umi_len $umi_len_max $count_rep01 $count_rep02 $count_rep03 $count_rep04 $count_rep05 $count_rep06 $count_rep07 $count_rep08 >> $IN_DIR/cluster_barcode"$barcode"/cluster_barcode"$barcode"_replicates.txt
done
done




for barcode in {1..8}; do
  cd "$IN_DIR"/cluster_barcode"$barcode"
  for coverage in {5,15}; do
    ## Subset to UMI consensus sequences with min read coverage
gawk -v UBS="$coverage" '
  /^>/{
    match($0,/;ubs=([0-9]+)/, s)
    if (s[1]+0 >= UBS+0){
      print
      getline
      print
    }
  }
' "$run_ID".barcode"${barcode}".replicates.crossoverremoved.single.trimmed.fa \
> "$run_ID".barcode"${barcode}".replicates.crossoverremoved.single_${coverage}.trimmed.fa
done
done

THREADS=20




# align megabins & duo
for barcode in {1..8}; do
cd "$IN_DIR"/cluster_barcode"$barcode"
cat /data/thesisguest1/Pipe_3SEP20/db/blast/HXB2.fasta "$run_ID".barcode"${barcode}".replicates.crossoverremoved_consensus.fa "$run_ID".barcode"${barcode}".replicates.crossoverremoved.duo.trimmed.fa > "$run_ID".barcode"${barcode}"_consensus_megabins_duo.ref.fa ;
mafft --ep 0 --genafpair --thread $THREADS --reorder --adjustdirection --maxiterate 10 "$run_ID".barcode"${barcode}"_consensus_megabins_duo.ref.fa > "$run_ID".barcode"${barcode}"_consensus_megabins_duo.ref.aln.fa;
done

# align megabins & duo & singles
for barcode in {1..8}; do
cd "$IN_DIR"/cluster_barcode"$barcode"
cat /data/thesisguest1/Pipe_3SEP20/db/blast/HXB2.fasta "$run_ID".barcode"${barcode}".replicates.crossoverremoved_consensus.fa "$run_ID".barcode"${barcode}".replicates.crossoverremoved.duo.trimmed.fa "$run_ID".barcode"${barcode}".replicates.crossoverremoved.single_15.trimmed.fa > "$run_ID".barcode"${barcode}"_consensus_megabins_duo_single_15.ref.fa ;
mafft --ep 0 --genafpair --thread $THREADS --reorder --adjustdirection --maxiterate 10 "$run_ID".barcode"${barcode}"_consensus_megabins_duo_single_15.ref.fa > "$run_ID".barcode"${barcode}"_consensus_megabins_duo_single_15.ref.aln.fa;
done

# align megabins & duo & singles
for barcode in {1..8}; do
cd "$IN_DIR"/cluster_barcode"$barcode"
cat /data/thesisguest1/Pipe_3SEP20/db/blast/HXB2.fasta "$run_ID".barcode"${barcode}".replicates.crossoverremoved_consensus.fa "$run_ID".barcode"${barcode}".replicates.crossoverremoved.duo.trimmed.fa "$run_ID".barcode"${barcode}".replicates.crossoverremoved.single_5.trimmed.fa > "$run_ID".barcode"${barcode}"_consensus_megabins_duo_single_5.ref.fa ;
mafft --ep 0 --genafpair --thread $THREADS --reorder --adjustdirection --maxiterate 10 "$run_ID".barcode"${barcode}"_consensus_megabins_duo_single_15.ref.fa > "$run_ID".barcode"${barcode}"_consensus_megabins_duo_single_5.ref.aln.fa;
done