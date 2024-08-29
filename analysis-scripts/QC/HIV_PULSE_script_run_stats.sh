#!/bin/bash 

source activate longread_umi_HIV || conda activate longread_umi_HIV 

#set working directory
wrk=/data/thesisguest1/work/UMI

##set the following parameters

# name of umi output folders
umi_dir="umi_out_all_length"
# create UMI run stat folder
mkdir $wrk/UMI_run_stats

# function for calculation operations
calc(){ awk "BEGIN { print "$*" }"; }

today=$(date +"%Y_%m_%d")

 ## Actual code
 # set up file
 rm $wrk/UMI_run_stats/umi_binning_stats_$today.txt
 echo -e "Run_ID\tBarcode\tTotal_reads\tHIV_reads\tHIV_reads_sum_len\tHIV_reads_avg_len\tnonHIV_reads\tnonHIV_reads_sum_len\tnonHIV_reads_avg_len\t%_HIV_reads/Total_reads\tDiscovered_number_of_bins\tTotal_number_of_bins\tTotal_number_of_OK_bins\t%_Total_number_of_OK_bins/Total_number_of_bins\tTotal_number_of_binned_reads\tTotal_number_of_OK_binned_reads\t%_binned_ok_reads/all_binned_reads\t%_binned_reads/tot_HIV_reads\t%_ok_binned_reads/tot_HIV_reads" >> $wrk/UMI_run_stats/umi_binning_stats_$today.txt
 
 guppy_version=6.5.7
 for run in {55..55}; 
 do prj=$wrk/data_RUN"$run"; 
 sample_id=RUN"$run"_guppy_v"$guppy_version"_sup_BC;  
 for bc in {01..24};  
 do sample_umi_dir=$prj/$sample_id"$bc"/"$umi_dir"_"$sample_id""$bc"/umi_binning/read_binning; 
sample_umi_ref_dir=$prj/$sample_id"$bc"/"$umi_dir"_"$sample_id""$bc"/umi_binning/umi_ref; 

disc_bins=$(seqkit stat $sample_umi_ref_dir/umi12c.fa -T | awk 'NR!=1 {print $4}')

hiv_reads=$(seqkit stats $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq -T | awk 'NR!=1 {print $4}'); 
nonhiv_reads=$(seqkit stats $prj/$sample_id"$bc"/$sample_id"$bc".nonHIV.fastq -T | awk 'NR!=1 {print $4}'); 

tot_bins=$(cat $sample_umi_dir/umi_binning_stats.txt | awk -F" " 'NR!=1 {print NR-1}' | tail -1);

tot_bins=$(cat $sample_umi_dir/umi_binning_stats.txt | awk -F" " 'NR!=1 {print NR-1}' | tail -1);
tot_bins_ok=$(cat $sample_umi_dir/umi_binning_stats.txt | awk -F" " 'NR!=1 {sum += gsub(/ok/,"",$17)} END{print sum} ');
all_binned_reads=$(cat $sample_umi_dir/umi_binning_stats.txt | awk -F" " 'NR!=1 {sum+= $2} END{print sum}');
all_binned_reads_ok=$(cat $sample_umi_dir/umi_binning_stats.txt | awk -F" " '$17=="ok" {sum+= $2} END{print sum}');  

tot_reads=$(calc $hiv_reads+$nonhiv_reads)
hiv_reads_per=$(calc $hiv_reads/$tot_reads*100)
reads_ok_per=$(calc $all_binned_reads_ok/$all_binned_reads*100)
bins_ok_per=$(calc $tot_bins_ok/$tot_bins*100)
reads_binned_per=$(calc $all_binned_reads/$hiv_reads*100)
reads_binned_OK_per=$(calc $all_binned_reads_ok/$hiv_reads*100)

sum_len_nonhiv=$(seqkit stats $prj/$sample_id"$bc"/$sample_id"$bc".nonHIV.fastq -T | awk 'NR!=1 {print $5}'); 
avg_len_nonhiv=$(seqkit stats $prj/$sample_id"$bc"/$sample_id"$bc".nonHIV.fastq -T | awk 'NR!=1 {print $7}');
sum_len_hiv=$(seqkit stats $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq -T | awk 'NR!=1 {print $5}'); 
avg_len_hiv=$(seqkit stats $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq -T | awk 'NR!=1 {print $7}');  

echo -e RUN"$run"'\t'BC"$bc"'\t'"$tot_reads"'\t'$hiv_reads'\t'$sum_len_hiv'\t'$avg_len_hiv'\t'$nonhiv_reads'\t'$sum_len_nonhiv'\t'$avg_len_nonhiv'\t'$hiv_reads_per'\t'$disc_bins'\t'$tot_bins'\t'$tot_bins_ok'\t'$bins_ok_per'\t'$all_binned_reads'\t'$all_binned_reads_ok'\t'$reads_ok_per'\t'$reads_binned_per'\t'$reads_binned_OK_per >> $wrk/UMI_run_stats/umi_binning_stats_$today.txt;  
done;  
done
