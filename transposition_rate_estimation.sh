#De novo TE insertion detection is highly influenced by read coverage and library size. For this reason, to estimate transposition rates, we re-ran TEd-seq while subsampling all datasets for each TE family and sample to match the sample with the lowest coverage using samples derived from the same batch of sequencing:
#For ATCOPIA93 we subsample 24131276 reads (maximum number of reads in smallest dataset (Index_13)
#For ATENSPM3 we subsample 4843888 reads (maximum number of reads in smallest dataset (Index_13)
#For VANDAL21 we subsample 800412 reads (maximum number of reads in smallest dataset (Index_13)

# We did this adding the following command to the TEd-seq short reads script (Short_read_TEDSeq_TE_germline.sh) after line 128. For example for ATCOPIA93:

reformat.sh in=$outTrim/sub_1_part1.fastq in2=$outTrim/sub_2_part1.fastq out1=$outTrim/sub_1.fastq out2=$outTrim/sub_2.fastq samplereadstarget=24131276 overwrite=true
#All the insertions were then concatenated:
cat Index_3_L1_ATCOPIA93_subsample/Index_3_L1_ATCOPIA93_subsample_TEDseq_germline_insertion_corrected.bed Index_13_L1_ATCOPIA93_subsample/Index_13_L1_ATCOPIA93_subsample_TEDseq_germline_insertion_corrected.bed Index_15_L1_ATCOPIA93_subsample/Index_15_L1_ATCOPIA93_subsample_TEDseq_germline_insertion_corrected.bed Index_17_L1_ATCOPIA93_subsample/Index_17_L1_ATCOPIA93_subsample_TEDseq_germline_insertion_corrected.bed Index_20_L1_ATCOPIA93_subsample/Index_20_L1_ATCOPIA93_subsample_TEDseq_germline_insertion_corrected.bed | sort -k 1.4n -k 2n |awk '{OFS="\t"}{print $1, $2-500,$3+500}' | bedtools merge -i - >all_insertions_ATCOPIA93.bed
#Then we intersect again to work in all the cases with the same coordinates
intersectBed -b Index_3_L1_ATCOPIA93_subsample/Index_3_L1_ATCOPIA93_subsample_TEDseq_germline_insertion.bed -a all_insertions_ATCOPIA93.bed -wa > G0_ATCOPIA93.bed
intersectBed -b Index_13_L1_ATCOPIA93_subsample/Index_13_L1_ATCOPIA93_subsample_TEDseq_germline_insertion.bed -a all_insertions_ATCOPIA93.bed -wa >G1_ATCOPIA93_present_P1.bed
intersectBed -b Index_15_L1_ATCOPIA93_subsample/Index_15_L1_ATCOPIA93_subsample_TEDseq_germline_insertion.bed -a all_insertions_ATCOPIA93.bed -wa >G1_ATCOPIA93_present_P2.bed
intersectBed -b Index_17_L1_ATCOPIA93_subsample/Index_17_L1_ATCOPIA93_subsample_TEDseq_germline_insertion.bed -a all_insertions_ATCOPIA93.bed -wa >G2_ATCOPIA93_present_P1.bed
intersectBed -b Index_20_L1_ATCOPIA93_subsample/Index_20_L1_ATCOPIA93_subsample_TEDseq_germline_insertion.bed -a all_insertions_ATCOPIA93.bed -wa >G2_ATCOPIA93_present_P2.bed
#Then we identify all the unique insertions of each file. The total number of unique insertions divided by 1000 (1000 seedlings) will give us a proxy of the transposition rate per each population. We repeated the analysis for each TE_family.

files=("G0_ATCOPIA93.bed" "G1_ATCOPIA93_present_P1.bed" "G1_ATCOPIA93_present_P2.bed" "G2_ATCOPIA93_present_P1.bed" "G2_ATCOPIA93_present_P2.bed")
for file in "${files[@]}"; do
    temp_file="combined_others.bed"
    cat "${files[@]/$file}" > "$temp_file"
    unique_count=$(grep -Fxv -f "$temp_file" "$file" | wc -l)
    result=$(echo "scale=4; $unique_count /1000" | bc | awk '{printf "%.4f", $0}')
    echo "$file $result" >> "table_transposition_rate.txt"
    rm -r -f "$temp_file"
done

