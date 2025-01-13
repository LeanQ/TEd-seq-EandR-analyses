#This is an example file of the analysis followed to analyze TE allele frequency for ATCOPIA93, but can be extended to all the other families replacing the TE name, the reference TE copies and the threshold to call a TE insertion (30 for ATCOPIA93, 20 for ATENSPM3 and 5 reads for VANDAL21).
#First we merge all the bed files including the reference copy by family. (Attach to the file the different bed files for the reference copies):


cat ../reference_insertions_ATCOPIA93.bed ../Index_*_L1_ATCOPIA93_TEDseq_germline_insertion_corrected.bed | sort -k1.4n -k2n | bedtools merge -i - | awk '{OFS="\t"} {print $1, $2-250, $2+250}' | bedtools merge -i - > list_putative_insertions_ATCOPIA93.bed

for file in Index_*_L1_ATCOPIA93_clip_disc-local.bam; do
    base=$(basename "$file" _clip_disc-local.bam)
    samtools view -b -h -f 128 "$file" > "${base}_clipped.bam"
done


bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_3_L1_ATCOPIA93_clipped.bam -counts > G0R1.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_4_L1_ATCOPIA93_clipped.bam -counts > G0R2.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_7_L1_ATCOPIA93_clipped.bam -counts > G1NTR1.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_8_L1_ATCOPIA93_clipped.bam -counts > G1NTR2.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_13_L1_ATCOPIA93_clipped.bam -counts > G1NTR3.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_5_L1_ATCOPIA93_clipped.bam -counts > G1NBR1.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_6_L1_ATCOPIA93_clipped.bam -counts > G1NBR2.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_15_L1_ATCOPIA93_clipped.bam -counts > G1NBR3.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_17_L1_ATCOPIA93_clipped.bam -counts > G1NTR1.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_18_L1_ATCOPIA93_clipped.bam -counts > G2NTR2.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_19_L1_ATCOPIA93_clipped.bam -counts > G2NTR3.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_20_L1_ATCOPIA93_clipped.bam -counts > G2NBR1.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_21_L1_ATCOPIA93_clipped.bam -counts > G2NBR2.counts
bedtools coverage -a list_putative_insertions_ATCOPIA93.bed  -b Index_22_L1_ATCOPIA93_clipped.bam -counts > G2NBR3.counts


###And with this we build a provisional matrix:

paste G0R1.counts G0R2.counts G1NTR1.counts G1NTR2.counts G1NTR3.counts G1NBR1.counts G1NBR2.counts G1NBR3.counts G2NTR1.counts  G2NTR2.counts  G2NTR3.counts  G2NBR1.counts  G2NBR2.counts  G2NBR3.counts  |awk '{ printf "%s\t%s\t%s", $1, $2, $3; for (i=4; i<=NF; i+=4) printf "\t%s", $i; printf "\n" }'> provisional_matrix_ATCOPIA93.txt

##Then we binarize the matrix, we use the same criteria more than 30 reads to be called for ATCOPIA93, more than 20 for ATENSPM3 and more than 5 for VANDAL21:
awk '{    for (i=4; i<=NF; i++) {if ($i < 30) {$i = 0;} else if ($i >=30) {$i = 1;} else {$i = "NA";}}print;}' provisional_matrix_ATCOPIA93.txt > binary_matrix_ATCOPIA93.txt

## And we build the final matrix:

awk '{sum=0; for(i=4; i<=NF; i++) sum+=$i; if(sum > 2) print}' binary_matrix_ATCOPIA93.txt | awk '{OFS="\t"}{print $1,$2,$3}' | intersectBed -b - -a provisional_matrix_ATCOPIA93.txt -wa| awk '{OFS="\t"}{print $1":"$2"-"$3,"ATCOPIA93",$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}'> matrix_ATCOPIA93.txt
##Then we estimate the allele frequency of each column dividing the number of reads at each position by the average total of reads at the reference positions (p.e. for ATCOPIA93: Chr3:8674540-8675098,Chr5:5629733-5630233, Chr3:10181732-10182633, Chr1:12754815-12755315). In the case of ATENSPM3 the reference positions used to normalized are (Chr1:12942357-12942857,Chr2:4394299-4394799, Chr2:4900314-4900814) and for VANDAL21 (Chr2:10001255-10001765)
# We recommend manually checking the data here as due to the high coverage in these positions it's possible that some insertions are detected at the nearby sections, that are false positives and may distortionate the analysis. These insertions needs to be taken out from the matrix. Once done it generate the table defnitive_matrix_ATCOPIA93.txt. In this case after doing it and sorting the table by coverage (from more to less) we generate the final table. 

## After doing this for each TE family, we build the allele frequency TE

cat definitive_matrix_ATCOPIA93.txt definitive_matrix_ATENSPM3.txt definitive_matrix_VANDAL21.txt > allele_frequency_matrix.txt

#This table is then split in two for the two populations:
awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$11,$12,$13}' allele_frequency_matrix.txt > allele_freq_population1.txt
awk '{OFS="\t"}{print $1,$2,$3,$4,$8,$9,$10,$14,$15,$16}' allele_frequency_matrix.txt > allele_freq_population2.txt
#Then we proceed to the analysis in R
