PROGRAM_DIRECTORY='/public/home/chaohe/ChIA-PET/ChIA-PET_Tool_V3-0/program'
INPUT_INFO_FILE='/public/home/yyliu/ChIA-PET/V28-H3K4me3.txt'
#INPUT_ANCHOR_FILE='/home/data/givenAnchor.cluster.txt' ### The path and file name of given anchors' clusters. If you don't have this file, please input 'null'. 
INPUT_ANCHOR_FILE='/public/home/yyliu/ChIA-PET/ChIA-PET_Tool_V3/merged_submit_3kb_middleR5.bed'
OUTPUT_DIRECTORY='/public/home/yyliu/ChIA-PET/ChIA-PET_Tool_V3/V28-H3K4me3R2'
OUTPUT_PREFIX='V28-H3K4me3R2'

#SPECIES='1' ##1:human ; 2:mouse
CYTOBAND_DATA='/public/home/chaohe/ChIA-PET/ChIA-PET_Tool_V3/chromInfo/wheat_cytoBandIdeo.txt'
#CHROM_SIZE_INFO='mh63.chromSize.txt'
#GENOME_LENGTH='3.6E8'  ### human
#GENOME_LENGTH='2.7E9'  ### mouse
CHROM_SIZE_INFO='genome_table.txt'
GENOME_LENGTH='17E9'  ###wheat
GENOME_COVERAGE_RATIO='0.8' ### the proportion of the genome covered by the reads
GENOME_INDEX='/public/home/chaohe/db/IWGSC_v1'

BWA='/public/home/software/opt/bio/software/BWA/0.7.17//bwa'
NTHREADS='8' ### number of threads used in mapping reads to a reference genome
SAMTOOLS='/public/home/software/opt/bio/software/SAMtools/1.9/bin/samtools'
BAM2BEDPE='/public/home/software/opt/bio/software/BEDTools/2.27.0//bin/bamToBed -bedpe'

MAPPING_CUTOFF='20' ### cutoff of mapping quality score for filtering out low-quality or multiply-mapped reads

MERGE_DISTANCE='2'
SELF_LIGATION_CUFOFF='8000'
EXTENSION_LENGTH='500'
PVALUE_CUTOFF_PEAK='0.00001'
PVALUE_CUTOFF_INTERACTION='0.05'

PEAK_MODE='1' ### value: 1: the peak is an enriched region; 2: the peak is a local summit
MIN_DISTANCE_BETWEEN_PEAK='500' ### minimum distance between two different peaks
MIN_COVERAGE_FOR_PEAK='5' ### minimum coverage for peaks by extended reads

####### interaction calling
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.file.Pet2Cluster1 ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.txt ${EXTENSION_LENGTH} ${PROGRAM_DIRECTORY}/${CHROM_SIZE_INFO}
LANG=C sort -k1,1 -k4,4  < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.txt >  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.sorted
if [ ${INPUT_ANCHOR_FILE} != 'null' ];then
java -Xmx10G -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.PetClusterWithGivenAnchors ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.sorted ${INPUT_ANCHOR_FILE} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster 1
else
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.PetCluster2 ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.sorted ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster
mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.cluster.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster  #改
fi
awk '{if($13>=2)print}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered
cut -f1-3,7-9,13-15 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.txt

###### calculation of p-values
cut -f1-3 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln
cut -f4-6 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet >>${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln

cut -f1-3,13 <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor1
cut -f7-9,13 <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor2
for y in anchor1 anchor2
do
    java -Xmx10G -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.shortReads.TagCountInGivenRegions ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.${y} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.${y}.tagCount.txt 1 2
done

## generate the global tag count for global density
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln | sed 's/ /\t/g' |  cut -f1 >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nTags.txt
## calculate p-value
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nTags.txt nTags.txt
paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor1.tagCount.txt  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor2.tagCount.txt |  cut -f 4,5,10  >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt data.txt
R --vanilla < ${PROGRAM_DIRECTORY}/hypergeometric.r
mv -f result.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.hypergeo.txt
paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.hypergeo.txt | cut -f 1-3,7-9,13-15,19-24 > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.withpvalue.txt
awk -v pvalue_cutoff=${PVALUE_CUTOFF_INTERACTION} '{if($13<pvalue_cutoff+0.0)print}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.withpvalue.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt
date '+%s' >> ${OUTPUT_DIRECTORY}/time.txt
awk '{print $7"\t"$8}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt |LANG=C sort |uniq -c | awk '{if($2>=10){a+=$1;b+=$1*$3;f+=$1;g+=$1*$3}else{c[$2]+=$1;d[$2]+=$1*$3;f+=$1;g+=$1*$3}}END{for(i in c){print c[i]"\t"i"\t"d[i]};print a"\t10\t"b;print f"\t11\t"g}' |LANG=C sort -k2,2n | awk 'BEGIN{print "PET counts\tNo. of clusters\tNo.intra-chrom clusters\tNo.inter-chrom clusters\tPercent of intra-chrom clusters"}{if($2==10){print ">="$2"\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}else if($2==11){print "Total\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}else{print $2"\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}}'> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.PET_count_distribution.txt

############################################################
####### peak calling
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.BindingSitesFromPETs  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.spet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak  ${EXTENSION_LENGTH}  ${SELF_LIGATION_CUFOFF}  ${MIN_COVERAGE_FOR_PEAK}  ${PEAK_MODE}  ${MIN_DISTANCE_BETWEEN_PEAK}

###### p-value calculation for peaks
# generate local tag counts for local density
awk '{center=int(($2+$3)/2); print $1"\t"$2"\t"$3"\t"center-5000"\t"center+5000"\t"$6}'   < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak | awk '{if($4 < 0){a=0}else{a=$4};if($2<0){b=0}else{b=$2}{print $1"\t"b"\t"$3"\t"a"\t"$5"\t"$6}}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K.temp
awk '{center=int(($2+$3)/2); print $1"\t"$2"\t"$3"\t"center-10000"\t"center+10000"\t"$6}' < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak | awk '{if($4 < 0){a=0}else{a=$4};if($2<0){b=0}else{b=$2}{print $1"\t"b"\t"$3"\t"a"\t"$5"\t"$6}}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.10K_10K.temp
awk '{if(ARGIND==1){a[$1]=$2}else{b=a[$1];if($5>b){c=b}else{c=$5};if($3>b){d=b}else{d=$3}{print $1"\t"$4"\t"c"\t"$2"\t"d"\t"$6}}}' ${PROGRAM_DIRECTORY}/${CHROM_SIZE_INFO} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K.temp > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K
awk '{if(ARGIND==1){a[$1]=$2}else{b=a[$1];if($5>b){c=b}else{c=$5};if($3>b){d=b}else{d=$3}{print $1"\t"$4"\t"c"\t"$2"\t"d"\t"$6}}}' ${PROGRAM_DIRECTORY}/${CHROM_SIZE_INFO} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.10K_10K.temp > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.10K_10K
rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K.temp
rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.10K_10K.temp

for y in 5K_5K 10K_10K
do
    java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.spetCountForPeaks ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.${y}  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.spet  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.${y}
done

# generate the global spet count for global density
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.spet | sed 's/ /\t/g' |  cut -f1 >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nSpets.txt
# calculate p-value
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nSpets.txt nSpets.txt
paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.5K_5K  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.10K_10K |  cut -f5,6,13  >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.txt
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.txt data.txt
R --no-save --no-readline --args genomeLengthStr=${GENOME_LENGTH}  genomeCoverageRatioStr=${GENOME_COVERAGE_RATIO} extensionLengthStr=${EXTENSION_LENGTH} < ${PROGRAM_DIRECTORY}/pois.r
mv -f result.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.pois.txt

# intra- and inter-chromosomal pet counts for peaks
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.TagCountForPeaks ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.withTagCounts
awk '{print $6"\t"$7"\t"$6+$7}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.withTagCounts >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.temp.txt

paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.pois.txt | awk -v pvalue_cutoff=${PVALUE_CUTOFF_PEAK} '{if($8<pvalue_cutoff+0.0) print}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.tsv
awk '{print $1"\t"$2"\t"$3"\t.\t"0-100*log($8)/log(10)"\t.\t"$4"\t"$5}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.tsv >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.bed
awk '{print $1"\t"int(($2+$3)/2)"\t"$6}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.tsv  >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.aln

awk '{print $1":"$2"-"$3"\t"$6}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.compact
paste  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.compact  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.pois.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.temp.txt | awk -v pvalue_cutoff=${PVALUE_CUTOFF_PEAK} '{if($4<pvalue_cutoff+0.0) print}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.long
cat ${PROGRAM_DIRECTORY}/peakHeader.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.long >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.xls

paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.pois.txt |awk -v pvalue_cutoff=${PVALUE_CUTOFF_PEAK} '{if($8<pvalue_cutoff+0.0){print $1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.FDRfiltered.txt

wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.FDRfiltered.txt | sed 's/ /\t/g' | awk '{print "Peaks from self-ligation\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt | sed 's/ /\t/g' | awk '{print "Interacting pairs\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt
date '+%s' >> ${OUTPUT_DIRECTORY}/time.txt

awk '{printf $0"\t"}' ${OUTPUT_DIRECTORY}/time.txt |awk '{printf "Linker filtering\t";printf("%.2f",($2-$1)/60);printf "\nMapping\t";printf ("%.2f",($3-$2)/60); printf "\nRemoving redundancy\t"; printf ("%.2f",($4-$3)/60); printf "\nCategorization of PETs\t"; printf ("%.2f",($5-$4)/60); printf "\nPeak calling\t"; printf ("%.2f",($7-$6)/60); printf "\nClustering\t"; printf ("%.2f",($6-$5)/60); printf "\nTotal\t"; printf ("%.2f",($7-$1)/60); printf "\n"}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.running_time.txt
rm ${OUTPUT_DIRECTORY}/time.txt

## summary
sed -n '2p' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.linker_composition_distribution.txt |awk '{print "same-linker PETs/Total PETs\t"($2+$3+$4+$6+$7+$8)/$11}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
awk -F"\t" '{printf $2"\t"}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt |awk '{print "Unique mapped same-linker PETs/Total PETs\t"$3/$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
awk -F"\t" '{printf $2"\t"}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt |awk '{print "PETs after removing redundancy/Total PETs\t"$5/$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
awk -F"\t" '{printf $2"\t"}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt |awk '{print "inter-ligation PETs/PETs after removing redundancy\t"$7/$5}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
awk '{if($8==1){print}}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt |wc -l | awk '{printf $1"\t"}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
awk '{if($8==1 && $9<1000000){print}}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt |wc -l | awk '{printf $1"\t"}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt |awk '{print $1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
sed -n '5p' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt |awk '{print "intra-chromosomal inter-ligation PETs/inter-ligation PETs\t"$1/$3"\nintra-chromosomal inter-ligation PETs within 1Mb/intra-chromosomal inter-ligation PETs\t"$2/$3}' >>  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
sed -i '5d' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt

#move.sh
sh ${PROGRAM_DIRECTORY}/move.sh ${OUTPUT_DIRECTORY} ${OUTPUT_PREFIX}

#modify linker_composition_distribution.txt, delete 0

# Statistics Report
cp -Rf ${PROGRAM_DIRECTORY}/ChIA-PET_Tool_Report ${OUTPUT_DIRECTORY}/files_for_report
cp /disk2/users/che/ChIA-PET/ChIA-PET_Tool_V3/program/Rscript_and_genome_data/${CYTOBAND_DATA} ${OUTPUT_DIRECTORY}/files_for_report/${CYTOBAND_DATA}