#INPUT_INFO_FILE='/public/home/yyliu/ChIA-PET/V0-H3K4me3.txt'
#INPUT_ANCHOR_FILE='/home/data/givenAnchor.cluster.txt' ### The path and file name of given anchors' clusters. If you don't have this file, please input 'null'. 
INPUT_ANCHOR_FILE='/public/home/yyliu/ChIA-PET/ChIA-PET_Tool_V3/V0_V26N6_1.5kb_merged_submit_3kb_middleR7.bed'
OUTPUT_DIRECTORY='/public/home/yyliu/ChIA-PET/ChIA-PET_Tool_V3/V0-H3K4me3R3'
OUTPUT_PREFIX='V0-H3K4me3R3'

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
awk '{if($13>=2)print}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.cluster.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered
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

