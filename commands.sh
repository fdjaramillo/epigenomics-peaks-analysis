######################################################
# 4. EN‐TEx ATAC‐seq data: downstream analyses
######################################################

## Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files

# move to epigenomics_uvic

# run docker container
docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course

# move to folder ATAC-seq
cd ATAC-seq

# create folders to store bigBed data files and peaks analysis files
mkdir -p data/bigBed.files analyses/peaks.analysis

## Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections.

# download the metadata file
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_slims=DNA+accessibility&assay_title=ATAC-seq&biosample_ontology.term_name=stomach&biosample_ontology.term_name=sigmoid+colon&type=Experiment"

# retrieve the corresponding IDs for ATAC-seq peaks
grep -F "bigBed_narrowPeak" metadata.tsv |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11}' |\
sort -k2,2 -k1,1r > analyses/bigBed.peaks.ids.txt

# Download the bigBed files
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do 
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

## Check the integrity of the downloaded bigBed files

# retrieve original MD5 hash from the metadata
../bin/selectRows.sh <(cut -f1 analyses/bigBed.*.ids.txt) metadata.tsv | cut -f1,46 > data/bigBed.files/md5sum.txt

# compute MD5 hash on the downloaded files 
cat data/bigBed.files/md5sum.txt |\
while read filename original_md5sum; do 
  md5sum data/bigBed.files/"$filename".bigBed |\
  awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
done > tmp 
mv tmp data/bigBed.files/md5sum.txt

# make sure there are no files for which original and computed MD5 hashes differ
awk '$2!=$3' data/bigBed.files/md5sum.txt


## For each tissue, run an intersection analysis using BEDTools and report:

# Convert bigBed files to BED files
mkdir data/bed.files

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done

### 1) the number of peaks that intersect promoter regions
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -u > analyses/peaks.analysis/peaks.in.promoters."$tissue".txt
done

### 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions)
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -v -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.gene.body.bed > analyses/peaks.analysis/peaks.not.in.gene.body."$tissue".txt
done

######################################################
#5. Distal regulatory activity
######################################################

#create dir and move to it
mkdir regulatory_elements

# create other dirs
mkdir -p data/bigBed.files analyses/peaks.analysis

# Retrieve bigBed peak calling IDs for H3K4me1 H3K27ac following same procedure as for H3K4me3
for histon_mod in H3K4me1 H3K27ac; do
grep -F ${histon_mod} ../ChIP-seq/metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u >> analyses/bigBed.peaks.ids.txt; done

# Download the bigBed files
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do 
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

## Check the integrity of the downloaded bigBed files

# Retrieve MD5 hashes of the files from the metadata
../bin/selectRows.sh <(cut -f1 analyses/bigBed.peaks.ids.txt) ../ChIP-seq/metadata.tsv | cut -f1,46 > data/bigBed.files/tmp

# Compute MD5 hashes on the downloaded files
cat data/bigBed.files/tmp |\
while read filename original_md5sum; do 
  md5sum data/bigBed.files/"$filename".bigBed |\
  awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}'
done > data/bigBed.files/md5sum.txt
rm data/bigBed.files/tmp

# Make sure there are no files for which original and computed MD5 hashes differ
awk '$2!=$3' data/bigBed.files/md5sum.txt


# Convert bigBed files to BED files
mkdir data/bed.files

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do 
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done

# select open regions (for each tissue) that overlap with H3K27ac AND H3K4me1 peaks
for tissue in sigmoid_colon stomach; do
awk -v tissue=${tissue} '$2 == tissue' analyses/bigBed.peaks.ids.txt |\
sed 's/-human//g' |\
cut -f1,3 |\
paste -s -d '\t' |\

while read id_1 name_1 id_2 name_2; do
bedtools intersect -a ../ATAC-seq/analyses/peaks.analysis/peaks.not.in.gene.body.${tissue}.txt \
                   -b data/bed.files/${id_1}.bed > peaks.in.open.regions.with.${name_1}.${tissue}.tmp

bedtools intersect -a peaks.in.open.regions.with.${name_1}.${tissue}.tmp \
                   -b data/bed.files/${id_2}.bed -u > analyses/peaks.analysis/peaks.in.open.regions.with.${name_1}.and.${name_2}.${tissue}.txt

done
done

rm peaks.in.open.regions.with.*.tmp

# parse regulatory elements located on chr1 and generate a file that contains the name of the regulatory region and the start (5') coordinate of the region
for tissue in sigmoid_colon stomach; do
awk 'BEGIN{FS=OFS="\t"} $1 =="chr1" {print $4, $2}' analyses/peaks.analysis/peaks.in.open.regions.with.H3K4me1.and.H3K27ac."$tissue".txt > analyses/regulatory.elements.starts."$tissue".tsv
done

# extract protein-coding genes located on chr1 and generate a file which will store the name of the gene in the first column, and the start coordinate of the gene on the second column 
# (for genes located on the minus strand, the start coordinate will be at the 3')
awk 'BEGIN{FS=OFS="\t"} $1 =="chr1" {if ($6=="+"){start=$2} else {start=$3}; print $4, start}' ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed > analyses/gene.starts.tsv


# download python script and complete it
wget -P ../bin "https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/get.distance.py"

# For each regulatory element retrieve the closest gene and the distance to the closest gene 
for tissue in stomach sigmoid_colon; do
    cat analyses/regulatory.elements.starts.${tissue}.tsv | while read element start; do
        result=$(python ../bin/get.distance.py --input analyses/gene.starts.tsv --start "$start")
        echo -e "${element}\t${result}" 
    done > analyses/regulatoryElements.genes.distances.${tissue}.tsv
done

# Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv
for tissue in stomach sigmoid_colon; do
    echo "$tissue:"
    awk '{print $4}' analyses/regulatoryElements.genes.distances.${tissue}.tsv | \
    Rscript -e 'x <- scan("stdin", quiet=TRUE); cat("Mean:", mean(x), "\nMedian:", median(x), "\n")'
    echo ""
done



#Experiment_target (column23) = h3k4me3

#filtering:
#Assay type: DNA accessibility
#Assay title: ATAC-seq
#Biosample: sigmoid colon and stomach

