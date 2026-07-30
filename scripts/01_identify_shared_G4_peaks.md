##################################################################################
### IDENTIFY G4 peaks shared by at least two samples using bedtools multiinter ###
##################################################################################
# Author: Dr. Angelika Lahnsteiner
# Date: 2026-07-30
##################################################################################


# use bedtools multiinter to compare all and report all peaks which are shared between at least two samples

#set directory
cd ~/Desktop/G4_Analysis/SEACR_output_G4s
pwd
ls # list all files

# Put all BED files (shortened names) in an array
FILES=(SRR*.bed)
bedtools multiinter -i "${FILES[@]}" > multiinter_results.bed

# multiinter reports the amount of overlaps per base pair. 
# now I use this information and create a bed file which regions where a peak starts and has >=2 samples overlap to the end 0f >= two samples overlap

awk '
$4 >= 2 {
  if (!active) {
    chr=$1; start=$2; end=$3; active=1
  } else if ($1 == chr && $2 == end) {
    end=$3
  } else {
    print chr, start, end
    chr=$1; start=$2; end=$3
  }
}
$4 < 2 && active {
  print chr, start, end
  active=0
}
END {
  if (active) print chr, start, end
}
' OFS="\t" multiinter_results.bed > HEK293_G4s_hg38.bed
