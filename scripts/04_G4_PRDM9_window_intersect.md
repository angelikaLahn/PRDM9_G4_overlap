# intersect G4 Cut&Tag file with PRDM9 peak file to identify shared binding sites

cd ~/Desktop/G4_PRDM9_Analysis/Analysis/ChIPseeker
pwd
ls *.bed #list only bed files

bedtools window -a HEK293_G4s_hg38.bed -b HEK293_PRDM9-B_HA_MACS2_hg38.bed -w 25 > G4_PRDM9_25bp_window.bed
