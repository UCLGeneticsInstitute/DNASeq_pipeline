
cat chr1.fam > all.fam
for chr in {1..22} X Y
do
    cat chr$chr.bim
done > all.bim
(echo -en "\x6C\x1B\x01"; for chr in {1..22} X Y; do tail -c +4 chr$chr.bed; done) > all.bed

