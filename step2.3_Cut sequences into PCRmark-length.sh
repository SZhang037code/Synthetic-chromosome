cd ./intron_exon_intron/
for i in *
  do
f=${i/.bed/_win.bed}
echo $i $f
bedtools makewindows -w 24 -s 12 -b $i -i srcwinnum > $f
done