filenames=$(find *.bedGraph)
for f in $filenames
do

	echo ${f} 
	cat $f | grep "chr" | grep -v "chrM" | sort -k1,1 -k2,2n > "./bedGraph_chr/$f"
	# awk 'BEGIN {OFS="\t"}; {print $1,$2,$2,$4}'
	./bedGraphToBigWig "./bedGraph_chr/$f" hg19.chrom.sizes "./bw/${f}.bw"

done 

