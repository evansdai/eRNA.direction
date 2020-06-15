filenames=$(find *.fa)

for f in $filenames



do



echo ${f} 

RepeatMasker -pa 10 -engine rmblast -species human -alignments -dir ./masked -html -source -gff ${f}



done