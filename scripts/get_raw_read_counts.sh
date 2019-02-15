infl=$1
bn=`basename ${infl} | sed 's/_L00[1-4]_R1_001.*//'`

zcat ${infl} |\
awk -v bn=${bn} 'FNR%4==1 {
					tot++
					} END{
						split(bn, a, "_");
						sm=bn
						sub("_"a[length(a)], "", sm)
						 print bn"\t"sm"\t"a[length(a)]"\t"tot
						}' > ${infl/.fastq.gz}"_input_readcount.txt"
