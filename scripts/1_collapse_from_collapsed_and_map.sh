# 2018-08-07
# The full Biofragmenta pipeline has already been run. 
# This folder 02b-star_alignment_all was manually created in the same biofragmenta ouptut directory
# We want to demonstrate the impact our various filtering steps have on the output
# We will take the adapter-trimmed reads and align to hg38, using the same parameters as before.
# Then we can associate the read IDs with the points in the pipeline where the read was filtered.
# To do, this, we will simply modify the MAP_STAR.sh script already generated to point to the trimmed read files

map_script="../03-star_alignment/MAP_STAR.sh"

map_script_out="MAP_STAR.sh"
outdir=`pwd`
orig_fls="orig_files.tmp"
new_fls="new_files.tmp"
new_fls_bn="new_file_bn.tmp"
new_file_basedir="../02-map_contaminates"
new_file_suffix="Processed.fa"

if [ -f ${orig_fls} ]; then
	rm ${orig_fls}
fi
if [ -f ${new_fls} ]; then
	rm ${new_fls}
fi
if [ -f ${new_fls_bn} ]; then
	rm ${new_fls_bn}
fi


awk '$0 ~ "--readFilesIn" {print $(NF-1)}' ../03-star_alignment/MAP_STAR.sh |\
	sed 's/,/\n/g' > ${orig_fls}


find ${new_file_basedir} -name "*"${new_file_suffix} |\
	tee ${new_fls} |\
	sed "s/_$new_file_suffix//g" |\
	xargs -n 1 -I {} basename {} > ${new_fls_bn}

missing_fls=`awk '{if(FNR==1 && NR>1){RN[ARGIND-1]=rn} rn=FNR;} END{for(x in RN){ RN[x]==rn ? i+=0 : i++ } print i}' ${new_fls_bn} ${new_fls} ${orig_fls}`
if [ ${missing_fls} -ne 0 ]; then
	echo "Check input fasta files. Not 1:1 match between STAR mapping file" > /dev/stderr
	exit ${missing_fls}
fi

newFLs_list=(`grep -f ${new_fls_bn} -n ${orig_fls} |\
	awk '{
		if(ARGIND==1){
			split($1, a, ":")
			orderv[FNR]=a[1]			
		} else if(ARGIND==2){
			print $1"\t"orderv[FNR]
		}
	}' /dev/stdin ${new_fls} |\
	sort -k2n,2 |\
	awk '{FNR==1 ? fmt="%s" : fmt=" %s"; printf(fmt, $1)} END{printf("\n")}'`)

awk '{
	if($1 ~ /^>/){
		split($1, a, "-")
		this_count=a[2]
	} else{
		sample_count[$1]++
		tot_reads[$1]+=this_count
	}
} END{
	i=0
	for(x in sample_count){
		i++
		print ">"x"-"sample_count[x]"-"tot_reads[x]
		print x
	}
}' ${newFLs_list[@]} > "/dev/shm/"${new_file_suffix} && mv "/dev/shm/"${new_file_suffix} ./
	
STAR \
  --runThreadN 20 \
  --genomeDir /data/genomic_data/biofragmenta_reference/Homo_sapiens/UCSC/GRCh38/STARIndex-2.5.2b \
  --parametersFiles /data/genomic_data/biofragmenta_reference/STAR_parameter_file \
  --outFileNamePrefix "All_Files_Processed_" \
  --readFilesCommand - \
  --readFilesIn ${new_file_suffix}  2>&1 | tee STAR.log 

exit

awk -F " " -v thisdir=${outdir} -v newList=${newFLs_list} '{
			if(FNR==1){
				print "#"$0
			} else if($0 ~ "--outFileNamePrefix"){
				split($0, a, " ")
				print "  --outFileNamePrefix "thisdir"/ "a[length(a)]
			} else if($0 ~ "--readFilesIn"){
				$(NF-1)=newList
				print $0
			} else if($0 ~ "--outSAMattrRGline"){
				gsub("_CONTAM_READS_Unmapped", "", $0)
				gsub("./biofragmenta-v", "../../biofragmenta-v", $0)
				print $0
			} else if($0 ~ "parametersFiles"){
				print $0
				print "--readFilesCommand - "$NF
			} else{
				print $0
			}
		}' ${map_script} |\
		sed 's/tee \(.*\)03-star_alignment.tmp/tee \102b-star_alignment_all/' > ${map_script_out}
