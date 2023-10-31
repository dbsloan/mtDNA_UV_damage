import sys
import csv
import re
import os
import copy

def xr_seq():

	'''
	This function is design to filter sam file alignments to remove PCR duplicates - which are defined by this criteria:
		identical start position
		identical length
		identical sequence
		
		this script also removes reads which don't align (column 2 of sam file ==4)
		
		Other scripts in my pipeline remove reads with mismatches - so identical start and identical length would be sufficient - but this seems more general
	
	it reads in:
	
		xr_seq.sam file

		
	
	It is called with:
	python3 custom_dedup.py
		followed by these 6 arguments (all separated with spaces)
		1. input path and name of sam file 
		2. output file (filtered sam file)
	
		example call: python3 custom_dedup.py xr_human_NHF1_CPD_1a.trimmed.short.sam  xr_human_NHF1_CPD_1a.trimmed.short.filtered.sam
		
		
		and in a loop form: 
		
		for file in *.sam
		do
		python3 xr_seq_region_sum_210729.py $file NC_000932.plastid.noIR.txt,NC_037304.mito.txt plastid,mito ${file%sam}coordlike.txt ${file%sam}readsum.txt ${file%sam}region_cov_sum.txt ${file%sam}gene_cov_sum.txt ${file%.unique.dedup.mt_pt.sam}
		done

		note that the 210729 version of this script swapped template and codin strand coverage for genes on the R strand.
		
		
		updated all files on 210810:
		
		python3 xr_seq_region_sum_210811.py xr_sacc_05min_rep_1.trimmed.dedup.unique.mt.sam sac_NC_001224.1.coords.txt chrM xr_sacc_05min_rep_1.trimmed.dedup.unique.mt.coordlike.txt xr_sacc_05min_rep_1.trimmed.dedup.unique.mt.readsum.txt xr_sacc_05min_rep_1.trimmed.dedup.unique.mt.region_cov_sum.txt xr_sacc_05min_rep_1.trimmed.dedup.unique.mt.gene_cov_sum.txt 05min_rep_1 xr_sacc_05min_rep_1.trimmed.dedup.unique.mt.filter.sam

	
	'''


	sam_file = sys.argv[1]
	output_sam_file = sys.argv[2]

	
	
	
	try:
		outro_sam=open(output_sam_file,'w')
	except:
		return('unable to create filtered sam output file')
		

	try:
		sam_handle=open(sam_file)
	except:
		return('unable to open sam file')

		


# 	for chrom in seq_dict:
# 	
# 		print(chrom)
# 		if len(seq_dict[chrom])>100:
# 			print(seq_dict[chrom][0])
			
			
#	print(cut_site_dict)
	

	sam_reader = csv.reader(sam_handle, delimiter = '\t')
 	

	unique_read_dict={}
	
	
	for line in sam_reader:
# 		sam_count+=1
# 		if sam_count<15:
# 			print(line)

		if line[0][0]!='@':
			if line[1]!='4':
				unique_read_dict[line[2]]=unique_read_dict.get(line[2],{}) #chrom
				unique_read_dict[line[2]][line[3]]=unique_read_dict[line[2]].get(line[3],{}) #start pos
				unique_read_dict[line[2]][line[3]][len(line[9])]=unique_read_dict[line[2]][line[3]].get(len(line[9]),{}) #read length	
				unique_read_dict[line[2]][line[3]][len(line[9])][line[9]]=unique_read_dict[line[2]][line[3]][len(line[9])].get(line[9],line) #read seq
				

# 					
		else:
			for col in line:
				if col == line[-1]:
					outro_sam.write(col + '\n')
				else:
					outro_sam.write(col + '\t')
					
					
	for chrom in unique_read_dict:
		for start_pos in unique_read_dict[chrom]:
			for length in unique_read_dict[chrom][start_pos]:
				for seq in unique_read_dict[chrom][start_pos][length]:
					for col in 	unique_read_dict[chrom][start_pos][length][seq]:
						if col == 	unique_read_dict[chrom][start_pos][length][seq][-1]:
							outro_sam.write(col + '\n')
						else:
							outro_sam.write(col+'\t')
							


	
	
	
	
	
	
	
	
if __name__=='__main__':
	print(xr_seq())	