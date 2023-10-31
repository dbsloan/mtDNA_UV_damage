import sys
import csv
import re
import os
import copy

def xr_seq():

	'''
	This function is design to report the nucleotide composition of all reads in terms of % of each base and frequency of all di-pyrimidines (two files). It requires the aligned reads file. 
	Given the multiple filtering steps I've been utilizing in this xr-seq pipeline - it is best to use the filtered sam files. 
	
		xrseq.filtered.sam file
	
	It is called with:
	python3 xr_seq_region_sum.py
		followed by these 3 arguments (all separated with spaces)
		1. input path and name of sam file 
		2. list of genomes/chromosomes. The chromosome should be a unique string that is identical to the chrom name used in the alignment: example NC_000932.1_plastid,NC_037304.1_mito
		3. rep_ID
		4. read composition output path and file name. 
		5. di-pyrimidine freq outputpath and ffile name.

		example call: python3 xr_seq_read_nuc_comp.py ../read_location_sum/unique_mt_pt_reads.sam NC_000932.1_plastid,NC_037304.1_mito test_lib test_lib_read_nuc_comp.txt 	
		
		'''


	sam_file = sys.argv[1]
	listofchroms_temp = sys.argv[2]
	rep_ID=sys.argv[3]
	output_file = sys.argv[4]
	output_file2 = sys.argv[5]
	
	try:
		outer_handle=open(output_file,'w')
	except:
		return('unable to create  output file')
	
	try:
		di_py_handle=open(output_file2, 'w')
	except:
		return('unable to create di-py output file')

	#make matrices out of the cord map files
	
	chrom_list=listofchroms_temp.split(',')
	
	chrom_dict={}
	di_py_dict={}
	
	for chrom in chrom_list:
		chrom_dict[chrom]=chrom_dict.get(chrom,{})
		di_py_dict[chrom]=di_py_dict.get(chrom,{})
		
	try: 
		sam_handle=open(sam_file)
	except:
		return ('unable to open/find sam file')
		
	
	sam_reader = csv.reader(sam_handle, delimiter = '\t')
	comp={'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	di_py_finder={'CC':0,'CT':1,'TC':2,'TT':3}
	for line in sam_reader:
		
		if line[0][0]!='@':
			if 'N' not in line[9]:
				if line[1] == '16':
					seq= line[9][::-1]
					if str(len(seq)) in chrom_dict[line[2]]:
						string_count=0
						for char in seq:
							string_count+=1
							chrom_dict[line[2]][str(len(line[9]))][str(string_count)][comp[char]]+=1
							chrom_dict[line[2]][str(len(line[9]))][str(string_count)]['tot']+=1
							
							
							if seq[string_count-1:string_count+1] in di_py_finder:
								di_py_dict[line[2]][str(len(seq))][str(string_count)][di_py_finder[seq[string_count-1:string_count+1]]]+=1
								di_py_dict[line[2]][str(len(seq))][str(string_count)][4]+=1
							else:
								di_py_dict[line[2]][str(len(seq))][str(string_count)][4]+=1
						
							
					else:
						chrom_dict[line[2]][str(len(seq))]=chrom_dict[line[2]].get(str(len(seq)),{})
						di_py_dict[line[2]][str(len(seq))]=di_py_dict[line[2]].get(str(len(seq)),{})
						string_count=0
				
						for char in seq:
							string_count+=1
							chrom_dict[line[2]][str(len(seq))][str(string_count)]=chrom_dict[line[2]][str(len(seq))].get(str(string_count),{'A':0,'C':0,'G':0,'T':0,'tot':0})
							di_py_dict[line[2]][str(len(seq))][str(string_count)]=di_py_dict[line[2]][str(len(seq))].get(str(string_count),[0,0,0,0,0])
							
							
							chrom_dict[line[2]][str(len(seq))][str(string_count)][comp[char]]+=1
							chrom_dict[line[2]][str(len(seq))][str(string_count)]['tot']+=1
							
							if seq[string_count-1:string_count+1] in di_py_finder:
								di_py_dict[line[2]][str(len(seq))][str(string_count)][di_py_finder[seq[string_count-1:string_count+1]]]+=1
								di_py_dict[line[2]][str(len(seq))][str(string_count)][4]+=1
							else:
								di_py_dict[line[2]][str(len(seq))][str(string_count)][4]+=1
					
		
								
							
							
				
				
				
				else:
					if str(len(line[9])) in chrom_dict[line[2]]:
						string_count=0
						for char in line[9]:
							string_count+=1
							chrom_dict[line[2]][str(len(line[9]))][str(string_count)][char]+=1
							chrom_dict[line[2]][str(len(line[9]))][str(string_count)]['tot']+=1
							
							if line[9][string_count-1:string_count+1] in di_py_finder:
								di_py_dict[line[2]][str(len(line[9]))][str(string_count)][di_py_finder[line[9][string_count-1:string_count+1]]]+=1
								di_py_dict[line[2]][str(len(line[9]))][str(string_count)][4]+=1
							else:
								di_py_dict[line[2]][str(len(line[9]))][str(string_count)][4]+=1
					else:
						chrom_dict[line[2]][str(len(line[9]))]=chrom_dict[line[2]].get(str(len(line[9])),{})
						di_py_dict[line[2]][str(len(line[9]))]=di_py_dict[line[2]].get(str(len(line[9])),{})

						string_count=0
				
						for char in line[9]:
							string_count+=1
							chrom_dict[line[2]][str(len(line[9]))][str(string_count)]=chrom_dict[line[2]][str(len(line[9]))].get(str(string_count),{'A':0,'C':0,'G':0,'T':0,'tot':0})
							di_py_dict[line[2]][str(len(line[9]))][str(string_count)]=di_py_dict[line[2]][str(len(line[9]))].get(str(string_count),[0,0,0,0,0])
							chrom_dict[line[2]][str(len(line[9]))][str(string_count)][char]+=1
							chrom_dict[line[2]][str(len(line[9]))][str(string_count)]['tot']+=1
	
							if line[9][string_count-1:string_count+1] in di_py_finder:
								di_py_dict[line[2]][str(len(line[9]))][str(string_count)][di_py_finder[line[9][string_count-1:string_count+1]]]+=1
								di_py_dict[line[2]][str(len(line[9]))][str(string_count)][4]+=1
							else:
								di_py_dict[line[2]][str(len(line[9]))][str(string_count)][4]+=1
					
		
# 	print(di_py_dict)
	
	for chrom in di_py_dict:
		for length in di_py_dict[chrom]:
			for pos in di_py_dict[chrom][length]:
				if pos != length:
					di_py_handle.write(rep_ID + '\t' + chrom + '\t' + length + '\t' + pos+ '\t'+ pos+'-'+ str((int(pos)+1)) + '\t' + str(di_py_dict[chrom][length][pos][0]/di_py_dict[chrom][length][pos][4]) + '\t' 	+ str(di_py_dict[chrom][length][pos][1]/di_py_dict[chrom][length][pos][4]) + '\t' 	+ str(di_py_dict[chrom][length][pos][2]/di_py_dict[chrom][length][pos][4]) + '\t' 	+ str(di_py_dict[chrom][length][pos][3]/di_py_dict[chrom][length][pos][4]) + '\n' 	)
	
	for chrom in chrom_dict:
		for length in chrom_dict[chrom]:
			for pos in chrom_dict[chrom][length]:
				for nt in chrom_dict[chrom][length][pos]:
					if nt != 'tot':
						outer_handle.write(rep_ID + '\t' + chrom + '\t' + length + '\t' + pos + '\t' + nt +  '\t' + str((chrom_dict[chrom][length][pos][nt])/(chrom_dict[chrom][length][pos]['tot']))  +'\n'  )	
		
			

	
	
	
	
	
if __name__=='__main__':
	print(xr_seq())	