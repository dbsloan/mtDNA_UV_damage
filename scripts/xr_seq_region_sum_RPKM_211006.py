import sys
import csv
import re
import os
import copy

def xr_seq():

	'''
	This function is design to report the FPKM of xr_seq reads given a aligned sam files and a coordmap file
	
	it reads in:
	
		xr.mt.sam file
		xr.nuc.sam file
		list of coordmap files, comma separated, no spaces
		list of genomes/chromosomes of interest, comma separated, no spaces
		
	
	It is called with:
	python3 xr_seq_region_sum.py
		followed by these 6 arguments (all separated with spaces)
		1. input path and name of mt sam file 
		2. nuc sam file
		3. list of coordmap files, comma separated, no spaces
		4. list of genomes/chromosomes. The chromosome should be a unique string that is shared between the fasta name of the chromosome used in the alignment and the name of the coordmap file (example mito or plastid)
		5. coord_like output path and file name. Essentially identical to coorddmap - but with two additional columns reporting xr seq coverage on coding and template strands for genes or F and R strands for intergenic regions
		6. read sum output path and file name. Ironically, this output was never created with any of the previous versions of the script (I forgot about it). I will write the most common read sequences here.
		7. region pivot like output path and file name. This output will summarize the coverage by region and strand
		8. gene pivot like output path and file name. This output will summarize the coverage by gene and strand
		9. rep_ID. This will come in handy when I concatenate all files togethere and sift through them in a pivot table
		10. filtered sam file, where only perfect matching reads are retained
		11. nuc genome length
		
		outputs are tab delimited txt files
		
		NOTE - the order of list items for arguments 3 and 4 must be consistent - for example I put plastid first in both cases.
	
	
		example call: python3 xr_seq_region_sum_RPKM_211006.py xr_sacc_60min_rep_2.trimmed.dedup.unique.mt.sam xr_sacc_60min_rep_2.trimmed.dedup.unique.nuc.short.sam sac_NC_001224.1.coords.chrM.txt chrM xr_sacc_60min_rep_2.trimmed.dedup.unique.mt.coordlike.txt xr_sacc_60min_rep_2.trimmed.dedup.unique.mt.readsum.txt xr_sacc_60min_rep_2.trimmed.dedup.unique.mt.cov_sum_region.txt xr_sacc_60min_rep_2.trimmed.dedup.unique.mt.cov_sum_gene.txt xr_sacc_60min_rep_2

		and in a loop form: 
		
		for file in *.sam
		do
		python3 xr_seq_region_sum_210729.py $file NC_000932.plastid.noIR.txt,NC_037304.mito.txt plastid,mito ${file%sam}coordlike.txt ${file%sam}readsum.txt ${file%sam}region_cov_sum.txt ${file%sam}gene_cov_sum.txt ${file%.unique.dedup.mt_pt.sam}
		done


	'''


	mt_sam_file = sys.argv[1]
	nuc_sam_file = sys.argv[2]
	listofcoordmap_temp= sys.argv[3]
	listofchroms_temp = sys.argv[4]
	output_handle_coordlike = sys.argv[5]
	output_hanlde_nuc_vs_org = sys.argv[6]
	output_handle_cov_sum1=sys.argv[7]
	output_handle_cov_sum2=sys.argv[8]
	rep_ID=sys.argv[9]
	filtered_sam=sys.argv[10]
	nuc_genome_length=sys.argv[11]
	
	
	
	try:
		coord_outer=open(output_handle_coordlike,'w')
	except:
		return('unable to create coord like output file')
		
	try:
		nuc_vs_org_outer=open(output_hanlde_nuc_vs_org,'w')
	except:
		return('unable to create nuc_vs_mito output file')
		
	try:
		cov_sum_outer_region=open(output_handle_cov_sum1, 'w')
		cov_sum_outer_gene=open(output_handle_cov_sum2, 'w')
	except:
		return('unable to create cov sum output file(s)')
	
	
	try:
		outro_sam=open(filtered_sam, 'w')
	except:
		return ('unable to create filtered sam file')
	
	#make matrices out of the cord map files
	
	chrom_list=listofchroms_temp.split(',')
	
	coord_list=listofcoordmap_temp.split(',')
	
	
	handle_list=[]
	reader_list=[]
	counter=0
	for chrom in chrom_list:
		counter+=1
		handle_list.append('handle'+str(counter))
		reader_list.append('reader'+str(counter))
	
	dict_of_mats={}
	
	for chrom in chrom_list:
		for coordmap in coord_list:
			if chrom in coordmap:
				dict_of_mats[chrom]=dict_of_mats.get(chrom,[])
				
				try:
					handle_list[chrom_list.index(chrom)]=open(coord_list[chrom_list.index(chrom)])
				except:
					return ('unable to open/find coord map file(s): ')
				
				reader_list[chrom_list.index(chrom)]=csv.reader(handle_list[chrom_list.index(chrom)], delimiter = '\t')
#				print (reader_list[chrom_list.index(chrom)])
			#	print (dict_of_mats)
				count=-0
				for line in reader_list[chrom_list.index(chrom)]:
					count+=1
					if count > 1: 
						#print(count)
						dict_of_mats[chrom].append([0]*8)
						dict_of_mats[chrom][count-2][0]=line[0] #col 0 is the ref pos
						dict_of_mats[chrom][count-2][1]=line[1] #col 1 is the bp
						dict_of_mats[chrom][count-2][2]=line[2] #col 2 is the type (intergenic, tRNA ext)
						dict_of_mats[chrom][count-2][3]=line[3] #col 3 is the gene name
						dict_of_mats[chrom][count-2][4]=line[4] #col 4 is the gene strand (if applicable, . for intergenic)
						dict_of_mats[chrom][count-2][5]=0 #col 5 is the xr_seq coverage on the + strand (coding strand), for intergenic sequence the F strand as reported in genbank
						dict_of_mats[chrom][count-2][6]=0 #col 6 is the xr_seq coverage on the - strand (template strand), for intergenic sequence the R strand as reported by genbank
						dict_of_mats[chrom][count-2][7]=rep_ID #col 7 is the repID - this will be useful when comparing multiple files
					
	try: 
		mt_sam_handle=open(mt_sam_file)
	except:
		return ('unable to open/find mt sam file')
		
	try: 
		nuc_sam_handle=open(nuc_sam_file)
	except:
		return ('unable to open/find nuc sam file')
		
	nuc_map_count = 0
	nuc_sam_reader=csv.reader(nuc_sam_handle,delimiter='\t')
	for line in nuc_sam_reader:
		if line[0][0]!='@':
			if line[-2][5:].isnumeric():
				nuc_map_count +=1
			









	mt_sam_reader = csv.reader(mt_sam_handle, delimiter = '\t')
	
	sam_count=0
	
	strand_dict={'F':5,'R':6}
	
	read_dict={}

	
	comp={'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	
	
	
	
	total_organelle_count={}
	for chrom in chrom_list:
		total_organelle_count[chrom]=total_organelle_count.get(chrom,0)
	
	for line in mt_sam_reader:
# 		sam_count+=1
# 		if sam_count<15:
# 			print(line)
		if line[0][0]!='@':
		
			if line[-2][5:].isnumeric():
				for col in line:
					if col == line[-1]:
						outro_sam.write(col + '\n')
					else:
						outro_sam.write(col + '\t')
				
				
				
				
				
				readlength=''
				for char in line[-2][5:]:
					readlength += char
				if line[1]=='0':
					strand='F'
				elif line[1]=='16':
					strand='R'
				
				for chrom in chrom_list:
					if chrom in line[2]:
						
						

						if int(line[3]) + len(line[9]) < int(dict_of_mats[line[2]][-1][0]):
							if dict_of_mats[line[2]][int(line[3])][3] == dict_of_mats[line[2]][int(line[3])+len(line[9])][3]:
								
								total_organelle_count[line[2]]+=1	
								strand_indexer=strand
								dict_of_mats[chrom][int(line[3])][strand_dict[strand_indexer]]+=1


	
		else:
			for col in line:
				if col == line[-1]:
					outro_sam.write(col + '\n')
				else:
					outro_sam.write(col + '\t')
				
	summary_dict_region={}
	for key in dict_of_mats:
		summary_dict_region[key]=summary_dict_region.get(key,{})
		for row in dict_of_mats[key]:
			if row[2] in summary_dict_region[key]:
				if row[2]=='intergenic':
				
					summary_dict_region[key][row[2]][0]+=row[5]
					summary_dict_region[key][row[2]][1]+=row[6]
					summary_dict_region[key][row[2]][2]+=1
				else:
					if row[4]=='1':
						summary_dict_region[key][row[2]][0]+=row[5]
						summary_dict_region[key][row[2]][1]+=row[6]
						summary_dict_region[key][row[2]][2]+=1
					elif row[4]=='-1':
						summary_dict_region[key][row[2]][0]+=row[6]
						summary_dict_region[key][row[2]][1]+=row[5]
						summary_dict_region[key][row[2]][2]+=1
	
				
			else:
				summary_dict_region[key][row[2]]=summary_dict_region[key].get(row[2],[0,0,0])
				if row[2]=='intergenic':
					summary_dict_region[key][row[2]][0]+=row[5]
					summary_dict_region[key][row[2]][1]+=row[6]	
					summary_dict_region[key][row[2]][2]+=1	
				else:
					if row[4]=='1':
						summary_dict_region[key][row[2]][0]+=row[5]
						summary_dict_region[key][row[2]][1]+=row[6]
						summary_dict_region[key][row[2]][2]+=1
					elif row[4]=='-1':
						summary_dict_region[key][row[2]][0]+=row[6]
						summary_dict_region[key][row[2]][1]+=row[5]
						summary_dict_region[key][row[2]][2]+=1
						
	summary_dict_gene={}
	for key in dict_of_mats:
		summary_dict_gene[key]=summary_dict_gene.get(key,{})
		for row in dict_of_mats[key]:
			if row[2]!= 'intergenic':
				if row[3] in summary_dict_gene[key]:
					if row[4]=='1':
						summary_dict_gene[key][row[3]][0]+=row[5]
						summary_dict_gene[key][row[3]][1]+=row[6]
						summary_dict_gene[key][row[3]][2]+=1
						summary_dict_gene[key][row[3]][3]=row[2]
					elif row[4]=='-1':
						summary_dict_gene[key][row[3]][0]+=row[6]
						summary_dict_gene[key][row[3]][1]+=row[5]
						summary_dict_gene[key][row[3]][2]+=1
						summary_dict_gene[key][row[3]][3]=row[2]
				else:
					summary_dict_gene[key][row[3]]=summary_dict_gene[key].get(row[3],[0,0,0,''])	
					if row[4]=='1':
						summary_dict_gene[key][row[3]][0]+=row[5]
						summary_dict_gene[key][row[3]][1]+=row[6]
						summary_dict_gene[key][row[3]][2]+=1
						summary_dict_gene[key][row[3]][3]=row[2]
					elif row[4]=='-1':
						summary_dict_gene[key][row[3]][0]+=row[6]
						summary_dict_gene[key][row[3]][1]+=row[5]
						summary_dict_gene[key][row[3]][2]+=1
						summary_dict_gene[key][row[3]][3]=row[2]
						
			
			
	for chrom in total_organelle_count:
		nuc_vs_org_outer.write(rep_ID + '\t' + chrom + '\t' +dict_of_mats[line[2]][-1][0] + '\t' + nuc_genome_length +'\t' +  str(total_organelle_count[chrom]) +'\t' + str(nuc_map_count) +  '\t' + str((total_organelle_count[chrom]*1000*1000000)/((total_organelle_count[chrom]+nuc_map_count)*int(dict_of_mats[line[2]][-1][0]))) + '\t' +  str((nuc_map_count*1000*1000000)/((total_organelle_count[chrom]+nuc_map_count)*(int(nuc_genome_length))) ) +'\n')					
	
	for chrom in summary_dict_region:
		for region in summary_dict_region[chrom]:
			cov_sum_outer_region.write(rep_ID + '\t' + chrom + '\t' + region + '\t' + str(summary_dict_region[chrom][region][2]) + '\t' +str(summary_dict_region[chrom][region][0] +summary_dict_region[chrom][region][1]) + '\t' + str( ((summary_dict_region[chrom][region][0] +summary_dict_region[chrom][region][1])*1000*1000000)/((total_organelle_count[chrom])*summary_dict_region[chrom][region][2]))      + '\t'+ str(summary_dict_region[chrom][region][0])+ '\t' + str(summary_dict_region[chrom][region][1]) + '\t'+ str(summary_dict_region[chrom][region][1]/summary_dict_region[chrom][region][0])+ '\t'+ str( ((summary_dict_region[chrom][region][0])*1000*1000000)/((total_organelle_count[chrom])*summary_dict_region[chrom][region][2])) + '\t'+ str( ((summary_dict_region[chrom][region][1])*1000*1000000)/((total_organelle_count[chrom])*summary_dict_region[chrom][region][2])) +    '\n' )
		cov_sum_outer_region.write(rep_ID + '\t' + chrom + '\t' + 'total' + '\t' +dict_of_mats[line[2]][-1][0] +'\t' +  str(total_organelle_count[chrom]) + '\t' + str( ((total_organelle_count[chrom])*1000*1000000)/((total_organelle_count[chrom])*int(dict_of_mats[line[2]][-1][0]))) )
	for chrom in summary_dict_gene:
		for gene in summary_dict_gene[chrom]:
	#		print (str(summary_dict_gene[chrom][gene]))

			if summary_dict_gene[chrom][gene][0] == 0:
	#			print ('HEY')
				cov_sum_outer_gene.write( rep_ID + '\t' + chrom + '\t' +str(summary_dict_gene[chrom][gene][3]) + '\t' + gene + '\t' + str(summary_dict_gene[chrom][gene][2]) + '\t' +str(summary_dict_gene[chrom][gene][0] +summary_dict_gene[chrom][gene][1]) + '\t' + str((summary_dict_gene[chrom][gene][0] +summary_dict_gene[chrom][gene][1])/summary_dict_gene[chrom][gene][2])      + '\t'+ str(summary_dict_gene[chrom][gene][0])+ '\t' + str(summary_dict_gene[chrom][gene][1]) + '\t'+ "NA"+ '\n' )
				
			else:
				#print (str(summary_dict_gene[chrom][gene]))
				cov_sum_outer_gene.write( rep_ID + '\t' + chrom + '\t' +str(summary_dict_gene[chrom][gene][3]) + '\t' + gene + '\t' + str(summary_dict_gene[chrom][gene][2]) + '\t' +str(summary_dict_gene[chrom][gene][0] +summary_dict_gene[chrom][gene][1]) + '\t' + str((summary_dict_gene[chrom][gene][0] +summary_dict_gene[chrom][gene][1])/summary_dict_gene[chrom][gene][2])      + '\t'+ str(summary_dict_gene[chrom][gene][0])+ '\t' + str(summary_dict_gene[chrom][gene][1]) + '\t'+ str(summary_dict_gene[chrom][gene][1]/summary_dict_gene[chrom][gene][0])+ '\n' )
	
							
	for key in dict_of_mats:	
		for row in dict_of_mats[key]:
			if row[2]== "intergenic":
				coord_outer.write(row[7]+ '\t'+key+ '\t'+row[0]+ '\t'+row[1]+ '\t'+row[2]+ '\t'+row[3]+ '\t'+row[4]+ '\t'+str(row[5])+ '\t'+str(row[6])+ '\t' + 'na' + '\t' + 'na' + '\n' )
			else:
				if row[4]== '-1':
					coord_outer.write(row[7]+ '\t'+key+ '\t'+row[0]+ '\t'+row[1]+ '\t'+row[2]+ '\t'+row[3]+ '\t'+row[4]+ '\t'+'na'+ '\t'+'na'+ '\t' + str(row[6]) + '\t' + str(row[5]) + '\n' )
				elif row[4] == '1':
					coord_outer.write(row[7]+ '\t'+key+ '\t'+row[0]+ '\t'+row[1]+ '\t'+row[2]+ '\t'+row[3]+ '\t'+row[4]+ '\t'+'na'+ '\t'+'na'+ '\t' + str(row[5]) + '\t' + str(row[6]) + '\n' )

# 	print(nuc_map_count)	
# 	print(total_organelle_count)
# 					
# 							

	
	
	
	
	
	
	
	
if __name__=='__main__':
	print(xr_seq())	