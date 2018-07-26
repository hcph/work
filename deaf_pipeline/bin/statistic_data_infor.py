import sys
import os

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})
    return thedict

def read_infor2dict(infor_file):
	infor_dict = {}
	f= open(infor_file,'r')
	for line in f:
		if len(line.strip("\n")) == 0 : #remove null row
			continue
		if len(line.strip("\n").split('\t')) == 2 :	# for bam.stat.xls
			infor_dict[line.split("\t")[0].strip()] = line.split("\t")[-1].strip()
		elif len(line.strip("\n").split('\t')) == 3 :	# for Basic_Statistics_of_Sequencing_Quality.txt
			infor_dict[line.split("\t")[0].strip()] = [line.split("\t")[1].strip(),line.split("\t")[2].strip()]
		elif len(line.strip("\n").split('\t')) == 5 : # for sample or index splitRate.txt
			infor_dict[line.split("\t")[0].strip()] = line.split("\t")[-1].strip()
	return infor_dict

root_dirs = sys.argv[1]
spbarcodedir = root_dirs+'/01.split_samplebarcode'
idbarcodedir = root_dirs+'/02.split_index_barcode'
filterdir = root_dirs+'/03.soapnuke'
mappingdir = root_dirs+'/04.bwa'

barcode_list = os.listdir(filterdir)
# print barcode_list
sample_split_xls = spbarcodedir +'/'+'splitRate.txt'
if not os.path.exists(sample_split_xls):
	print "Sample\tIndex\tIndex_split_rate\tRaw_Reads\tClean_Reads\tRaw_Q20\tClean_Q20\tRaw_Q30\tClean_Q30\tfilter_Pct%\tMapped_reads\tMapping_rate\tUniq_reads\tUnique_rate\tTarget_Mapped_reads\tTarget_Mapping_rate\tTarget_Uniq_reads\tTarget_Unique_rate"
	for barcode in sorted(barcode_list):
		index_list = os.listdir(filterdir+'/'+barcode)
		# print index_list
		index_split_xls = idbarcodedir+'/'+barcode+'/'+'splitRate.txt'
		index_split_dict = read_infor2dict(index_split_xls)
		# print sample,barcode,
		for index in sorted(index_list):
			## barcode_25
			tmp_barcode = barcode[:-2]+'0'+barcode[-2:]
			#
			Basic_infor_txt = filterdir+'/'+barcode+'/'+index+'/Basic_Statistics_of_Sequencing_Quality.txt'
			# print Basic_infor_txt
			Mapping_infor_xls = mappingdir+'/'+barcode+'/'+index+'.5.final.bam.stat.xls'
			# Mapping_infor_xls = mappingdir+'/'+barcode+'/'+index[:-3]+index[-2:]+'.5.final.bam.stat.xls'
			target_Mapping_infor = mappingdir+'/'+barcode+'/'+index+'.6.posi.bam.stat.xls'
			# target_Mapping_infor = mappingdir+'/'+barcode+'/'+index[:-3]+index[-2:]+'.6.posi.bam.stat.xls'
			
			if not os.path.exists(Basic_infor_txt) or os.stat(Basic_infor_txt).st_size == 0 :
				filter_dict = {'Total number of reads':['NULL (NULL)','NULL (NULL)'],'Total number of reads':['NULL (NULL)','NULL (NULL)'],'Number of base calls with quality value of 20 or higher (Q20+) (%)':['NULL (NULL)','NULL (NULL)'],'Number of base calls with quality value of 30 or higher (Q30+) (%)':['NULL (NULL)','NULL (NULL)'],'Number of filtered reads (%)':['NULL (NULL)','NULL (NULL)'],}
			else:
				filter_dict = read_infor2dict(Basic_infor_txt)
			if not os.path.exists(Mapping_infor_xls) or os.stat(Mapping_infor_xls).st_size == 0 :
				Mapped_dict = {'Mapped reads':'NULL','Mapping rate':'NULL','Uniq reads':'NULL','Unique rate':'NULL'}
			else:
				Mapped_dict = read_infor2dict(Mapping_infor_xls)	
			if not os.path.exists(target_Mapping_infor) or os.stat(target_Mapping_infor).st_size == 0 :
				target_dict = {'Mapped reads':'NULL','Mapping rate':'NULL','Uniq reads':'NULL','Unique rate':'NULL'}
			else:
				target_dict = read_infor2dict(target_Mapping_infor)

			print barcode,"\t",index,"\t",index_split_dict.setdefault(index,0),"\t",filter_dict['Total number of reads'][0].split()[0],"\t",filter_dict['Total number of reads'][1].split()[0],"\t",\
			filter_dict['Number of base calls with quality value of 20 or higher (Q20+) (%)'][0].split()[1].rstrip(")").lstrip('('),"\t",\
			filter_dict['Number of base calls with quality value of 20 or higher (Q20+) (%)'][1].split()[1].rstrip(")").lstrip('('),"\t",\
			filter_dict['Number of base calls with quality value of 30 or higher (Q30+) (%)'][0].split()[1].rstrip(")").lstrip('('),"\t",\
			filter_dict['Number of base calls with quality value of 30 or higher (Q30+) (%)'][1].split()[1].rstrip(")").lstrip('('),"\t",\
			filter_dict['Number of filtered reads (%)'][0].replace('(','').split()[1].rstrip(')'),"\t",Mapped_dict['Mapped reads'],"\t",Mapped_dict['Mapping rate'],"\t",Mapped_dict['Uniq reads'],"\t",Mapped_dict['Unique rate'],"\t",target_dict['Mapped reads'],"\t",target_dict['Mapping rate'],"\t",target_dict['Uniq reads'],"\t",target_dict['Unique rate']

else:
	sample_split_dict = read_infor2dict(sample_split_xls)
	for barcode in sorted(barcode_list):
		index_list = os.listdir(filterdir+'/'+barcode)
		# print index_list
		index_split_xls = idbarcodedir+'/'+barcode+'/'+'splitRate.txt'
		index_split_dict = read_infor2dict(index_split_xls)
		# print sample,barcode,
		for index in sorted(index_list):
			## barcode_25
			tmp_barcode = barcode[:-2]+'0'+barcode[-2:]
			#
			Basic_infor_txt = filterdir+'/'+barcode+'/'+index+'/Basic_Statistics_of_Sequencing_Quality.txt'
			# print Basic_infor_txt
			Mapping_infor_xls = mappingdir+'/'+barcode+'/'+index+'.5.final.bam.stat.xls'
			# Mapping_infor_xls = mappingdir+'/'+barcode+'/'+index[:-3]+index[-2:]+'.5.final.bam.stat.xls'
			target_Mapping_infor = mappingdir+'/'+barcode+'/'+index+'.6.posi.bam.stat.xls'
			# target_Mapping_infor = mappingdir+'/'+barcode+'/'+index[:-3]+index[-2:]+'.6.posi.bam.stat.xls'
			
			if not os.path.exists(Basic_infor_txt) or os.stat(Basic_infor_txt).st_size == 0 :
				filter_dict = {'Total number of reads':['NULL (NULL)','NULL (NULL)'],'Total number of reads':['NULL (NULL)','NULL (NULL)'],'Number of base calls with quality value of 20 or higher (Q20+) (%)':['NULL (NULL)','NULL (NULL)'],'Number of base calls with quality value of 30 or higher (Q30+) (%)':['NULL (NULL)','NULL (NULL)'],'Number of filtered reads (%)':['NULL (NULL)','NULL (NULL)'],}
			else:
				filter_dict = read_infor2dict(Basic_infor_txt)
			if not os.path.exists(Mapping_infor_xls) or os.stat(Mapping_infor_xls).st_size == 0 :
				Mapped_dict = {'Mapped reads':'NULL','Mapping rate':'NULL','Uniq reads':'NULL','Unique rate':'NULL'}
			else:
				Mapped_dict = read_infor2dict(Mapping_infor_xls)	
			if not os.path.exists(target_Mapping_infor) or os.stat(target_Mapping_infor).st_size == 0 :
				target_dict = {'Mapped reads':'NULL','Mapping rate':'NULL','Uniq reads':'NULL','Unique rate':'NULL'}
			else:
				target_dict = read_infor2dict(target_Mapping_infor)

			print barcode,"\t",index,"\t",sample_split_dict[tmp_barcode],"\t",index_split_dict[index],"\t",filter_dict['Total number of reads'][0].split()[0],"\t",filter_dict['Total number of reads'][1].split()[0],"\t",\
			filter_dict['Number of base calls with quality value of 20 or higher (Q20+) (%)'][0].split()[1].rstrip(")").lstrip('('),"\t",\
			filter_dict['Number of base calls with quality value of 20 or higher (Q20+) (%)'][1].split()[1].rstrip(")").lstrip('('),"\t",\
			filter_dict['Number of base calls with quality value of 30 or higher (Q30+) (%)'][0].split()[1].rstrip(")").lstrip('('),"\t",\
			filter_dict['Number of base calls with quality value of 30 or higher (Q30+) (%)'][1].split()[1].rstrip(")").lstrip('('),"\t",\
			filter_dict['Number of filtered reads (%)'][0].replace('(','').split()[1].rstrip(')'),"\t",Mapped_dict['Mapped reads'],"\t",Mapped_dict['Mapping rate'],"\t",Mapped_dict['Uniq reads'],"\t",Mapped_dict['Unique rate'],"\t",target_dict['Mapped reads'],"\t",target_dict['Mapping rate'],"\t",target_dict['Uniq reads'],"\t",target_dict['Unique rate']



























