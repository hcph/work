# coding: utf-8
import sys,os
import random
import gzip
import argparse
from Bio import SeqIO
from contextlib import nested
def Usage():
    __version__ = "1.0.0"
    print """
Program: %s
Version: 1.0
Contact: dushiyi@genomics.cn
updated: March 16 2018 created.

Description:
  Trim fastq reads.
  Attention:the first position is 1 not 0 when use -s parameter.
  If some contents is not correct,please check the script or contact with me.

Usage: python %s [options]
Options:
  -a <int>     fasta file
  -b <int>     sequence type [SE50/SE100/PE50/PE100/PE200 ...]   
  -o <dir>     outdir directory
  -f <dir>     fastq file path
  -h help

Example:
    python %s -f /yourpath/01.split_virusbarcode/barcode_16/barcode_03.fq.gz -s 17 -l  30 -o /yourpath/02.bwa/barcode_16

"""% (sys.argv[0],sys.argv[0],sys.argv[0])

parser = argparse.ArgumentParser(description='Short sample app')

parser.add_argument('-a', action="store", dest="ref",help="fasta file")
parser.add_argument('-b', action="store", dest="seq_type",type=str,help='sequence type [SE50/SE100/PE50/PE100/PE200 ...]')
parser.add_argument('-c', action="store", dest="outputname1", help='outputname1')
parser.add_argument('-d', action="store", dest="outputname2", help='outputname2')
parser.add_argument('-v', action="version", help='show version' )
argument = parser.parse_args()

def random_err(seq,sequence_length):
	i = int(0.1*float(sequence_length))
	# i 表示错误碱基个数：现规定测序错误率为10%，则测序错误个数为 10% * 测序长度。比如：测序长度为100bp，那么10%的错误率则为10bp，即每100bp中最多可以有10bp错误。
	err_num = random.randrange(i) # 这里表示随机生成一个数，用该数值表示read中有多少个错误的碱基，但是当该随机数值为0时，表示read中没有错误碱基。
	print("err_num: %s"%err_num)
	if err_num == 0: # 当read中没有错误碱基时直接返回原始序列。
		return seq
	# err_num表示本条read的中将会有多少个错误碱基，
	err_pos = [random.randrange(int(sequence_length)-1)+1 for each in range(err_num)] # 这里sequence_length -1 表示从测序长度的最末端-1位开始随机 而最后 整体加1是为了不让其实位置从0开始，而直接起始位置从1开始。
	# 从错误碱基中随机出一个列表，每个值表示错误碱基的位置。! 有重复，所以错误的位置是从0开始的	
	err_pos.sort()
	 # 对错误的位置进行从小到大排序
	print("err_pos: %s"%err_pos)
	print("rawseq: %s"%seq) 
	# 打印原始的read
	random_base = [random.choice("ATCG") for i in range(err_num)] 
	# 随机出每个位置的错误的碱基，用于替换原来的。
	print("random_base: %s"%random_base)
	for index,value in enumerate(err_pos): 
	# 开始替换错误碱基，
		newseq = seq[:value]+random_base[index]+seq[value+1:] #如果起始位置为0，会报错！，所以上面设置了err_pos的每个位置都+1，这样截取的时候就是当前位置。
		seq = newseq
		print("newseq: %s"%newseq) 
	#最后得到新的seq
	return newseq
for seq_record in SeqIO.parse(argument.ref,"fasta"):
	fasta_seq = seq_record.seq
	fasta_length = len(seq_record.seq)
# print(type(fasta_seq))

sequence_type = argument.seq_type[:2]
sequence_length = argument.seq_type[2:]

result_list = [random.randrange(fasta_length) for i in range(100)] # 表示随机生成100个（0到fasta_length）之间的随机数，并且随机数之间有可能重复。
# print result_list

if sequence_type == 'SE':
	outputname1 = os.path.basename(argument.outputname1)
	with gzip.open(outputname1,"w") as f:
		for i in result_list:
			f.write("@read_%s%s_%s/1\n"%(sequence_type,sequence_length,i))
			seq = fasta_seq[i:i+int(sequence_length)]
			seq_with_err = random_err(seq,sequence_length)
			# print(seq)
			# print(seq_with_err)
			f.write("%s\n"%seq_with_err)
			f.write("+\n")
			f.write("%s\n"%("I"*int(sequence_length))) # 造假质量值：I
elif sequence_type == 'PE':
	outputname1 = os.path.basename(argument.outputname1)
	outputname2 = os.path.basename(argument.outputname2)
	with nested(gzip.open(outputname1,"wb"),gzip.open(outputname2,"wb")) as (f1,f2):
		for i in result_list:
			f1.write("@read_%s%s_%s/1\n"%(sequence_type,sequence_length,i))
			f2.write("@read_%s%s_%s/2\n"%(sequence_type,sequence_length,i))
			seq = fasta_seq[i:i+int(sequence_length)]
			seq_with_err = random_err(seq,sequence_length)
			# print(seq)
			# print(seq_with_err)
			f1.write("%s\n"%(seq_with_err))
			f2.write("%s\n"%(seq_with_err))
			f1.write("+\n")
			f2.write("+\n")
			f1.write("%s\n"%("I"*int(sequence_length)))
			f2.write("%s\n"%("I"*int(sequence_length)))
	# print("+")
	# print("qulity")
