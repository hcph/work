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
parser.add_argument('-e', action="store", dest="err_rate", type=float,default=0.1,help='err_rate')
parser.add_argument('-o', action="store", dest="outputpath", help='out put path')
parser.add_argument('-n', action="store", dest="data_amount",type=int,default=5000000, help='data_amount')
parser.add_argument('-v', action="version",dest="1.0",help='show version' )
argument = parser.parse_args()

def random_err(seq,sequence_length,err_rate):
	i = int(err_rate*sequence_length)
	# i 表示错误碱基个数 = 错误率 * 测序长度。比如：测序长度为100bp，那么10%的错误率则为10bp，即每100bp中最多可以有10bp错误。
	err_num = random.randrange(i+1) 
	# 这里表示随机生成一个数，用该数值表示read中有多少个错误的碱基，但是当该随机数值为0时，表示read中没有错误碱基。此处最大值为i
	# print("err_num: %s"%err_num)
	if err_num == 0: # 当read中没有错误碱基时直接返回原始序列。
		return seq,err_num
	# err_num表示本条read的中将会有多少个错误碱基，
	err_pos = [random.randrange(sequence_length)+1 for each in range(err_num)]
	# err_pos 表示错误碱基的位置，需要随机出上面的err_num每个错误碱基的位置
	# 这里sequence_length -1 表示从测序长度的最末端-1位开始随机 而最后 整体加1是为了不让其实位置从0开始，而直接起始位置从1开始。
	# 从错误碱基中随机出一个列表，每个值表示错误碱基的位置。! 有重复,	
	err_pos = list(set(err_pos))
	 # 对错误的位置列表进行去重并从小到大排序
	# print("err_pos: %s"%err_pos)
	# print("rawseq: %s"%seq) 
	# 打印原始的read
	random_base = [random.choice("ATCG") for i in range(err_num)] 
	# 随机出每个位置的错误的碱基，用于替换原来的。这里有可能随机出来的还是原来正确的碱基。
	# print("random_base: %s"%random_base)
	for index,value in enumerate(err_pos): 
	# 开始替换错误碱基，
		newseq = seq[:value-1]+random_base[index]+seq[value:] 
		# 这里需要测试两个极端，min和max，err_pos中的两个极端值为1和50（假设sequence_length=50）因为切片时的前端值是索引，但后端值切不到，所以真实位置需要-1，
		#如果起始位置为1，那么索引中的0应该成为后端值，所以上面设置了value-1，后面部分value又代表了索引。
		# 如果起始位置为50，那么，value-1=49,49切不上，只切了前48个，+随机值+seq[50:],50又不能成为索引,因为总共的索引是0-49。
		seq = newseq
		# print("newseq: %s"%newseq) 
	#最后得到新的seq
	return newseq,err_num
for seq_record in SeqIO.parse(argument.ref,"fasta"):
	fasta_seq = seq_record.seq
	fasta_length = len(seq_record.seq)
# print(type(fasta_seq))

sequence_type = argument.seq_type[:2]
sequence_length = int(argument.seq_type[2:])

# result_list = [random.randrange(10) for i in range(100)] 
# 表示随机生成100个（0到10）之间的随机数，其中0包含在里面，10不包含在里面，并且随机数之间会有重复。从0开始也表示索引，
result_list = [random.randrange(fasta_length+1-sequence_length+1) for i in range(argument.data_amount)]
# 这里应该减去测序长度才行，否则，一旦随机出来的数值落在最后的测序长度范围内时，就没办法截取足够长的片段。
# 表示随机生成5000000个（0到fasta_length+1-sequence_length）之间的随机数，并且随机数之间可能会有重复。随机出来的值的最小值为0，最大值为fasta_length+1-sequence_length+1，推理过程见下面
# 但是如何才能让随机数具有均一性，或者能平均的覆盖到整个ref上？

if sequence_type == 'SE':
	outputname1 = os.path.abspath(argument.outputpath)+"/"+os.path.basename(argument.outputname1)
	with gzip.open(outputname1,"w") as f:
		for i in result_list:
			seq = fasta_seq[i:i+int(sequence_length)]
			# 上一行代码表示：字符串切片是截取首尾字符的索引。冒号前一位包含，冒号后一位不包含，比如：[2:8]，2截取时第二位包含，8截取时第8位不包含在里面。
			# 那么就需要 i+sequence_length的最大值为fasta_length+1，这样才能保证，ref的最后1bp不被漏掉。那么就需要i=fasta_length+1-sequence_length.
			# 此处的i是随机出来的，但是产生随机数时，random.randrange(x),此处的x是不包括x本身的，那么就需要产生随机数时的x = i+1 =fasta_length+1-sequence_length+1.
			seq_with_err,err_num = random_err(seq,sequence_length,argument.err_rate)

			f.write("@read_%s%s_%s_E%s/1\n"%(sequence_type,sequence_length,i,err_num))
			f.write("%s\n"%seq_with_err)
			f.write("+\n")
			f.write("%s\n"%("I"*int(sequence_length))) # 造假质量值：I
elif sequence_type == 'PE':
	outputname1 = os.path.abspath(argument.outputpath)+ "/" + os.path.basename(argument.outputname1)
	outputname2 = os.path.abspath(argument.outputpath)+ "/" + os.path.basename(argument.outputname2)
	with nested(gzip.open(outputname1,"wb"),gzip.open(outputname2,"wb")) as (f1,f2):
		for i in result_list:
			seq = fasta_seq[i:i+int(sequence_length)]
			seq_with_err,err_num = random_err(seq,sequence_length,argument.err_rate)
			f1.write("@read_%s%s_%s_E%s/1\n"%(sequence_type,sequence_length,i,err_num))
			f2.write("@read_%s%s_%s_E%s/2\n"%(sequence_type,sequence_length,i,err_num))
			f1.write("%s\n"%(seq_with_err))
			f2.write("%s\n"%(seq_with_err))
			f1.write("+\n")
			f2.write("+\n")
			f1.write("%s\n"%("I"*int(sequence_length)))
			f2.write("%s\n"%("I"*int(sequence_length)))
