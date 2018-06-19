#!/ldfssz1/MGI_BIT/Public/Software/miniconda2/bin/python
# coding:utf-8
import os,sys
import gzip
import argparse
from Bio import SeqIO
import random
from threading import Thread
# import Queue 

class MyThread(Thread):

    def __init__(self, fastq):
        Thread.__init__(self)
        self.fastq = fastq

    def run(self):
        self.seq_record, self.total_read = deal_fastq(self.fastq)
        sys.stdout.write("I am dealing with %s \n"%self.fastq)

    def get_result(self):
        return self.seq_record,self.total_read

class PrintSeq(Thread):
	"""docstring for PrintSeq"""
	def __init__(self, position,number):
		Thread.__init__(self)
		self.position = position
		self.number = number
	def run(self):
		print_reads(self.position,self.number)

		
def deal_start(ratio,total_number,fq1_total_num,fq2_total_num):
	# ratio = 10:90
	first,second = ratio.split(":")	
	fastq1_reads = int(int(first)* 0.01 * total_number)
	fastq2_reads = int(int(second)* 0.01 * total_number)
	# print ("select %s reads in first fastq"%(fastq1_reads))
	# print ("select %s reads in second fastq"%(fastq2_reads))
	fastq1_pos = [random.randrange(fq1_total_num) for each in range(fastq1_reads)]
	fastq2_pos = [random.randrange(fq2_total_num) for each in range(fastq2_reads)]
	# print("select positions in first fastq is : %s"%(fastq1_pos))
	# print("select positions in second fastq is : %s"%(fastq2_pos))
	return fastq1_pos,fastq2_pos

def get_handle(file):
    if os.path.basename(file).endswith("gz"):
        handle = gzip.open(file, "rU")
    elif os.path.basename(file).endswith("fq"):
        handle = open(file, "rU")
    return handle

def deal_fastq(fastq):
	handle = get_handle(fastq)
	# 采用字典生成式 生成fastq的字典。
	seq_record = {key:value for (key,value) in enumerate(SeqIO.parse(handle,"fastq"))}
	# seq_record = [record for record in SeqIO.parse(handle,"fastq")]
	total_read = len(seq_record)
	return seq_record,total_read

def start_thread(your_thread_list):
	for thr in your_thread_list:
		thr.start()
	for thr in your_thread_list:
		if thr.isAlive():
			thr.join()

parser = argparse.ArgumentParser(description="this is a script of mix two fastq")
parser.add_argument('-a',dest='fastq1',action="store",type=str,help='fastq1 file')
parser.add_argument('-b',dest='fastq2',action="store",type=str,help='fastq2 file')
parser.add_argument('-c',dest='ratio' ,action="store",type=str,help='ratio of fastq1 to fastq2')
parser.add_argument('-d',dest='number',action="store",type=int,help='read number that you want to output')
parser.add_argument('-o',dest='output',action="store",type=str,help='output file name')

argument = parser.parse_args()

# deal fastq 
threads1 = []
t1 = MyThread(argument.fastq1)
threads1.append(t1)
t2 = MyThread(argument.fastq2)
threads1.append(t2)
start_thread(threads1)
# 获取多线程返回值
fastq1_total_number = threads1[0].get_result()[1]
fastq2_total_number = threads1[1].get_result()[1]
sys.stdout.write("The First  fastq contains : %s number of reads.\n"%fastq1_total_number)
sys.stdout.write("The Second fastq contains : %s number of reads.\n"%fastq2_total_number)
fastq1_pos,fastq2_pos = deal_start(argument.ratio,argument.number,fastq1_total_number,fastq2_total_number)
sys.stdout.write("the ratio is %s,so first number is %s,second number is %s, total number is %s\n"%(argument.ratio,len(fastq1_pos),len(fastq2_pos),len(fastq1_pos)+len(fastq2_pos)))
sys.stdout.write("now dealing with writing ...\n")
#################### mult process
'''
threads = []
# t3=Thread(target=print_reads,args=(fastq1_pos,threads1[0].get_result()[0],))
t3 = PrintSeq(fastq1_pos,threads1[0].get_result()[0])
threads.append(t3)
# t4=Thread(target=print_reads,args=(fastq2_pos,threads1[1].get_result()[0],))
t4 = PrintSeq(fastq2_pos,threads1[1].get_result()[0])
threads.append(t4)
start_thread(threads)
'''
##########################################
def function():
	if argument.output.endswith("gz"):
		f =  gzip.open(argument.output,"wb")
	else:
		f =  open(argument.output,"wb")
	return f

def print_reads(fastq_pos,fastq_dict,f):
	for each in fastq_pos:
		SeqIO.write(fastq_dict[each],f,"fastq")

output_file_handle = function()

sys.stdout.write("writing first fastq \n")
print_reads(fastq1_pos,threads1[0].get_result()[0],output_file_handle)

sys.stdout.write("writing second fastq\n")
print_reads(fastq2_pos,threads1[1].get_result()[0],output_file_handle)
f.close()
sys.stdout.write("all finished\n")