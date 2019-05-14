模拟数据的脚本：

脚本地址：https://github.com/levinyi/work/tree/master/script/simulate_fasta2fastq

该脚本用于将基因组序列的reference打断，模拟成fastq格式。

用途说明：假设你有一个流程需要做准确性验证。你没有真实的样本，你可以在ncbi中下载该物种的基因组序列，将它模拟成测序的下机结果。再用你的流程去验证。

1，脚本可设定read的错误率，默认为10%，

2，可自己设定需要模拟的测序类型，支持 【SE50，PE100】这种格式的字样。

3，模拟fastq的质量值为最高质量值 I（Sanger Phred+33 quality score）

4，可自己设定需要模拟产出的数据量 单位为条reads

代码逻辑详解：

首先是将一条fasta的序列随机出起始的位点（数据量根据参数调整，将随机出的位置存在一个列表中），用这个列表中的值和总长度画个分布图，看看随机情况如何，画图功能后续更新。

然后用于截取read长度。

随机出的位点列表是从0开始的，到长度减去测序长度 结束。

seq = fasta_seq[i:i+int(sequence_length)]

1，随机出错误的个数，（不能超多10%的个数）

2，随机出错误的位置，

3，然后替换相应位置的碱基。

4，模拟出碱基质量值（此处为伪造为最高质量值），后续修改为随机某个范围内的质量值。

用法：

python simulate_fasta2fastq.py -a ref.fasta -b PE100 -n 5000000 -e 0.1 -o ./ -c CL100000100_L02_153_1.fq.gz -d CL100000100_L02_153_2.fq.gz

python simulate_fasta2fastq.py -a ref.fasta -b SE80 -n 5000000 -e 0.1 -o ./ -c CL100000100_L02_153_1.fq.gz

说明：
若为PE类型，则模拟出的read1和read2一模一样。

脚本在有生之年会持续更新。欢迎评论留言
https://www.jianshu.com/p/4181ac1b0c4e
