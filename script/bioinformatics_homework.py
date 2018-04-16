#bioinformatics homework 
import re 
class CDS2AA():
	pa = re.compile(r'/s+') 
	Pa = re.compile(r'[TCAG]TG')                 # 细菌起始密码子NTG 
	def __init__(self,fna,gff): 
		self.fna = fna 
		self.gff = gff 
	def N2M(self,sequence):
		hash = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
		sequence = ''.join([hash[i] for i in sequence])     #正负链转换
		return sequence[::-1]
	def Get_CDS_index(self,line):
		line = self.pa.split(line)
		CDS = (line[0], line[3], line[4], line[6])
		return CDS
	def Seq2AA(self,sequence,hash):
		AA = hash[sequence[:3]]
		if self.Pa.match(sequence[:3]):
			AA = 'M'                                 #起始密码子
			for i in range(3, len(sequence) - 3, 3):
				AA += hash[sequence[i:i + 3]]
		return AA
	def CDS2AA(self):
		fn = open(self.fna, 'r')
		gf = open(self.gff,'r')
		r = open('AA_sequence.txt', 'w')
		w = open('CDS.txt', 'w')
		hash_AA = {}  # AA hash,sequence2AA
		with open('AA.txt', 'r') as f:              
			for line in f:
				line = line.strip()
				if line:
					line = self.pa.split(line)
					hash_AA[line[0]] = line[1]      #AA hash
		hash_CDS = {}  # CDS hash,CDS2sequence
		# read fasta
		for line in fn:
			line = line.strip()
			if line.startswith('>'):
				A = self.pa.split(line)[0].replace('>', '')
				hash_CDS[A] = ''
			else:
				hash_CDS[A] += line
		fn.close()
		#############################
		for line in gf:
			line = line.strip()
			if 'CDS' in line:
				i = self.Get_CDS_index(line)
				sequence = hash_CDS[i[0]][int(i[1]) - 1:int(i[2])]
				if i[3] == '-':
					sequence = self.N2M(sequence)
				 #w.write(i[0] + '/n' + sequence + '/n')
				s = self.pa.split(line)
				p = re.compile(r'ID=(.*?);.*?Dbxref=(.*?);.*?Name=(.*?);.*?gbkey=(.*?);.*?product=(.*?);.*?protein_id=(.*?);')
				m = re.findall(p,line)
				s = s[0]+'_'+m[0][0]+m[0][2]+'/tdbxref='+m[0][1]+'/tprotein='+m[0][4]+'/tprotein_id='+m[0][5]+'/tgbkey='+m[0][3]
				w.write(s + '/n' + sequence + '/n')
				AA = self.Seq2AA(sequence, hash_AA)
				r.write(i[0] + '/n' + AA + '/n')
		w.close()
		r.close()
if __name__=='__main__':
	fna = 'GCA_000160075.2_ASM16007v2_genomic.fna'
	gff = 'GCA_000160075.2_ASM16007v2_genomic.gff'
	m = CDS2AA(fna,gff)
	m.CDS2AA()