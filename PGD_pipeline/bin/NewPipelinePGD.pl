#! usr/bin/perl -w 
use strict;
use Getopt::Long;
use File::Basename;
use Cwd;

my $usage =<<USAGE;
**************************************************************************************************************************
	Function:     Automatic Generate Shell Script of PGD-MD pipeline
	Version:      Ver 3.0
	Date:         2017-02-16
    EarlyModifier:    Li jinliang;Huang yile;Chen dayang
	Contact:      chendayang,<chendayang\@genomics.org.cn>,BGI-SZ
	               
	Usage:        perl $0 [options]
		        Step
			   		1 >> Raw FastQ file Filter(Quality,N-base,Adapter), BWA Alignment[-I -i 15 -L -k 2 -l 31 -t 4],Samtools operation;
					2 >> GATK: realignment; 
					3 >> GATK: SNP Indel Calling;
					4 >> GATK: SNP Indel Correcting
					5 >> Annotation Haplotyping and making report;

               		-QuaSym         1,stand for Solexa,which is the system of quality based on the Value [64], ASCII [@]
                	                2,stand for Sanger,which is the system of quality based on the Value [33], ASCII [i]
			-Memory <int>   qsub's vf options;
			-Rmdup  <int>	[Y] rm duplication/[N] don't rm duplication;
			-Gene   <STR>   Disease gene name. If undefined, it will not generate report file;
			-Target <Dir>   Directory of Target File,must content capture and gene region bed file;
					[Nemblegen: /hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/bed/PGD-Nemblegene/
					 BGI-V1:    /hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/bed/PGD-PGD-BGI-V1/
					 BGI-V2-ChipA:    /hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/bed/PGD-BGI-V2-Chip1/
					 BGI-V2-ChipB:    /hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/bed/PGD-BGI-V2-Chip2/]
       			-Help           Display this help information to screen;
			-FQ	<Dir>   Directory of original sequenceing data;
			-Lib	<File>  File which content the Sample and Library information; [Format: Sample	Lib]
			-Auto	<Con>	wheather automatic run the whole pipeline or not[N];
			-OutDir <Dir>   Output Directory;
			-Queue  <STR>   Provide queue information;

	Example:      perl $0 -Lib Lib -FQ FQ-Dir -Target Bed-Dir -Gene Disease-gene -Quasym [1|2] -Rmdup [Y|N] -Memory Int -Auto [Y|N] -Output OutDir -Queue queue -User chendayang  
**************************************************************************************************************************
USAGE

my($Depth,@Sample,$user,$h,$QuaSym,$Merge,$Lib,$FQ,$Tar,$Gene,$Memory,$Rmdup,$Output,$Auto,$queue);
GetOptions(
	'Help'    => \$h,
	'QuaSym=i'=> \$QuaSym,
	'Lib=s'   => \$Lib,
	'Target=s'=> \$Tar,
	'Gene=s'  => \$Gene,
	'Memory=i'=> \$Memory,
	'Rmdup=s' => \$Rmdup,
	'OutDir=s'=> \$Output,
	'Auto=s'  => \$Auto,
    'User=s' =>\$user,
	'Depth=i' => \$Depth,
	'Queue=s' => \$queue
);
if(defined $h){print $usage};
if(defined $h){die;}

#**************************************************************************************************************
#**********************************************    Global variate       ***************************************
#**************************************************************************************************************
$Depth  ||= 10;
$Memory ||= 2;
$Rmdup  ||= 'Y';
$QuaSym ||= 2;
$user ||='chendayang';
$Output ||= getcwd $Lib;
$Auto   ||= 'N';
$queue  ||= 'P17Z18000N0042';

my $Ref     = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/database/GATK/ucsc.hg19.fasta';
my $bin     = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/bin/New_Pipeline';
my $soft    = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft';
my $samtools= '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/samtools-1.3.1/samtools';
my $samtool = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/samtool';
my $enrich  = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/enrichstat';
my $cov     = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/soap.coverage';
my $BWA     = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/bwa/bwa';
my $GATK    = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/GenomeAnalysisTK.jar';
my $picard  = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/picard';
my $dbSNP   = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/database/GATK/dbsnp_138.hg19.vcf';
my $HapMap  = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/database/GATK/hapmap_3.3.hg19.sites.vcf';
my $Omni    = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/database/GATK/1000G_omni2.5.hg19.sites.vcf';
my $Mills   = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/database/GATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf';
my $Nregion = '/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/database/hg19/N_region.txt';


#**************************************************************************************************************
#************************************************     Executing  **********************************************
#**************************************************************************************************************
`mkdir $Output/sh` unless(-d $Output.'/sh'); 
open LIB,"<$Lib" or die "Can\'t open Lib file !!! Please cheak it\t $usage";
my(%Sample,%CleanFQ,@RawFQ);
while(<LIB>){
	chomp;
	my($Sample,$LB,@RawFQ)=split;
	$Sample{$Sample}=$LB;
	my$Label="\"\@RG\\tID:PGD-MD\\tSM:$Sample\\tLB:$LB\\tPL:BGI-SEQ500\"";
	my$Out=$Output.'/'.$Sample;
        chomp(@RawFQ);
  #	my$RawFQDir=dirname($RawFQ[0]);
	`mkdir $Out` unless(-d $Out);
#**************************************************************************************************************
#************************************************     Filtering  **********************************************
#**************************************************************************************************************
         
	my$alignment_shell=$Output.'/sh/step1-'.$Sample.'-Alignment.sh';
	open OUT,"> $alignment_shell";
        foreach my $n(0..$#RawFQ){
	print OUT "echo \"Start to Filter_Raw_FastQ $n !\" \n";
    print OUT "perl $bin/Filter_raw_data.pl -Fq1 $RawFQ[$n] -Q_cut 20 -Res $Out -N 5 -Q_system $QuaSym  -Rsta Yes && echo Filter_Raw_FastQ is completed at `date`! \n";
    if ($RawFQ[$n]=~/gz/){
    $CleanFQ{$n}[0] = basename($RawFQ[$n]),$CleanFQ{$n}[0] =~ s/_1\.fq\.gz/_1\.clean\.fq\.gz/;
	$CleanFQ{$n}[1] = basename($RawFQ[$n]),$CleanFQ{$n}[1] =~ s/_1\.fq\.gz/_2\.clean\.fq\.gz/;
    }else {
        $CleanFQ{$n}[0] = basename($RawFQ[$n]),$CleanFQ{$n}[0] =~ s/_1\.fq/_1\.clean\.fq\.gz/;
        $CleanFQ{$n}[1] = basename($RawFQ[$n]),$CleanFQ{$n}[1] =~ s/_1\.fq/_2\.clean\.fq\.gz/;
          }
         }
	#print OUT "perl $bin/Illumina_merge_Lines.pl\t";
        #map {print OUT "$Out/$CleanFQ{$_}[0] \t";}(0..$#RawFQ);
        #print OUT "\n";
        #print OUT "perl $bin/Illumina_merge_Lines.pl\t";
        #map {print OUT "$Out/$CleanFQ{$_}[1] \t";}(0..$#RawFQ);
        #print OUT "\n";
	#map {print OUT "rm $Out/$CleanFQ{$_}[0] $Out/$CleanFQ{$_}[1] \t"}(1..$#RawFQ);
        #print OUT "\n";

my @num=0..$#RawFQ;

        
        
#**************************************************************************************************************
#************************************************     Alignment  **********************************************
#**************************************************************************************************************

	foreach my$num(@num){
        if($QuaSym == 1){
                print OUT "echo \"Start to BWA alignment!round_$num   \"\n";
                 print OUT "$BWA aln -I -i 15 -L -k 2 -l 31 -t 12 $Ref $Out/$CleanFQ{$num}[0] > $Out/temp_$num\_1.sai && echo sai-1 finished!`date`\n";
                 print OUT "$BWA aln -I -i 15 -L -k 2 -l 31 -t 12 $Ref $Out/$CleanFQ{$num}[1] > $Out/temp_$num\_2.sai && echo sai-2 finished!`date`\n";
        }elsif($QuaSym == 2){
                print OUT "echo \"Start to BWA alignment round_$num  !\" \n";
                 print OUT "$BWA aln -i 15 -L -k 2 -l 31 -t 12 $Ref $Out/$CleanFQ{$num}[0] > $Out/temp_$num\_1.sai && echo sai1 finished!`date`\n";
                 print OUT "$BWA aln -i 15 -L -k 2 -l 31 -t 12 $Ref $Out/$CleanFQ{$num}[1] > $Out/temp_$num\_2.sai && echo sai2 finished!`date`\n";
        }
                my $Label="\"\@RG\\tID:PGD-MD-$Sample\\tSM:$Sample\\tLB:$LB\\tPL:BGISEQ-500\"";
                print OUT "echo \"Start to BWA sampe!\" \n";
                print OUT "$BWA sampe -a 1000 -r $Label $Ref $Out/temp_$num\_1.sai $Out/temp_$num\_2.sai $Out/$CleanFQ{$num}[0] $Out/$CleanFQ{$num}[1]  > $Out/temp_$num.sam && echo sam finished `date`!\n";
              #  print OUT "perl $bin/Sam_sta.V1.03.pl -sam $Out/temp_$num.sam -label Y -map PE -insert Y -dupd Y -out $Out && echo finished mapping statistical `date`\n";
                print OUT "perl $bin/Sam_sta.V1.03.pl -sam $Out/temp_$num\.sam  -map PE -insert Y -dupd Y -out $Out && echo finished mapping statistical `date`\n";
                print OUT "perl $bin/Filter_Sam_insertsize_uniq.pl -I $Out/temp_$num.sam -O $Out\/${Sample}_$num\.sam -minI 50 -maxI 500 && echo extract uniq is completed! `date`\n";
                # print OUT "rm $Out/1.sai $Out/2.sai $Out/$CleanFQ{$num}[0] $Out/$CleanFQ{$num}[1] \n";
                     }
         # print OUT "rm\t";
         #      map{print OUT "$Out/temp_$_.sam \t"}(@num);
         #      print OUT "\n";
           foreach $_(@num){
                print OUT "echo \"Start to Samtools process!\" \n";
                print OUT "$samtools view -ubSt $Ref.fai $Out/${Sample}_$_\.sam > $Out/temp_$_.bam && echo samtools view is completed!`date`\n";
                print OUT "$samtools sort -m 1000000000  $Out/temp_$_.bam > $Out/temp_$_.bam.sort.bam && echo samtools sort is completed!`date`\n";
                print OUT "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -jar $picard/SortSam.jar I=$Out/temp_$_.bam.sort.bam  O=$Out/temp_$_.bam.sort.picard.bam SO=coordinate VALIDATION_STRINGENCY=SILENT\n";
                print OUT "$samtools view $Out/temp_$_.bam.sort.picard.bam -H >$Out/temp_$_.bam.sort.header\n";
                print OUT "$samtools view $Out/temp_$_.bam.sort.picard.bam |awk '\$5>0' |cat $Out/temp_$_.bam.sort.header  - | $samtools view -Sb - > $Out/${Sample}_$_\.sort.header.bam\n";
                print OUT "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -jar $picard/MarkDuplicates.jar INPUT=$Out/${Sample}_$_\.sort.header.bam OUTPUT=$Out/${Sample}_$_\.sort.bam METRICS_FILE=$Out/${Sample}_$_\.sort.dedup.metrics VALIDATION_STRINGENCY=SILENT  PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates TMP_DIR=$Out\n"; 
#                print OUT "$samtools rmdup $Out/temp_$_.bam.sort.bam $Out/${Sample}_$_\.sort.bam && echo samtools rmdup is completed!`date`\n" if($Rmdup=~/Y/i);
                #print OUT " rm -rf $Out/${Sample}_$_\.sam  $Out/temp_$_.bam $Out/temp_$_.bam.sort.bam $Out/temp_$_.bam.sort.picard.bam $Out/${Sample}_$_\.sort.header.bam\n";
                }
                print OUT "PATH=$soft:\$PATH\n";
                print OUT "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx5g -Djava.io.tmpdir=/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/bin/javatmp -jar $picard/MergeSamFiles.jar\t";
                map{print OUT "INPUT=$Out/${Sample}_$_\.sort.bam \t"}(@num);
                print OUT "OUTPUT=$Out/$Sample\.bam VALIDATION_STRINGENCY=SILENT\n";
                #print OUT "$samtools rmdup $Out/temp.bam.sort.bam $Out/$Sample\.bam && echo samtools rmdup is completed!`date`\n" if($Rmdup=~/Y/i);
                print OUT "$samtools index $Out/$Sample\.bam && echo $Sample samtools index has completed at `date`! >> $Output/Project.log\n";
                print OUT "$enrich $Out/$Sample\.bam $Tar/Target_region.bed $Out \n";
                print OUT "$samtools view -h $Out/$Sample\.bam > $Out/$Sample\.sam\n";
                # print OUT "perl $bin/Sam_sta.V1.03.pl -sam $Out/$Sample\.sam  -map PE -insert Y -dupd Y -out $Out && echo finished mapping statistical `date`\n";
                print OUT "$cov -cvg -i $Out/$Sample\.sam -sam -refsingle $Ref -o $Out/WholeChr.cov -addn $Nregion -cdsinput $Tar/Target_region.bed -cdsdetail $Out/TargetRegion.cov -onlyuniq && echo Statistic of coverage is finished! `date`\n";
                print OUT "cp $Out/$Sample\.map $Out/temp\.map\n";
                print OUT "rm -rf  $Out/$Sample\.sam\t";
                map{print OUT "$Out/${Sample}_$_\.sort.bam  $Out/${Sample}_$_\.sam\t"}(@num);
                print OUT "\n";
#**************************************************************************************************************
#************************************************     Realignment  ********************************************
#**************************************************************************************************************

		`mkdir $Out/chr` unless(-e $Out.'/chr');
		my @CHR=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY');
		my$outfile=$Output.'/sh/step2-'.$Sample.'-Realignment.sh';
		open GATK,"> $outfile";
		foreach my$chr(@CHR)
		{
			my$ChrTar=$Tar.'/chr/'.$chr.'.bed';
			print GATK "$samtools view -hb $Out/$Sample\.bam $chr > $Out/chr/$chr.bam && echo split $chr bam done`date`\n";
        		print GATK "$samtools view -H $Out/chr/$chr.bam > $Out/chr/$chr.reheader.sam\n";
        		print GATK "$samtools reheader $Out/chr/$chr.reheader.sam $Out/chr/$chr.bam > $Out/chr/$chr.reheader.bam && echo Reheader-is-Completed`date`\n";
        		print GATK "$samtools sort  $Out/chr/$chr.reheader.bam > $Out/chr/$chr.sort.bam\n";
        		print GATK "$samtools index $Out/chr/$chr.sort.bam\n";
                #	print GATK "rm $Out/chr/$chr.bam $Out/chr/$chr.reheader.sam $Out/chr/$chr.reheader.bam\n";
	
			print GATK "echo \"Start to GATK intervals!\"\n";
			print GATK "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar $GATK -l INFO -T RealignerTargetCreator -R $Ref -L $ChrTar -I $Out/chr/$chr.sort.bam -o $Out/chr/$chr.intervals && echo Making Inervals is finished!`date`\n";
			print GATK "echo \"Start to GATK Realigner!\" \n";
			print GATK "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar $GATK -l INFO -T IndelRealigner -R $Ref -L $ChrTar -I $Out/chr/$chr.sort.bam -o $Out/chr/$chr.realign.bam -targetIntervals $Out/chr/$chr.intervals && echo $Sample $chr Realigning has finished!`date` >> $Output/Project.log\n";
            #		print GATK "rm $Out/chr/$chr.intervals $Out/chr/$chr.sort.bam $Out/chr/$chr.sort.bam.bai\n";
			#print ">Recal" ,`date '+ %X'`;
			#print GATK "/ifs2/S2PH/PGS-MD/soft/samtool sort -m 200000000 $out/GATK/temp.realign.bam $out/GATK/temp.realign.bam.sort && echo samtools sort is completed!\n";
			#print GATK "/ifs2/S2PH/PGS-MD/soft/samtool index $out/temp.realign.bam.sort.bam \n";
			#`/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar /ifs2/S2PH/PGS-MD/soft/GenomeAnalysisTK.jar -l INFO -T CountCovariates -nt 8 -R $Ref1 -knownSites $dbSNP -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -I $out\/temp\.realign\.bam -recalFile $out\/temp\.recal`;
			#`/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar /ifs2/S2PH/PGS-MD/soft/GenomeAnalysisTK.jar -l INFO -R $Ref1 -noOQs -T TableRecalibration -I $out\/temp\.realign\.bam -recalFile $out\/temp\.recal -o $out\/temp\.realign\.recal\.bam`;
			#print ">Reduce reads" ,`date '+ %X'`;
			#`java -Xmx3g -jar /ifs2/S2PH/PGS-MD/soft/GenomeAnalysisTK.jar -l INFO -R $Ref1 -T ReduceReads -I $out\/temp\.realign\.bam -o $out\/temp\.reduce\.bam`;
			#print ">Base recalibrator" ,`date '+ %X'`;
			#`java -Xmx3g -jar /ifs2/S2PH/PGS-MD/soft/GenomeAnalysisTK.jar -l INFO -R $Ref1 -T BaseRecalibrator -I $out\/temp\.realign\.bam -knownSites $dbSNP -knownSites -HapMap -knownSites $Omni -knownSites $Mills -o $out\/temp\.recal\.grp`;
		}
		close GATK;
}

################step3	
#**************************************************************************************************************
#************************************************   Variate Detecting  ****************************************
#**************************************************************************************************************

	my @CHR=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY');
	`mkdir $Output/chr` unless(-d $Output.'/chr');
	foreach my$chr(@CHR)
              	{
		my$ChrTar=$Tar.'/chr/'.$chr.'.bed';
		my$outfile=$Output.'/sh/step3-Genotyper-'.$chr.'.sh';
		open  GATK,"> $outfile";
		print GATK "echo \"Start to GATK SNP calling\" \n";
		#print GATK "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar $GATK -l INFO -T HaplotypeCaller -R $Ref -L $ChrTar --dbsnp $dbSNP -stand_call_conf 20 -stand_emit_conf 4 -ip 200 -o $Output/chr/$chr.vcf -A IndelType -A QualByDepth -A HaplotypeScore -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A DepthOfCoverage";
		print GATK "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar $GATK -l INFO -T UnifiedGenotyper -R $Ref -L $ChrTar --dbsnp $dbSNP  -o $Output/chr/$chr.vcf -glm BOTH";
		map{print GATK " -I $Output/$_/chr/$chr.realign.bam "} keys %Sample;
		print GATK " && echo $chr finished HaplotypeCaller `date`>> $Output/Project.log\n";
		#print GATK "java -Xmx10g -jar /ifs2/S2PH/PGS-MD/soft/GenomeAnalysisTK.jar -l INFO -R $Ref -T UnifiedGenotyper -glm BOTH -nt 8 -I $out/temp\.realign\.bam\.sort\.bam -o $out/temp\.vcf -A AlleleBalance -A DepthOfCoverage -A HaplotypeScore -stand_call_conf 20 -stand_emit_conf 4 --dbsnp $dbSNP -dcov 30 && echo SNP calling is finished!\n";
		#print ">Merge all vcfs",`date '+ %X'`;
	}

#**************************************************************************************************************
#************************************************     Correcting  **********************************************
#**************************************************************************************************************

	open GATK,">$Output/sh/step4-Correct.sh";
	print GATK "echo Start to GATK combine Variantion call `date` \n";
	print GATK "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -cp $GATK org.broadinstitute.gatk.tools.CatVariants  -R $Ref  -V $Output/chr/chr1.vcf -V $Output/chr/chr2.vcf -V $Output/chr/chr3.vcf -V $Output/chr/chr4.vcf -V $Output/chr/chr5.vcf -V $Output/chr/chr6.vcf -V $Output/chr/chr7.vcf -V $Output/chr/chr8.vcf -V $Output/chr/chr9.vcf -V $Output/chr/chr10.vcf -V $Output/chr/chr11.vcf -V $Output/chr/chr12.vcf -V $Output/chr/chr13.vcf -V $Output/chr/chr14.vcf -V $Output/chr/chr15.vcf -V $Output/chr/chr16.vcf -V $Output/chr/chr17.vcf -V $Output/chr/chr18.vcf -V $Output/chr/chr19.vcf -V $Output/chr/chr20.vcf -V $Output/chr/chr21.vcf -V $Output/chr/chr22.vcf -V $Output/chr/chrX.vcf -V $Output/chr/chrY.vcf -out  $Output/temp.vcf  -assumeSorted && echo finished combined! \n";
	print GATK "echo Start to GATK SNP Variantion call `date`\n";
	print GATK "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar $GATK -l INFO -R $Ref -T VariantRecalibrator -input $Output/temp.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HapMap -resource:omni,known=false,training=true,truth=false,prior=12.0 $Omni -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $dbSNP -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ --maxGaussians 3 -recalFile $Output/temp.SNP.recal -tranchesFile $Output/temp.SNP.tranches -rscriptFile $Output/temp.SNP.plot.R -mode SNP && echo SNP Variation calling is finished!\n";
	print GATK "echo Start to GATK SNP filter `date` \n";
	print GATK "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar $GATK -l INFO -R $Ref -T ApplyRecalibration -input $Output/temp.vcf --ts_filter_level 99.0 -recalFile $Output/temp.SNP.recal -tranchesFile $Output/temp.SNP.tranches -mode SNP -o $Output/temp.SNP-VQSR.vcf && echo SNP filter is finished!\n";
	print GATK "echo Start to GATK Indel Variantion call `date` \n";
	print GATK "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar $GATK -l INFO -R $Ref -T VariantRecalibrator -input $Output/temp.SNP-VQSR.vcf -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $Mills -an QD -an FS  -an ReadPosRankSum --maxGaussians 2 -std 10.0   -recalFile $Output/temp.INDEL.recal -tranchesFile $Output/temp.INDEL.tranches -rscriptFile $Output/temp.INDEL.plot.R -mode INDEL && echo Indedl variation calling is finished!\n";
	print GATK "echo Start to GATK Indel filter `date`\n";
	print GATK "/hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/jre1.8.0_121/bin/java -Xmx6g -jar $GATK -l INFO -R $Ref -T ApplyRecalibration -input $Output/temp.SNP-VQSR.vcf --ts_filter_level 99.0 -recalFile $Output/temp.INDEL.recal -tranchesFile $Output/temp.INDEL.tranches -mode INDEL -o $Output/Family.SNP.INDEL.vcf && echo Indel filter if finished!`date`\n";
	print GATK "rm $Output/temp.SNP.recal $Output/temp.SNP.recal.idx $Output/temp.vcf.idx $Output/temp.SNP.tranches $Output/temp.SNP-VQSR.vcf $Output/temp.SNP-VQSR.vcf.idx $Output/temp.INDEL.recal $Output/temp.INDEL.recal.idx $Output/temp.INDEL.tranches && echo Project SNP-INDEL-VQSR has done!>>$Output/Project.log \n"; 

#**************************************************************************************************************
#************************************************     Annotation **********************************************
#**************************************************************************************************************

	`mkdir $Output/anno` unless(-d "$Output/anno");
	open ANNO,">$Output/sh/step5-Anno.sh" or die $!;
	print ANNO "perl $soft/convert2annovar.pl --includeinfo --format vcf4 $Output/Family.SNP.INDEL.vcf  > $Output/anno/Family.SNP.INDEL.anno \n";
	print ANNO "perl $soft/summarize_annovar.pl --remove --buildver hg19 --verdbsnp 138 $Output/anno/Family.SNP.INDEL.anno /hwfssz1/ST_MCHRI/REHEAL/USER/chendayang/00.rawdata/ifs1/bin/PGD/soft/humandb  && echo Mutation annovar if finished!`date`\n";
	print ANNO "perl $bin/Exome_SNP_in_target_gene.pl -v $Output/Family.SNP.INDEL.vcf -G $Tar/target.gene -I $Output/anno/Family.SNP.INDEL.anno.exome_summary.csv -O $Output/anno/Exome_mutation.csv \n";
	print ANNO "perl $bin/Haplotype_based_merged_GATK_V3.pl -R $Tar/target.gene -I $Output/Family.SNP.INDEL.vcf -O $Output/Hap -DT $Depth \n";
        print ANNO "perl $bin/Origin_SNP.pl  $Output/Family.SNP.INDEL.vcf  $Tar/target.gene $Output/Hap/Family.snp.xls\n";
	if(defined $Gene){print ANNO "perl $bin/PGD_MD_summary_report.pl -r $Tar/target.gene -l $Lib -a $Output/anno/Exome_mutation.csv -g $Gene -s $Output/Hap/Family.score.xls -o $Output/PGD-MD-Project-case-report.xls && echo Congratelatopns, all progresses were done at `date`! >> $Output/Project.log\n";}

#**************************************************************************************************************
#*********************************************  Submit all mission    *****************************************
#**************************************************************************************************************

if($Auto=~/Y/){
	my$record,my@job,my$job_aln,my$job_realn,my$job_gt,my$job_var;
	map{$record=`cd $Output/sh/ && qsub -cwd -P $queue -l p=5  -l  num_proc=1  -q mgi.q -l vf=5g step1-$_-Alignment.sh`;push @job,(split /\s+/,$record)[2]} keys %Sample;
	$job_aln=join ",",@job;@job=(); 
	map{$record=`cd $Output/sh/ && qsub -cwd -P $queue -l p=5  -l  num_proc=1 -q mgi.q  -l vf=15g -hold_jid $job_aln step2-$_-Realignment.sh`;push @job,(split /\s+/,$record)[2]} keys %Sample;
	$job_realn=join ",",@job;@job=(); 
	map{$record=`cd $Output/sh/ && qsub -cwd -P $queue -l p=5 -l  num_proc=1 -q mgi.q  -l vf=14g -hold_jid $job_realn step3-Genotyper-chr$_.sh`;push @job,(split /\s+/,$record)[2]}(1..22,'X','Y');
	$job_gt=join ",",@job;@job=(); 
	$record=`cd $Output/sh/ && qsub -cwd -l vf=13g -P $queue  -l  p=5  -l  num_proc=2 -q mgi.q  -hold_jid $job_gt step4-Correct.sh`;
	push @job , (split /\s+/,$record)[2];
	$job_var=join ",",@job;
	$record=`cd $Output/sh/ && qsub -cwd -l vf=3g -P $queue  -l   p=5  -l  num_proc=2 -q mgi.q  -m abes -M $user\@genomics.cn  -hold_jid $job_var step5-Anno.sh`;
}
