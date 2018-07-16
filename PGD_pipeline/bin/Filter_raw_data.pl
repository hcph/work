#!/usr/bin/perl 
use strict;
use Getopt::Long;
use File::Basename;
my $usage=<<USAGE;
**************************************************************************************************************************
# Function:	Filter Solexa/Sanger system raw data	
# Version:	V1.02
# Modify1:	2011-11-28
# Modify2:	2011-12-1 1)--> add Rscript statistic,the rate of validated data. 
# Last modify:	2011-12-9 1)-->deal with parameters bugs
# Contact:	Aillen Wang,<wangyu\@genomics.org.cn>,BGI-SZ
# Usage:	perl $0	[options]
#		-Fq1	<file>	raw_read1.fq(or fq.gz) pathway
#		-Res	<dir>	result output directary,contained  clean_read.fq, filter mes and so on
#		-N	<int>	filter the total number(>=x) of N in the PE read, default [5]
#		-Q_cut	<int>	filter the PE reads whose the average of base quality <=x, default [20]
#		-Q_system	1,stand for Solexa,which is the system of quality based on the Value [64], ASCII [@] 
#				2,stand for Sanger,which is the system of quality based on the Value [33], ASCII [i]
#		-Ada	<file>	raw_read1.adapter.list(or list.gz) pathway
#		-Rsta	<str>	Draw the pie chart,the rate of validated data,[Yes] for do,[No] for don't,default[No]
#		-h		output help information to screen
# Example1:perl $0 -Fq1 *_1.fq.gz -Res result_dir -N 5 -Q_cut 20 -Q_system 1 -Ada *1.adapter.list.gz
# Example2:perl $0 -Fq1 *_1.fq    -Res result_dir -N 5 -Q_cut 20 -Q_system 2 
**************************************************************************************************************************
USAGE
my ($Fq1,$Res,$N,$Q,$system,$Ada,$Rsta,$h,%opt);
%opt=GetOptions(
		'Fq1=s'=>\$Fq1,
		'Res=s'=>\$Res,
		'N=i'=>\$N,
		'Q_cut=i'=>\$Q,
		'Q_system=i'=>\$system,
		'Ada=s'=>\$Ada,
		'Rsta=s'=>\$Rsta,
		'h'=>\$h
		);
die $usage if (!$Fq1 or !$Res or !$system or $h);
$N=5 if(!$N);
$Q=20 if(!$Q);
$Rsta="No" if(!$Rsta);
my ($base_value,$adapter1,$adapter2,$ID);
$base_value=64 if($system==1);$base_value=33 if($system==2);
my $fq1=$Fq1;
my $fq2=$Fq1;
$ID=(split /\//,$fq1)[-1];
$ID=(split /\_1\./,$ID)[-2];
my %mark;
if($Ada)
{
	$adapter1=$Ada;
	$adapter2=$Ada;
	$adapter2=~s/1\.adapter\.list/2\.adapter\.list/ if ($Ada=~/\.list$/); $adapter2=~s/1\.adapter\.list\.gz/2\.adapter\.list\.gz/ if ($Ada=~/\.gz$/);
	if($adapter1=~/gz$/)
	{
		open IN1, "gzip -dc $adapter1|" or die $!;
	}
	else
	{
		open IN1,$adapter1 or die $!;
	}
	if($adapter2=~/gz$/)
	{
		open IN2, "gzip -dc $adapter2|" or die $!;
	}
	else
	{
		open IN2,$adapter2 or die $!;
	}
	while (<IN1>)
	{
		next if ($_=~/^#reads_id/);
		chomp;
		my ($name,$start)=(split /\s+/,$_)[0,2];
		my $lame="\@".$name;
		my $Name=(split /\//,$lame)[0];
		$mark{$Name}=$start;
	}
	close IN1;
	while (<IN2>)
	{
		next if ($_=~/^#reads_id/);
		chomp;
		my ($name,$start)=(split /\s+/,$_)[0,2];
		my $lame="\@".$name;
		my $Name=(split /\//,$lame)[0];
		$mark{$Name}=$start;
	}
	close IN2;
}
$fq2=~s/_1\.fq\.gz/_2\.fq\.gz/ if($fq2=~/\.gz$/);$fq2=~s/_1\.fq/_2\.fq/ if($fq2=~/\.fq$/);
if($fq1=~/gz$/)
{
	open IN3,"gzip -dc $fq1|" or die $! ;
}
else
{
	open IN3,$fq1 or die $!;
}
if($fq2=~/gz$/)
{
	open IN4,"gzip -dc $fq2|" or die $!;
}
else
{
	open IN4,$fq2 or die $!;
}
open OU1,"| gzip >$Res/$ID\_1.clean.fq.gz" or die $!;
open OU2,"| gzip >$Res/$ID\_2.clean.fq.gz" or die $!;
open OU3,">$Res/$ID\.sta";
my ($pre_r1a,$pre_r1t,$pre_r1c,$pre_r1g,$r1a,$r1t,$r1c,$r1g,$pre_r2a,$pre_r2t,$pre_r2c,$pre_r2g,$r2a,$r2t,$r2c,$r2g,$pre_r1q20,$pre_r1q25,$pre_r1q30,$pre_r2q20,$pre_r2q25,$pre_r2q30,$r1q20,$r1q25,$r1q30,$r2q20,$r2q25,$r2q30,$adap,$pre_total,$total,$pre_N_1,$pre_N_2,$pre_N_3,$N_1,$N_2,$N_3);
while (<IN3>)
{
	chomp $_;
	my $read1_name=$_;
	my $mark1=(split /\//,$read1_name)[0];
	my $read2_name=<IN4>;
	chomp $read2_name;
	my $mark2=(split /\//,$read2_name)[0];
	my $read1_base_line=<IN3>;
	my $read2_base_line=<IN4>;
	my $read1_ori=<IN3>;
	my $read2_ori=<IN4>;
	my $read1_quality=<IN3>;
	my $read2_quality=<IN4>;
	chomp $read1_base_line;
	chomp $read2_base_line;
	chomp $read1_quality;
	chomp $read2_quality;
	chomp $read1_ori;
	chomp $read2_ori;	
	my @Qua1= unpack ("C*",$read1_quality);
	my @Qua2= unpack ("C*",$read2_quality);
	my $pre_r1A=$read1_base_line=~tr/A/A/;
	my $pre_r1T=$read1_base_line=~tr/T/T/;
	my $pre_r1C=$read1_base_line=~tr/C/C/;
	my $pre_r1G=$read1_base_line=~tr/G/G/;
	$pre_r1a+=$pre_r1A;
	$pre_r1t+=$pre_r1T;
	$pre_r1c+=$pre_r1C;
	$pre_r1g+=$pre_r1G;
	my $pre_r2A=$read2_base_line=~tr/A/A/;
	my $pre_r2T=$read2_base_line=~tr/T/T/;
	my $pre_r2C=$read2_base_line=~tr/C/C/;
	my $pre_r2G=$read2_base_line=~tr/G/G/;
	$pre_r2a+=$pre_r2A;
	$pre_r2t+=$pre_r2T;
	$pre_r2c+=$pre_r2C;
	$pre_r2g+=$pre_r2G;
	my $N1=$read1_base_line=~tr/N/N/;
	my $N2=$read2_base_line=~tr/N/N/;
	my $M=$N1+$N2;
	 if($M==1)
	 {
	 	$pre_N_1++;
	 }
	 elsif($M==2)
	 {
	 	$pre_N_2++;
	 }
	 elsif($M==3)
	 {
	 	$pre_N_3++;
	 }
	my $Qua1_sum=0;
	my $Qua2_sum=0;
	my $len1=length($read1_quality);
	my $len2=length($read2_quality);
	for(0..@Qua1-1)
	{
		$Qua1_sum+=($Qua1[$_]-$base_value);
	}
	my $Qua1_ave=$Qua1_sum/$len1;
	$pre_r1q20++  if ($Qua1_ave>=20);
	$pre_r1q25++  if ($Qua1_ave>=25);
	$pre_r1q30++  if ($Qua1_ave>=30);
	for(0..@Qua2-1)
	{
		$Qua2_sum+=($Qua2[$_]-$base_value);		
	}
	my $Qua2_ave=$Qua2_sum/$len2;
	$pre_r2q20++ if ($Qua2_ave>=20);
	$pre_r2q25++ if ($Qua2_ave>=25);
	$pre_r2q30++ if ($Qua2_ave>=30);
	$pre_total++;
	if($Ada)
	{
		if(exists $mark{$mark1} or exists $mark{$mark2})
		{
			$adap++;
			next;
		}
		else
		{
			if($Qua1_ave>=$Q && $Qua2_ave>=$Q  && $M<=$N)
			{
				$total++;
				$r1a+=$pre_r1A;
				$r1t+=$pre_r1T;
				$r1c+=$pre_r1C;
				$r1g+=$pre_r1G;
				$r2a+=$pre_r2A;
				$r2t+=$pre_r2T;
				$r2c+=$pre_r2C;
				$r2g+=$pre_r2G;
				$r1q20++  if ($Qua1_ave>=20);
				$r1q25++  if ($Qua1_ave>=25);
				$r1q30++  if ($Qua1_ave>=30);
				$r2q20++ if ($Qua2_ave>=20);
				$r2q25++ if ($Qua2_ave>=25);
				$r2q30++ if ($Qua2_ave>=30);
	 			if($M==1)
	 			{
	 				$N_1++;
	 			}
	 			elsif($M==2)
	 			{
	 				$N_2++;
	 			}
	 			elsif($M==3)
	 			{
	 				$N_3++;
	 			}
				print OU1 "$read1_name\n$read1_base_line\n$read1_ori\n$read1_quality\n";
				print OU2 "$read2_name\n$read2_base_line\n$read2_ori\n$read2_quality\n";
			}
		}
	}
	else
	{
			if($Qua1_ave>=$Q && $Qua2_ave>=$Q  && $M<=$N)
			{
				$total++;
				$r1a+=$pre_r1A;
				$r1t+=$pre_r1T;
				$r1c+=$pre_r1C;
				$r1g+=$pre_r1G;
				$r2a+=$pre_r2A;
				$r2t+=$pre_r2T;
				$r2c+=$pre_r2C;
				$r2g+=$pre_r2G;
				$r1q20++  if ($Qua1_ave>=20);
				$r1q25++  if ($Qua1_ave>=25);
				$r1q30++  if ($Qua1_ave>=30);
				$r2q20++ if ($Qua2_ave>=20);
				$r2q25++ if ($Qua2_ave>=25);
				$r2q30++ if ($Qua2_ave>=30);
	 			if($M==1)
	 			{
	 				$N_1++;
	 			}
	 			elsif($M==2)
	 			{
	 				$N_2++;
	 			}
	 			elsif($M==3)
	 			{
	 				$N_3++;
	 			}
				print OU1 "$read1_name\n$read1_base_line\n$read1_ori\n$read1_quality\n";
				print OU2 "$read2_name\n$read2_base_line\n$read2_ori\n$read2_quality\n";
			}
		
	}
}
close IN3;
close IN4;
close OU1;
close OU2;
my $pre_base=$pre_r1a+$pre_r1t+$pre_r1g+$pre_r1c+$pre_r2a+$pre_r2t+$pre_r2g+$pre_r2c;
my $N1_1=100*$pre_N_1/$pre_total;
my $N2_1=100*$pre_N_2/$pre_total;
my $N3_1=100*$pre_N_3/$pre_total;
my $R1AT_1=100*($pre_r1a-$pre_r1t)/($pre_r1a+$pre_r1t+$pre_r1g+$pre_r1c);
my $R1GC_1=100*($pre_r1g-$pre_r1c)/($pre_r1a+$pre_r1t+$pre_r1g+$pre_r1c);
my $R1per_1=100*($pre_r1g+$pre_r1c)/($pre_r1a+$pre_r1t+$pre_r1g+$pre_r1c);
my $R1q20_1=100*$pre_r1q20/$pre_total;
my $R1q25_1=100*$pre_r1q25/$pre_total;
my $R1q30_1=100*$pre_r1q30/$pre_total;
my $R2AT_1=100*($pre_r2a-$pre_r2t)/($pre_r2a+$pre_r2t+$pre_r2g+$pre_r2c);
my $R2GC_1=100*($pre_r2g-$pre_r2c)/($pre_r2a+$pre_r2t+$pre_r2g+$pre_r2c);
my $R2per_1=100*($pre_r2g+$pre_r2c)/($pre_r2a+$pre_r2t+$pre_r2g+$pre_r2c);
my $R2q20_1=100*$pre_r2q20/$pre_total;
my $R2q25_1=100*$pre_r2q25/$pre_total;
my $R2q30_1=100*$pre_r2q30/$pre_total;
my $base=$r1a+$r1t+$r1g+$r1c+$r2a+$r2t+$r2g+$r2c;
my $N1_2=100*$N_1/$total;
my $N2_2=100*$N_2/$total;
my $N3_2=100*$N_3/$total;
my $R1AT_2=100*($r1a-$r1t)/($r1a+$r1t+$r1g+$r1c);
my $R1GC_2=100*($r1g-$r1c)/($r1a+$r1t+$r1g+$r1c);
my $R1per_2=100*($r1g+$r1c)/($r1a+$r1t+$r1g+$r1c);
my $R1q20_2=100*$r1q20/$total;
my $R1q25_2=100*$r1q25/$total;
my $R1q30_2=100*$r1q30/$total;
my $R2AT_2=100*($r2a-$r2t)/($r2a+$r2t+$r2g+$r2c);
my $R2GC_2=100*($r2g-$r2c)/($r2a+$r2t+$r2g+$r2c);
my $R2per_2=100*($r2g+$r2c)/($r2a+$r2t+$r2g+$r2c);
my $R2q20_2=100*$r2q20/$total;
my $R2q25_2=100*$r2q25/$total;
my $R2q30_2=100*$r2q30/$total;
my $filter_rate=100-(100*$total/$pre_total);
my $adapter_per=0;	
$adapter_per=100*$adap/$pre_total if($Ada);
my $clean_per=100-$filter_rate;

if($Rsta =~ /Yes/)
{	
	`echo 'Dirty_data_rate	$filter_rate' >$Res/temp.txt`;
	`echo 'Clean_data_rate  $clean_per' >>$Res/temp.txt`;
my $Rscript=<<RS;
png("$ID.pie.percent.png",width=800,height=600)
c=read.table("$Res/temp.txt")
ratio=sprintf("%0.2f",c[,2])
ratio=paste(ratio,"%",sep="")
label=paste(c[,1],ratio,sep="\n")
pie(c[,2],col=c("red","green"),labels=label,font=2)
legend("topleft",legend=c("Dirty_data_rate","Clean_data_rate"),pch=16,col=c("red","green"))
RS
	open O1,">$Res/Rscript.R" || die "Can't open the output file:$!\n";
	print O1 "$Rscript";
	close O1;
	system("/opt/blc/genome/biosoft/R/bin/Rscript $Res/Rscript.R");
	system("rm $Res/Rscript.R $Res/temp.txt");
	system("mv $ID.pie.percent.png  $Res");
}

printf OU3 "Raw_data_statistics\|\tBases(bp):$pre_base\tReads(PE_num):$pre_total\tAdapter_contaminate_rate(\%):%.2f\n",$adapter_per;
printf OU3 "Raw_N_statistics\|\t(1N)\%:%.2f\t(2N)\%:%.2f\t(3N)\%:%.2f\n",$N1_1,$N2_1,$N3_1;
printf OU3 "Raw_Read1_base_statistics\|\t(A-T)\%:%.2f\t(G-C)\%:%.2f\t(GC)\%:%.2f\n",$R1AT_1,$R1GC_1,$R1per_1;
printf OU3 "Raw_Read2_base_statistics\|\t(A-T)\%:%.2f\t(G-C)%:%.2f\t(GC)\%:%.2f\n",$R2AT_1,$R2GC_1,$R2per_1;
printf OU3 "Raw_Read1_qual_statistics\|\t(q20)\%:%.2f\t(q25)\%:%.2f\t(q30)\%:%.2f\n",$R1q20_1,$R1q25_1,$R1q30_1;
printf OU3 "Raw_Read2_qual_statistics\|\t(q20)\%:%.2f\t(q25)\%:%.2f\t(q30)\%:%.2f\n",$R2q20_1,$R2q25_1,$R2q30_1;
printf OU3 "#########################################################################################################\n";
printf OU3 "# Filter_condition:| 1.filter the PE reads which have been contaminated by adapter.\n";
printf OU3 "#                  | 2.filter the total number(>=$N) of N in the PE reads.\n";
printf OU3 "#                  | 3.filter the PE reads whose the average of base quality <=$Q.\n";
printf OU3 "#########################################################################################################\n";
printf OU3 "Clean_data_statistics\|\tBases(bp):$base\tReads(PE_num):$total\tFilter_rate(\%):%.2f\n",$filter_rate;
printf OU3 "Clean_N_statistics\|\t(1N)\%:%.2f\t(2N)\%:%.2f\t(3N)\%:%.2f\n",$N1_2,$N2_2,$N3_2;
printf OU3 "Clean_Read1_base_statistics\|\t(A-T)\%:%.2f\t(G-C)\%:%.2f\t(GC)\%:%.2f\n",$R1AT_2,$R1GC_2,$R1per_2;
printf OU3 "Clean_Read2_base_statistics\|\t(A-T)\%:%.2f\t(G-C)\%:%.2f\t(GC)\%:%.2f\n",$R2AT_2,$R2GC_2,$R2per_2;
printf OU3 "Clean_Read1_qual_statistics\|\t(q20)\%:%.2f\t(q25)\%:%.2f\t(q30)\%:%.2f\n",$R1q20_2,$R1q25_2,$R1q30_2;
printf OU3 "Clean_Read2_qual_statistics\|\t(q20)\%:%.2f\t(q25)\%:%.2f\t(q30)\%:%.2f\n",$R2q20_2,$R2q25_2,$R2q30_2;
close OU3;
