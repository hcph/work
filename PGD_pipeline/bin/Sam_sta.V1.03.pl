#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;
my $usage=<<USAGE;
*********************************************************************************************************
# Function:	Get alingments'statistics messsage (Applied in raw standard SAM,Producted by BWA or else)
# Version:	V1.02
# Last Modify:	2013-1-11 
# Contact:	Yu Wang,<wangyu0309\@gmail.com>
# Usage:	perl $0 [options]
#		
#		-sam	<file>	the sam result with bwa aligment or else
#		-map	<str>	output,contain the mismatch rate or/and mismap ,unmap,SGmap,PEmap,gaps rate etc.Paired-end mapping[PE] or Single-end mapping[SE].
#		-max	<str>	the Max statistics number[NM] of Paired-end(Single-end) reads.default[25M]
#		-insert <str>  	whether product the file, containing the statistics of the library's insertion,[Y] for Yes,[N] for No.default[N]
#		-dupd	<str>	whether product the file, containing the duplication distribution of library,[Y] for Yes,[N] for No.default[N]
#		-out	<dir>	the output directory
#		-label  <str>	Alignment Label[N];
#		-h		<Mes>	help information to screen		
# 
# Example1:	perl $0 -sam file.sam  -map PE -insert Y -dupd Y
#
***********************************************************************************************************
USAGE
my ($sam,$insert,$map,$max,$dupd,$label,$out,%opt);
%opt=GetOptions(
		'sam=s'=>\$sam,
		'insert=s'=>\$insert,
		'map=s'=>\$map,
		'max=s'=>\$max,
		'dupd=s'=>\$dupd,
		'label=s'=>\$label,
		'out=s'=>\$out,
		);
die $usage if (!$sam  or !$map or !$out); 
$insert="N" if(!$insert);
$dupd="N" if(!$dupd);
$label="N" if(!$label);
$max="25M" if(!$max);
$max=1000000*(split /[Mm]/,$max)[0];
my %dup;
my $filename=basename($sam);
$filename=(split /\.sam/,$filename)[0];
mkdir $out unless(-d $out);
$out=~s/\/$//;
if($sam=~/gz$/)
{
	open SAM,"gzip -dc $sam|" or die "Can't open $sam file!\tmark1\n";
}
else
{
	open SAM,$sam or die "Can't open $sam file!\tmark1\n";
}
if($map eq "PE")
{	
	my ($Total,$Total_base,$Umap,$SE,$PE,$Mismap,$Uniq_map);
	my ($Mismat,$Gap,$Bases);
	my %insert;
	while(<SAM>)
	{
		next if ($_=~/^\@/);
		chomp;
		$Total++;
		my $bline1=$_;
		my @bline1=split /\s+/,$bline1;
		if($label eq 'Y')
		{
			$Uniq_map++ if($bline1[12]=~/XT:A:U/);
		}
		elsif($label eq 'N')
		{
			$Uniq_map++ if($bline1[11]=~/XT:A:U/);
		}
		my $bline2=<SAM>;
		chomp $bline2;
		my @bline2=split /\s+/,$bline2;
		if($label eq 'Y')
                {
                        $Uniq_map++ if($bline2[12]=~/XT:A:U/);
                }
                elsif($label eq 'N')
                {
                        $Uniq_map++ if($bline2[11]=~/XT:A:U/);
                }
		$Uniq_map++ if($bline2[11]=~/XT:A:U/);
		$Total_base+=length($bline1[9]);
		if($bline1[2]=~/\*/ && $bline2[2]=~/\*/) ### record unmamp
		{
			$Umap++;
		}
		elsif($bline1[5]=~/\*/ && $bline2[5]!~/\*/) ### record semap
		{
			$SE++;
			my ($XM,$XO);
			foreach(0..@bline2-1)
			{
				$XM=(split /\:/,$bline2[$_])[-1] if($bline2[$_]=~/XM:/);
				$XO=(split /\:/,$bline2[$_])[-1] if($bline2[$_]=~/XO:/);
			}
			$Mismat+=$XM;
			$Gap++ if($XO!=0);
			sam_dup($bline2[2],$bline2[3],$bline2[7]);
		}
		elsif($bline1[5]!~/\*/ && $bline2[5]=~/\*/) ### record semap
		{
			$SE++;
			my ($XM,$XO);
			foreach(0..@bline1-1)
			{
				$XM=(split /\:/,$bline1[$_])[-1] if($bline1[$_]=~/XM:/);
				$XO=(split /\:/,$bline1[$_])[-1] if($bline1[$_]=~/XO:/);
			}
			$Mismat+=$XM;
			$Gap++ if($XO!=0);
			sam_dup($bline1[2],$bline1[3],$bline1[7]);
		}
		elsif($bline1[2]!~/\*/ && $bline2[2]!~/\*/) ### record pemap
		{
			my ($XM1,$XO1,$XM2,$XO2);
			foreach(0..@bline1-1)
			{
				$XM1=(split /\:/,$bline1[$_])[-1] if($bline1[$_]=~/XM:/);
				$XO1=(split /\:/,$bline1[-3])[-1] if($bline1[$_]=~/XO:/);
			}
			foreach(0..@bline2-1)
			{
				$XM2=(split /\:/,$bline2[-4])[-1] if($bline2[$_]=~/XM:/);
				$XO2=(split /\:/,$bline2[-3])[-1] if($bline2[$_]=~/XO:/);
			}
			$Mismat+=($XM1+$XM2);
			$Gap++ if($XO1!=0);
			$Gap++ if($XO2!=0);
			if($bline1[2] ne $bline2[2])
			{
				$Mismap++;
			}else
			{
				$PE++;
				my $insertion=abs($bline1[8]);
				$insert{$insertion}++;
				if($bline1[3]<$bline1[7])
				{
					sam_dup($bline1[2],$bline1[3],$bline1[7]);
				}else{
					sam_dup($bline1[2],$bline1[7],$bline1[3]);
				}

			}
		}
		else
		{
			print "Eorror: Please check your sam format!!!\n";
			last;
		}
		last if($Total==$max);
	}
close SAM;
my %time;
foreach(keys%dup)
{
	$time{$dup{$_}}++;
}
my $dup_read=0;
my $dup_one=0;
foreach(sort {$a<=>$b} keys%time)
{
		$dup_read+=$_*$time{$_};
		$dup_one+=$time{$_} ;
}
my $dup_rate=sprintf("%.3f",(100*($dup_read-$dup_one)/$dup_read));
if($dupd eq "Y")
{
	open DUP,">$out/$filename.dup" or die $!;
	foreach(sort {$a<=>$b} keys%time)
	{
		print DUP "$_\t$time{$_}\n";
	}
}
if($insert eq "Y")
{
	open INSERT,">$out/$filename.insert" or die $!;
	foreach(sort {$a<=>$b} keys%insert)
	{
		print INSERT "$_\t$insert{$_}\n";
	}
}
my $Umap_rate=sprintf ("%.3f",100*$Umap/$Total);
my $SE_rate=sprintf ("%.3f",100*$SE/$Total);
my $PE_rate=sprintf ("%.3f",100*$PE/$Total);
my $Mismap_rate=sprintf ("%.3f",100*$Mismap/$Total);
my $Uniq_rate=sprintf ("%.3f",100*$Uniq_map/($Total*2));
my $Mismat_rate=sprintf ("%.6f",100*$Mismat/$Total_base);
my $Gap_rate=sprintf ("%.6f",100*$Gap/($SE+2*$PE+2*$Mismap));
open MAP,">$out/$filename.map" or die $!;
print MAP "Total_Paired_Reads:\t$Total\n";
print MAP "Total_Paired_Bases:\t$Total_base\n";
print MAP "Unmap_rate:\t$Umap_rate\t###Unmap_rate=100*Unmap_Paired_Reads_num/Total_Paired_Reads\n";
print MAP "SE_rate:\t$SE_rate\t###SE_rate=100*SE_Paired_Reads_num/Total_Paired_Reads\n";
print MAP "PE_rate:\t$PE_rate\t###PE_rate=100*PE_Paired_Reads_num/Total_Paired_Reads\n";
print MAP "Mismap_rate:\t$Mismap_rate\t###Mismap_rate=100*Mismap_Paired_Reads_num/Total_Paired_Reads\n";
print MAP "Uniq_rate:\t$Uniq_rate\t###Uniq_rate=100*Uniq_Reads_num/Total_Reads\n";
print MAP "Mismat_rate:\t$Mismat_rate\t###Mismat_rate=100*Mismatch_Bases/Total_Bases\n";
print MAP "Gap_rate:\t$Gap_rate\t###Gap_rate=The_Gap_reads/Total_map_Reads\n";
print MAP "Duplication_rate:\t$dup_rate\n";
}
if($map eq "SE")
{
	my ($Total,$SE,$Unmap,$Total_base,$Gap,$Uniq_map,$Mismat,%dup);
	while(<SAM>)	
	{
		next if($_=~/^\@/);
		chomp;
		$Total++;
		my @line=split /\s+/,$_;
		if($label eq 'Y')
                {
                        $Uniq_map++ if($line[12]=~/XT:A:U/);
                }
                elsif($label eq 'N')
                {
                        $Uniq_map++ if($line[11]=~/XT:A:U/);
                }
		$Total_base+=length($line[9]);
		if($line[2]=~/\*/) ###  record unmap
		{
			$Unmap++;
			next;		
		}
		else ### record semap
		{
			$SE++;
			my ($XM,$XO);
			foreach(0..@line-1)
			{
				$XM=(split /\:/, $line[$_])[-1] if($line[$_]=~/XM:/);
				$XO=(split /\:/, $line[$_])[-1] if($line[$_]=~/XO:/);
			}
			$Mismat+=$XM;
			$Gap++ if($XO!=0);
			$dup{$line[2]."-".$line[3]}++;
		}
		last if($Total==$max);

	}
	close SAM;
	my %time;
	foreach(keys%dup)
	{
		        $time{$dup{$_}}++;
	}
	my $dup_read=0;
	my $dup_one=0;
	foreach(sort {$a<=>$b} keys%time)
	{
		                $dup_read+=$_*$time{$_};
						$dup_one+=$time{$_} ;
	}
	my $dup_rate=sprintf("%.3f",(100*($dup_read-$dup_one)/$dup_read));
	my $Unmap_rate=sprintf("%.3f",100*$Unmap/$Total);
	my $SE_rate=sprintf("%.3f",100*$SE/$Total);
	my $Uniq_rate=sprintf ("%.3f",100*$Uniq_map/$Total);
	my $Mismat_rate=sprintf("%.6f",100*$Mismat/$Total_base);
	my $Gap_rate=sprintf("%.6f",100*$Gap/$SE);
	open MAP,">$out/$filename.map" or die $!;
	print MAP "Total_Singled_Reads:\t$Total\n";
	print MAP "Total_Singled_Bases:\t$Total_base\n";
	print MAP "Unmap_rate:\t$Unmap_rate\t###Unmap_rate=100*Unmap_Singled_Reads_num/Total_Singled_Reads\n";
	print MAP "SE_rate:\t$SE_rate\t###SE_rate=100*SE_Reads_num/Total_Single_Reads\n";
	print MAP "Uniq_rate:\t$Uniq_rate\t###Uniq_rate=100*Uniq_Reads_num/Total_Reads\n";
	print MAP "Mismat_rate:\t$Mismat_rate\t###Mismat_rate=100*Mismatch_Bases/Total_Bases\n";
	print MAP "Gap_rate:\t$Gap_rate\t###Gap_rate=The_Gap_reads/Total_map_Reads\n";
	print MAP "Duplication_rate:\t$dup_rate\n";
	if($dupd eq "Y")
	{
		        open DUP,">$out/$filename.dup" or die $!;
		        foreach(sort {$a<=>$b} keys%time)
		        {
	                print DUP "$_\t$time{$_}\n";
		        }
	}


}
sub sam_dup{
				my ($chr,$start,$end)=@_[0,1,2];
				my $bmark=$chr."-".$start."-".$end;
				$dup{$bmark}++;
			}

