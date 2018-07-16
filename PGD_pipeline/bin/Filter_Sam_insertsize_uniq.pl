#!usr/bin/perl -w 
use strict;
use Getopt::Long;
use File::Basename;

my($input,$output,$minI,$maxI);
GetOptions(
	'I=s'    => \$input,
	'O=s'    => \$output,
	'minI=i' => \$minI,
	'maxI=i' => \$maxI,
);
open IN,$input or die $!;
open OU,">$output" or die $!;
my($count_U,$count_R,$count,$count_ins);
while(<IN>){
	if(/^\@/){
		print OU "$_";
		next;
	}
	$count++;
	if(/XT:A:U/){
		$count_U++;	
		my@line=split;
		if(abs($line[8])>=$minI && abs($line[8])<=$maxI){
			print OU "$_";
			$count_ins++;
		}
	}else{
		$count_R++;
	}
}
my$dir=dirname($output);
open STAT,">$dir/uniq_insertsize.stat" or die $!;
print STAT "Total reads\t$count\n";
print STAT "Uniq reads\t$count_U\tPercentage\t",$count_U/$count,"\n";
print STAT "NonUniq reads\t$count_R\tPercentage\t",$count_R/$count,"\n";
print STAT "Left reads\t$count_ins\tPercentage\t",$count_ins/$count,"\n"; 
close STAT;
close OU;
close IN;
