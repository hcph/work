#!usr/bin/perl -w 
use strict;
use Getopt::Long;
use Test::More;

my $usage=<<USAGE;
        Author:     huangyile
        Version:    V3
        Usage:      -I      <File>  GATK vcf file;
		    -DE     <INT>   The mininum of Embryos' depth, Default[10];
		    -DT     <INT>   The mininum of Trios depth, Default[10];
                    -ID     <INT>   Whether use indel or not [Y|N] default [Y];
		    -C      <STR>   whether need to correct gatk default mutation[All|gDNA|None|Embryo], Default[None];
                    -O      <Dir>   Output path;
                    -T      <INT>   threshold of hom/ref mutation, Default[0.1];
		    -R      <File>  All capture genes' region
		    -H	    <help>  Print help message;
        Example:    perl $0 -I VCF -O OUT 
		    perl $0 -I VCF -O OUT -DE int -DT int -C [all|gDNA|none] -T int -G str -ID Y -R target.gene 

USAGE

my($Output,$Input,$EmbryoNnm,$gDNAMinDepth,$Correct,$Threshold,$EmbryoMinDepth,$help,$FilInDel,$target,$core);
GetOptions
(
        'DE=i'  => \$EmbryoMinDepth,
        'DT=i'  => \$gDNAMinDepth,
        'C=s'   => \$Correct,
        'O=s'   => \$Output,
	'I=s'   => \$Input,
        'T=f'   => \$Threshold,
	'R=s'   => \$target,
	'ID=s'  => \$FilInDel,
	'H'     => \$help
);
die $usage if($help);

my(%correct,%Hap_support_num,%sample,@sample,%region,%informative);
$EmbryoMinDepth||=10,$gDNAMinDepth||=10,$Threshold||=0.8,$FilInDel||="Y",$Correct||='None';
my$EmbryoNum=0,my@gDNA,my@Embryo,my%ROW,my@snp;

`mkdir $Output` unless(-d "$Output");
open VCF,"<$Input" or die "Can't open VCF file!!!  Pleast check it!\n$usage";
open HAP,">$Output\/Family.hap.xls" or die "Can't wirte to hap file! Pleast check it!\n$usage";
open LOW,">$Output\/Family.lowdepth.xls" or die "Can't wirte to lowdepth file! Pleast check it!\n$usage";
open NOV,">$Output\/Family.denove.xls" or die "Can't wirte to denove file! Pleast check it!\n$usage";
open SCR,">$Output\/Family.score.xls" or die "Can't wirte to score file! Pleast check it!\n$usage";
open SNP,">$Output\/Family.snp.xls" or die "Can't wirte to score file! Pleast check it!\n$usage";

#*********************************************************************************************************
#                                     Reading VCF file, and Correct SNP                                  *
#*********************************************************************************************************
open TAR,$target or die $!;
while(<TAR>){
	chomp;
        $region{(split)[0]}=[(split)[1..3]];
}
close TAR;

my$isPaternal=0,my$isMaternal=0;
START:while(<VCF>){
	chomp;
        my $lines=$_;
	next if(/\#\#/ || /LowQual/);
	my@line=split;
	if($line[0] =~ /\#CHROM/){
                map{push @snp,$line[$_]}(9..$#line);
                print SNP "Gene\tCHROM\tPOS\t";
                map{print SNP "$_\t"}(@snp);
                print SNP "\n";
		map{if($_=~/Embryo/){$EmbryoNum++}} @line;
		if(/Proband/){$core = 'Y'}else{$core = 'N'};
		if($core eq 'Y'){
			@sample=('Father','Mother','Proband');
	 		map{if($_ =~ /Paternal/i){push @sample,'Paternal';$isPaternal=1};if($_=~/Maternal/i){push @sample,'Maternal';$isMaternal=1}}@line;
                        foreach $_(@snp){
                        next if $_=~ /Embryo/;      	
                        push @gDNA,$_;}
                        print "@gDNA\n";
			map{push @Embryo,'Embryo'.$_;push @sample,'Embryo'.$_}(1..$EmbryoNum);
			print HAP "Gene\tchr\tPos\tF1\tF2\tM1\tM2\t";
			print HAP "Mut\t" if($isPaternal==1||$isMaternal==1);
			map{print HAP "Em$_\t"}(1..$EmbryoNum);
			map{print HAP "$_\t"}(@sample);
			map{print HAP "$_\t"}(@sample);
			print HAP "\n";
			print LOW "Chr\tPos\t";
			map{print LOW "$_\t"}(@sample);
			map{print LOW "$_\t"}(@sample);
			print LOW "\n";
			print NOV "Chr\tPos\t";
			map{print NOV "$_\t"}(@sample);
			map{print NOV "$_\t"}(@sample);
			print NOV "\n";
			print SCR "Embryo\tGene\tF1\tF2\tM1\tM2\tPaternal\tMaternal\n";
		}elsif($core eq 'N'){
=cut                   
			if($_ =~  /Paternal/i){@sample=('Father','Paternal','Mother');$isPaternal=1}
			elsif($_=~/Maternal/i){@sample=('Mother','Maternal','Father');$isMaternal=1}
=cut
                        if($_=~/Maternal/i){@sample=('Mother','Maternal','Father');$isMaternal=1}
			@gDNA=@sample;
                        map{push @Embryo,'Embryo'.$_;push @sample,'Embryo'.$_}(1..$EmbryoNum);
                        print HAP "Gene\tchr\tPos\tF1\tF2\t" if($isPaternal==1);
                        print HAP "Gene\tchr\tPos\tM1\tM2\t" if($isMaternal==1);
                        map{print HAP "Em$_\t"}(1..$EmbryoNum);
                        map{print HAP "$_\t"}(@sample);
                        map{print HAP "$_\t"}(@sample);
                        print HAP "\n";
                        print LOW "Chr\tPos\t";
                        map{print LOW "$_\t"}(@sample);
                        map{print LOW "$_\t"}(@sample);
                        print LOW "\n";
			print NOV "Chr\tPos\t";
			map{print NOV "$_\t"}(@sample);
			map{print NOV "$_\t"}(@sample);
			print NOV "\n";
			print SCR "Embryo\tGene\tF1\tF2\tM1\tM2\n";
		}
		map{$ROW{$line[$_]}=$_}(9..$#line);
		
		#map{$correct{$_}=0}@sample;
		#if($Correct=~/all/i){map{$correct{$_}=1}@sample}
		#elsif($Correct=~/gDNA/i){
                #map{$correct{$_}=1}@gDNA;
		#elsif($Correct=~/embryo/i){map{$correct{$_}=1}@Embryo}
		next START;
	}
	
	$line[4]=~/,/?next:(my$ref=$line[3],my$mut=$line[4]);
        next if($FilInDel =~ /Y/i && (length$line[3]>1||length$line[4]>1));
	
	my%DP=(),my%GT=();
	foreach my$val(keys %ROW){
		$DP{$val}=($line[$ROW{$val}]!~/:/||$line[$ROW{$val}]=~/:\.:/)?[0,0]:[split /,/,(split /:/,$line[$ROW{$val}])[1]];
		$GT{$val}= $line[$ROW{$val}]=~/:/?[split /\//,(split /:/,$line[$ROW{$val}])[0]]:[-1,-1];
		if(not defined $DP{$val}[1]){$DP{$val}=[0,0]}
        }
	
	foreach my$val(@gDNA){
		my$sign=0;
		$sign=1 if($DP{$val}[0]+$DP{$val}[1]<$gDNAMinDepth && $GT{$val}[0]==$GT{$val}[1]     );
                if ($GT{$val}[0] ne $GT{$val}[1]){
                             if ($DP{$val}[0]==1){
                                   $GT{$val}[0]=$GT{$val}[1] if ($DP{$val}[1]>40);
                                   $sign=1 if ($DP{$val}[1]<40);  
                                 }elsif($DP{$val}[1]==1){
                                   $GT{$val}[1]=$GT{$val}[0] if ($DP{$val}[0]>40);
                                   $sign=1 if ($DP{$val}[0]<40);}
                             elsif($DP{$val}[0]==2){
                                   $GT{$val}[0]=$GT{$val}[1] if ($DP{$val}[1]/$DP{$val}[0]>50);
                                   $sign=1 if ($DP{$val}[1]/$DP{$val}[0]>5 && $DP{$val}[1]/$DP{$val}[0]<50 );}elsif($DP{$val}[1]==2){
                                   $GT{$val}[1]=$GT{$val}[0] if ($DP{$val}[0]/$DP{$val}[1]>50);
                                   $sign=1 if ($DP{$val}[0]/$DP{$val}[1]>5 && $DP{$val}[0]/$DP{$val}[1]<50 );}elsif($DP{$val}[0]==3){
                                   $GT{$val}[0]=$GT{$val}[1] if ($DP{$val}[1]/$DP{$val}[0]>60);
                                   $sign=1 if ($DP{$val}[1]/$DP{$val}[0]>8 && $DP{$val}[1]/$DP{$val}[0]<60 );}elsif($DP{$val}[1]==3){
                                   $GT{$val}[1]=$GT{$val}[0] if ($DP{$val}[0]/$DP{$val}[1]>60);
                                   $sign=1 if ($DP{$val}[0]/$DP{$val}[1]>8 && $DP{$val}[0]/$DP{$val}[1]<60 );}elsif($DP{$val}[0]>3){
                                   $GT{$val}[0]=$GT{$val}[1] if ($DP{$val}[1]/$DP{$val}[0]>80);
                                   $sign=1 if ($DP{$val}[1]/$DP{$val}[0]>12 && $DP{$val}[1]/$DP{$val}[0]<80 );}elsif($DP{$val}[1]>3){
                                   $GT{$val}[1]=$GT{$val}[0] if ($DP{$val}[0]/$DP{$val}[1]>80);
                                   $sign=1 if ($DP{$val}[0]/$DP{$val}[1]>12 && $DP{$val}[0]/$DP{$val}[1]<80 );}
};   
                             
     
		#$sign=1 if($GT{$val}[0] ne $GT{$val}[1] && ${$DP{$val}}[0]>0 && ${$DP{$val}}[1]>0&&($DP{$val}[0]<2||$DP{$val}[1]<2));
		$sign=1 if($line[0]=~/chrX/&& $val eq 'Father' && $DP{$val}[0]>2 && $DP{$val}[1]>2);
		if($sign==1){
			print LOW "$line[0]\t$line[1]\t";
               		map{print LOW "$GT{$_}[0]\/$GT{$_}[1]\t"}(@sample);
			map{print LOW $DP{$_}[0]+$DP{$_}[1].":$DP{$_}[0]:$DP{$_}[1]\t"}(@sample);
               		print LOW "\n";
			next START;
		}
	}
            
        my$gene="none";
        map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$gene=$_}}keys %region;
        print SNP "$gene\t$line[0]\t$line[1]\t$line[3]\t$line[4]\t";
        map{print SNP ".$GT{$_}[0]\/$GT{$_}[1]\t"}(@snp);
        print SNP "\n";	
      
        map{@{$GT{$_}}[0,1]=&GT_correct($ref,$mut,$GT{$_}[0],$GT{$_}[1],$DP{$_}[0],$DP{$_}[1])}@snp;

      
	if($core eq 'Y'){
		my%base=(),my%Hap=(),my%score=(),my$Linkage_Provide="-";
        	map{$base{$_}++}(@{$GT{'Father'}},@{$GT{'Mother'}});
		unless(exists $base{$GT{'Proband'}[0]} && exists $base{$GT{'Proband'}[1]}){
			print NOV "$line[0]\t$line[1]\t";
			map{print NOV "$GT{$_}[0]/$GT{$_}[1]\t"}(@sample);
			map{print NOV $DP{$_}[0]+$DP{$_}[1].":$DP{$_}[0]:$DP{$_}[1]\t"}(@sample);
			print NOV "\n";
			next START;
		}
		#print $DP{'Father'}[0],"\t",$DP{'Father'}[1],"\t", $DP{'Mother'}[0],"\t",$DP{'Mother'}[1],"\n";
		if($GT{'Father'}[0] eq $GT{'Father'}[1] && $GT{'Mother'}[0] ne $GT{'Mother'}[1]){
                	my$tribase=""; map {$tribase = $_ if ($base{$_} == 3)} keys %base;
                	my$uniqbase="";$uniqbase=$GT{'Mother'}[0] eq $tribase ? $GT{'Mother'}[1] : $GT{'Mother'}[0];
                	$Hap{'F1'}=$Hap{'F2'}=$tribase;
                	if($GT{'Proband'}[0] eq $uniqbase || $GT{'Proband'}[1] eq $uniqbase){
                        	$Hap{'M1'} = $uniqbase,$Hap{'M2'} = $tribase;
				if($isMaternal==1 && $GT{'Maternal'}[0] eq $GT{'Maternal'}[1] && $GT{'Maternal'}[1] ne $GT{'Father'}[0]){
					if   ($uniqbase eq $GT{'Maternal'}[0]){$Linkage_Provide='M1'}
					elsif($uniqbase ne $GT{'Maternal'}[0]){$Linkage_Provide='M2'}
				}
				map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$informative{$_}{'M1'}++}
				}keys %region;
				map{my$num=$_;
                                my $eb=0;
                                $eb=&EB_correct($uniqbase,$GT{$num}[0],$GT{$num}[1],$DP{$num}[0],$DP{$num}[1]);
				if($eb==1){$Hap{$num}='M1';
				map{if($line[0]eq$region{$_}[0]&& $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$Hap_support_num{$_}{$num}{'M1'}++}
				}keys %region}else{$Hap{$num}='-'}}(@Embryo);
			}elsif($GT{'Proband'}[0] eq $tribase && $GT{'Proband'}[1] eq $tribase){
				$Hap{'M2'} = $uniqbase,$Hap{'M1'} = $tribase;
                                if($isMaternal==1 && $GT{'Maternal'}[0] eq $GT{'Maternal'}[1] && $GT{'Maternal'}[1] ne $GT{'Father'}[0]){
                                        if   ($uniqbase eq $GT{'Maternal'}[0]){$Linkage_Provide='M2'}
                                        elsif($uniqbase ne $GT{'Maternal'}[0]){$Linkage_Provide='M1'}
                                }
				map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$informative{$_}{'M2'}++}
                                }keys %region;
                        	map{my$num=$_;
			        my $eb=0;
                                $eb=&EB_correct($uniqbase,$GT{$num}[0],$GT{$num}[1],$DP{$num}[0],$DP{$num}[1]);
                                if($eb==1){$Hap{$num}='M2';
				map{if($line[0]eq$region{$_}[0]&& $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$Hap_support_num{$_}{$num}{'M2'}++}
				}keys %region}else{$Hap{$num}='-'}}(@Embryo);
			}
        	}	
	#Paternal hap;
		elsif($GT{'Father'}[0] ne $GT{'Father'}[1] && $GT{'Mother'}[0] eq $GT{'Mother'}[1]){
			my$tribase=""; map {$tribase = $_ if ($base{$_} == 3)} keys %base;
                	my$uniqbase="";$uniqbase=$GT{'Father'}[0] eq $tribase ? $GT{'Father'}[1] : $GT{'Father'}[0];
                	$Hap{'M1'}=$Hap{'M2'} = $tribase;
                	if($GT{'Proband'}[0] eq $uniqbase || $GT{'Proband'}[1] eq $uniqbase){
                	        $Hap{'F1'} = $uniqbase,$Hap{'F2'} = $tribase;
                                if($isPaternal==1 && $GT{'Paternal'}[0] eq $GT{'Paternal'}[1] && $GT{'Paternal'}[0] ne $GT{'Mother'}[0]){
                                        if   ($uniqbase eq $GT{'Paternal'}[0]){$Linkage_Provide='F1'}
                                        elsif($uniqbase ne $GT{'Paternal'}[0]){$Linkage_Provide='F2'}
                                }
				map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$informative{$_}{'F1'}++}
				}keys %region;	
				map{my$num=$_;
				my $eb=0;
                                $eb=&EB_correct($uniqbase,$GT{$num}[0],$GT{$num}[1],$DP{$num}[0],$DP{$num}[1]);
                                if($eb==1){$Hap{$num}='F1';
				map{if($line[0]eq$region{$_}[0]&& $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$Hap_support_num{$_}{$num}{'F1'}++}
				}keys %region}else{$Hap{$num}='-'}}(@Embryo);
                	        map{my$num=$_; my $eb=0;
                                 $eb=&EB_correct($uniqbase,$GT{$num}[0],$GT{$num}[1],$DP{$num}[0],$DP{$num}[1]);
                                 if($eb==1){$Hap{$num}='F1';map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$Hap_support_num{$_}{$num}{'F1'}++}}keys %region}else{$Hap{$num}='-'}}(@Embryo);
                	}elsif($GT{'Proband'}[0] eq $tribase && $GT{'Proband'}[1] eq $tribase){
                	        $Hap{'F2'} = $uniqbase,$Hap{'F1'} = $tribase;
                                if($isPaternal==1 && $GT{'Paternal'}[0] eq $GT{'Paternal'}[1] && $GT{'Paternal'}[0] ne $GT{'Mother'}[0]){
                                        if   ($uniqbase eq $GT{'Paternal'}[0]){$Linkage_Provide='F2'}
                                        elsif($uniqbase ne $GT{'Paternal'}[0]){$Linkage_Provide='F1'}
                                }
				map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$informative{$_}{'F2'}++}
                                }keys %region;
                                map{my$num=$_;
                                my $eb=0;
                                $eb=&EB_correct($uniqbase,$GT{$num}[0],$GT{$num}[1],$DP{$num}[0],$DP{$num}[1]);
                                if($eb==1){$Hap{$num}='F2';
                                map{if($line[0]eq$region{$_}[0]&& $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$Hap_support_num{$_}{$num}{'F2'}++}
                                }keys %region}else{$Hap{$num}='-'}}(@Embryo);
                	}
		}
        	if(defined $Hap{'F1'}){
			my$temp_gene="none";
			map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$temp_gene=$_}}keys %region;
			print HAP "$temp_gene\t$line[0]\t$line[1]\t";
			print HAP "$Hap{'F1'}\t$Hap{'F2'}\t$Hap{'M1'}\t$Hap{'M2'}\t";
			if($isPaternal==1||$isMaternal==1){print HAP "$Linkage_Provide\t"};
			map{print HAP "$Hap{$_}\t"}@Embryo;
			map{print HAP "$GT{$_}[0]\/$GT{$_}[1]\t";}(@sample);
			map{print HAP $DP{$_}[0]+$DP{$_}[1].":$DP{$_}[0]:$DP{$_}[1]\t";}(@sample);
			print HAP "\n";
		}
	}
	
        elsif($core eq 'N' && $isPaternal == 1 && $GT{'Father'}[0] ne $GT{'Father'}[1] && $GT{'Paternal'}[0] eq $GT{'Paternal'}[1] && $GT{'Mother'}[0] eq $GT{'Mother'}[1]){
		my%base=(),my%Hap=();
		map{$base{$_}++}(@{$GT{'Father'}},@{$GT{'Paternal'}});
                my$tribase=""; map {$tribase = $_ if ($base{$_} == 3)} keys %base;
                my$uniqbase="";$uniqbase=$GT{'Father'}[0] eq $tribase ? $GT{'Father'}[1] : $GT{'Father'}[0];
        	$Hap{'F1'}=$uniqbase,$Hap{'F2'}=$tribase;
		if($uniqbase ne $GT{'Mother'}[0]){
			map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$informative{$_}{'F1'}++}
                        }keys %region;
                        map{my$num=$_;
                        my $eb=0;
                        $eb=&EB_correct($uniqbase,$GT{$num}[0],$GT{$num}[1],$DP{$num}[0],$DP{$num}[1]);
                        if($eb==1){$Hap{$num}='F1';
                        map{if($line[0]eq$region{$_}[0]&& $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$Hap_support_num{$_}{$num}{'F1'}++}
                        }keys %region}else{$Hap{$num}='-'}}(@Embryo);
                }elsif($uniqbase eq $GT{'Mother'}[0]){
			map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$informative{$_}{'F2'}++}
                        }keys %region;
                        map{my$num=$_;
                        if(($GT{$num}[0] eq $tribase && $DP{$num}[0]>2) || ($GT{$num}[1] eq $tribase&& $DP{$num}[1]>2 )){$Hap{$num}='F2';
                        map{if($line[0]eq$region{$_}[0]&& $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$Hap_support_num{$_}{$num}{'F2'}++}
                        }keys %region}else{$Hap{$num}='-'}}(@Embryo);
                }
                if(defined $Hap{'F1'}){
                        my$temp_gene="none";
                        map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$temp_gene=$_}}keys %region;
                        print HAP "$temp_gene\t$line[0]\t$line[1]\t";
                        print HAP "$Hap{'F1'}\t$Hap{'F2'}\t";
                        map{print HAP "$Hap{$_}\t"}@Embryo;
                        map{print HAP "$GT{$_}[0]\/$GT{$_}[1]\t";}(@sample);
                        map{print HAP $DP{$_}[0]+$DP{$_}[1].":$DP{$_}[0]:$DP{$_}[1]\t";}(@sample);
                        print HAP "\n";
                }

	}

        elsif($core eq 'N' && $isMaternal == 1 && $GT{'Father'}[0] eq $GT{'Father'}[1] && $GT{'Maternal'}[0] eq $GT{'Maternal'}[1] && $GT{'Mother'}[0] ne $GT{'Mother'}[1]){
                my%base=(),my%Hap=();
                map{$base{$_}++}(@{$GT{'Mother'}},@{$GT{'Maternal'}});
                my$tribase=""; map {$tribase = $_ if ($base{$_} == 3)} keys %base;
                my$uniqbase="";$uniqbase=$GT{'Mother'}[0] eq $tribase ? $GT{'Mother'}[1] : $GT{'Mother'}[0];
                $Hap{'M1'}=$uniqbase,$Hap{'M2'}=$tribase;
                if($uniqbase ne $GT{'Father'}[0]){
                	map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$informative{$_}{'M1'}++}
                        }keys %region;
                        map{my$num=$_;
                        my $eb=0;
                        $eb=&EB_correct($uniqbase,$GT{$num}[0],$GT{$num}[1],$DP{$num}[0],$DP{$num}[1]);
                        if($eb==1){$Hap{$num}='M1';
                        map{if($line[0]eq$region{$_}[0]&& $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$Hap_support_num{$_}{$num}{'M1'}++}
                        }keys %region}else{$Hap{$num}='-'}}(@Embryo);
		}elsif($uniqbase eq $GT{'Father'}[0]){
                	map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$informative{$_}{'M2'}++}
                        }keys %region;
                        map{my$num=$_;
                        if(($GT{$num}[0] eq $tribase && $DP{$num}[0]>2) || ($GT{$num}[1] eq $tribase&& $DP{$num}[1]>2 )){$Hap{$num}='M2';
                        map{if($line[0]eq$region{$_}[0]&& $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$Hap_support_num{$_}{$num}{'M2'}++}
                        }keys %region}else{$Hap{$num}='-'}}(@Embryo);
		}        
        	if(defined $Hap{'M1'}){
                        my$temp_gene="none";
                        map{if($line[0]eq$region{$_}[0] && $line[1]>=$region{$_}[1] && $line[1]<=$region{$_}[2]){$temp_gene=$_}}keys %region;
                        print HAP "$temp_gene\t$line[0]\t$line[1]\t";
                        print HAP "$Hap{'M1'}\t$Hap{'M2'}\t";
                        map{print HAP "$Hap{$_}\t"}@Embryo;
                        map{print HAP "$GT{$_}[0]\/$GT{$_}[1]\t";}(@sample);
                        map{print HAP $DP{$_}[0]+$DP{$_}[1].":$DP{$_}[0]:$DP{$_}[1]\t";}(@sample);
                        print HAP "\n";
                }
	}
}
close VCF;

foreach my$gene(keys %region){
	print SCR "IF\t$gene\t";
	map{if(not defined $informative{$gene}{$_}){$informative{$gene}{$_}=0};print SCR "$informative{$gene}{$_}\t"}('F1','F2','M1','M2');
	print SCR "\n";
	foreach my$em(@Embryo){
		print SCR "$em\t$gene\t";
		map{if(not defined $Hap_support_num{$gene}{$em}{$_}){$Hap_support_num{$gene}{$em}{$_}=0};print SCR "$Hap_support_num{$gene}{$em}{$_}\t"}('F1','F2','M1','M2');
		my$FU=&HapIdentify($informative{$gene}{'F1'},$informative{$gene}{'F2'},$Hap_support_num{$gene}{$em}{'F1'},$Hap_support_num{$gene}{$em}{'F2'});
	#	my$FI=&HapIdentify($informative{$gene}{'F1-I'},$informative{$gene}{'F2-I'},$Hap_support_num{$gene}{$em}{'F1-I'},$Hap_support_num{$gene}{$em}{'F2-I'});
	#	my$FD=&HapIdentify($informative{$gene}{'F1-D'},$informative{$gene}{'F2-D'},$Hap_support_num{$gene}{$em}{'F1-D'},$Hap_support_num{$gene}{$em}{'F2-D'});
		my$MU=&HapIdentify($informative{$gene}{'M1'},$informative{$gene}{'M2'},$Hap_support_num{$gene}{$em}{'M1'},$Hap_support_num{$gene}{$em}{'M2'});
	#	my$MI=&HapIdentify($informative{$gene}{'M1-I'},$informative{$gene}{'M2-I'},$Hap_support_num{$gene}{$em}{'M1-I'},$Hap_support_num{$gene}{$em}{'M2-I'});
	#	my$MD=&HapIdentify($informative{$gene}{'M1-D'},$informative{$gene}{'M2-D'},$Hap_support_num{$gene}{$em}{'M1-D'},$Hap_support_num{$gene}{$em}{'M2-D'});
		print SCR "$FU\|$MU\n";
	}
}

sub HapIdentify{
	my($IF1,$IF2,$HP1,$HP2)=@_;
	my$rlt;
	if($IF1>=5&&$IF2>=5){
		$rlt = $HP1/$IF1>=0.2 && $HP2/$IF2<=0.05 ? '1' : $HP1/$IF1<=0.05 && $HP2/$IF2>=0.2 ? '2' : '-';
	}else{
		$rlt = '-';
	}
	return $rlt;
}

sub GT_correct{
	my($ref,$mut,$base1,$base2,$depth1,$depth2)=@_;
	#my$sum=$depth1+$depth2,
        my$GT1,my$GT2;
	$GT1=($base1==0?$ref:$base1==1?$mut:'-');
	$GT2=($base2==0?$ref:$base2==1?$mut:'-');
#	if($sign==1){
#		if($sum==0){
#			$GT1='-',$GT2='-';
#		}elsif($depth1/$sum>=$pst){
#			$GT1=$ref,$GT2=$ref;
#		}elsif($depth2/$sum>=$pst){
#			$GT1=$mut,$GT2=$mut;
#		}else{
#			$GT1=$ref,$GT2=$mut;
#		}
#		if($depth1<3){$GT1=$mut,$GT2=$mut}
 #               if($depth2<3){$GT1=$ref,$GT2=$ref}
#	}
	return $GT1,$GT2;
}



sub EB_correct{
       my $be=0;
       my($uniq,$GT1,$GT2,$depth1,$depth2)=@_;
       $be=1 if ($GT1 eq $uniq && $GT2 eq $uniq  && $depth1+$depth2>=4);
       if ($GT1 ne $GT2){
       $be=1 if ($GT1 eq $uniq &&  $depth1>$depth2);
       $be=1 if ($GT2 eq $uniq &&  $depth2>$depth1);
       if ($depth1 !=0){
       $be=1 if (($GT1 eq $uniq && $depth1<=$depth2  && $depth1==2||3 && $depth2/$depth1<10) || ($GT2 eq $uniq   && $depth1>=$depth2 && $depth2==2||3 && $depth1/$depth2<10));
       $be=0 if (($GT1 eq $uniq && $depth1<=$depth2  &&  $depth1==1) || ($GT2 eq $uniq && $depth1>=$depth2 && $depth2==1));
#       $eb=1 if (($GT1 eq $uniq && $depth1<$depth2  &&  $depth1==3 && $depth2/$depth1<10) || ($GT2 eq $uniq && $depth1>$depth2 && $depth2==3 && $depth1/$depth2<10));
       $be=1 if (($GT1 eq $uniq && $depth1<=$depth2  &&  $depth1>3 && $depth2/$depth1<40) || ($GT2 eq $uniq && $depth1>=$depth2 && $depth2>3 && $depth1/$depth2<40));
}
}
      return $be;
}

