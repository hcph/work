#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin $Script);
die "perl $0 <indir> <exp_file> <outdir> <name> <20pos>"unless @ARGV >= 4;

my (%pos_info, @pos_info);
my (%exp_order, %exp_info);
my (%detect_info);
my (%depth_info);
my $samtools  = "$Bin/samtools";
my $pos_file = (@ARGV==5) ? $ARGV[4] : "$Bin/../database/target_20variants.xls";

open IN,$pos_file;
while(<IN>){
    next if /^#/;
    chomp;
    my @a = split;
    push @pos_info, $a[4];
    my $ps = $a[1];
    my $pe = $a[4] =~ /del/ ? $a[1] + length($a[5]) : $a[1];
    $pos_info{$a[4]} = "$a[0]:$ps-$pe";
}
close IN;

open IN,$ARGV[1];
# Barcode_25  ID_001  16D0985274  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT  WT
while(<IN>){
    chomp;
    my @a = split /\t/;
    if( /^barcode\s+index/i ){
        for(my $i = 3; $i < @a; $i++){
            $exp_order{$i} = $a[$i];
        }
        next;
    }
    my $key1 = (split /\_/, $a[0])[1]; $key1 =~ s/^0+//;
    my $key2 = (split /\_/, $a[1])[1]; $key2 =~ s/^0+//;
    $key1 = ($key1 < 10) ? "0$key1" : $key1;
    #$key2 = ($key2 < 10) ? "00$key2" : "0$key2";
    $key2 = "$key2" if (length($key2)==3);
    $key2 = "0$key2" if (length($key2)==2);
    $key2 = "00$key2" if (length($key2)==1);
    $exp_info{$key1}{$key2} = $_;
}
close IN;

open IN,"$ARGV[2]/$ARGV[3].detect_sample_stat.xls";
# barcode_25      ID_0001 16D0985274      0       20      0       0       0        100
while(<IN>){
    chomp;
    next if /^Barcode/;
    my @a = split /\t/;
    my $key1 = (split /\_/, $a[0])[1]; $key1 =~ s/^0+//;
    my $key2 = (split /\_/, $a[1])[1]; $key2 =~ s/^0+//;
    $key1 = ($key1 < 10) ? "0$key1" : $key1;
    #$key2 = ($key2 < 10) ? "00$key2" : "0$key2";
    $key2 = "$key2" if (length($key2)==3);
    $key2 = "0$key2" if (length($key2)==2);
    $key2 = "00$key2" if (length($key2)==1);
    $detect_info{$key1}{$key2} = [@a];
}
close IN;

open OT,">$ARGV[2]/$ARGV[3].statistics.xls";
print OT "Barcode\tIndex\t#Exp_site\t#Detect_site\t#Missing_site\t%Detect_rate\t#TP\t#TN\t#FP\t#FN\t%Correct_rate";
for(my $i = 0; $i < @pos_info; $i++){
    print OT "\tDepth_$pos_info[$i]";
}
print OT "\n";
foreach my $barcode( `ls $ARGV[0] | grep barcode | grep -v xls` ){
    chomp $barcode;
    my $key1 = (split /\_/, $barcode)[1];
    $key1 =~ s/^0+//; $key1 = ($key1 < 10) ? "0$key1" : $key1;
    for(my $i = 1; $i <= 128; $i++){
        my @out = ();
        #my $key2 = ($i < 10) ? "00$i" : "0$i";
        my $key2;
	$key2 = "$i" if (length($i)==3);
        $key2 = "0$i" if (length($i)==2);
        $key2 = "00$i" if (length($i)==1);
        if(exists $exp_info{$key1}{$key2}){
            $out[0] = exp_count( $exp_info{$key1}{$key2} );
        }
        else{
            $out[0] = 0;
        }
        if( $detect_info{$key1}{$key2}[2] ne "-" ){
            $out[1] = $detect_info{$key1}{$key2}[3] + $detect_info{$key1}{$key2}[4] + $detect_info{$key1}{$key2}[5] + $detect_info{$key1}{$key2}[6];
        }
        else{
            $out[1] = "-";
        }
        $out[2] = $detect_info{$key1}{$key2}[7];
        $out[3] = ($out[0] == 0) ? 0 : int(($out[1]-$out[2]) / $out[0] * 10000 + 0.5) / 100;
        $out[4] = $detect_info{$key1}{$key2}[3];
        $out[5] = $detect_info{$key1}{$key2}[4];
        $out[6] = $detect_info{$key1}{$key2}[5];
        $out[7] = $detect_info{$key1}{$key2}[6];
        $out[8] = $detect_info{$key1}{$key2}[8];
        my $bam = "$ARGV[0]/barcode_$key1/index_$key2.5.final.bam";
        my $vcf = "$ARGV[0]/barcode_$key1/index_$key2.5.final.bam.vcf";
        next unless -f $vcf;
        my $depth = 0;
        for(my $j = 0; $j < @pos_info; $j++){
            open IN,$vcf;
            while(<IN>){
                chomp;
                next if /^#/;
                my @a = split /\t/;
                next unless $a[2] eq $pos_info[$j];
                if( $a[-1] =~ /\.\/\./){
                    $depth = 0;
                }else{
                    my @b = split /\:/, $a[-1];
                    my @c = split /\,/, $b[1];
                    $depth = $c[0] + $c[1];
                }
            }
            close IN;
=head
            system("$samtools depth $bam -r $pos_info{$pos_info[$j]} > $ARGV[2]/$ARGV[3].tmp");
            my @depth = (0) x 2;
            open TMP,"$ARGV[2]/$ARGV[3].tmp";
            while(<TMP>){
                chomp;
                next unless /\S+/;
                my @a = split;
                $depth[0] += $a[2];
                $depth[1] += 1;
            }
            close TMP;
=cut
            
            $out[$j+9] = $depth;
        }
        my $out = join "\t", @out;
        print OT "$key1\t$key2\t$out\n";
    }
}
close OT;

sub exp_count{
    my ($exp) = (@_);
    my $num = 0;
    my @array = (split /\t/, $exp);
    for(my $i = 3; $i < @array; $i++){
        $num++ unless $array[$i] eq "NA";
    }
    return $num;
}


