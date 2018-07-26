#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin $Script);
die "perl $0 <indir> <pos_file> <exp_file> <outfile>"unless @ARGV==4;
#eg: 

my (%pos_info, %wt_info, @header_list);
my ($barcode_id, $index_id, %barcode_id);
my (%info_vcf, @out);
my (%exp_info, %exp_order, %th, %tpr, $tp, $fp, $fn, $tpr, $tn, $sum, $miss, %error_stat, $error_sum);
my $samtools = "$Bin/samtools";
my $ref_file = "$Bin/../database/hg19.fa";
#my $bin_dir = "/zfssz3/MGI_BIT/RUO/dushiyi/script/Program/05.forQseq/HL/bin";
my $threshold = "$Bin/../database/threshold.txt";
my $cutoff = 30;
my $del_cutoff = 20;
my @flag_id = qw/Type GT AD GQ/;
my ($cutoff2, %roc);

open THE,$threshold;
while(<THE>){
    chomp;
    my @a = split;
    $th{$a[0]} = $_;
}
close THE;

open POS,$ARGV[1];
while(<POS>){
    chomp;
    next if /^#/;
    my @field = split /\s+/;
    $pos_info{$field[0]}{$field[1]} = [@field];
    $wt_info{$field[4]} = [@field];
    push @header_list, $field[4];
}
close POS;

open EXP,"$ARGV[2]";
while(<EXP>){
    chomp;
    my @field = split /\t/;
    if(/^barcode\s+index/i){
        for( my $k = 3; $k < @field; $k++ ){
            $exp_order{$k} = $field[$k];
        }
        next;
    }
    my $key1 = (split /\_/, $field[0])[1]; $key1 =~ s/^0+//;
    my $key2 = (split /\_/, $field[1])[1]; $key2 =~ s/^0+//;
    #my $key2 = (split /\_/, $field[1])[1]; $key2 =~ s/^\d//;
    $key1 = ($key1 < 10) ? "0$key1" : $key1;
    #$key2 = ($key2 < 10) ? "00$key2" : $key2;
    $key2 = "$key2" if (length($key2)==3);
    $key2 = "0$key2" if (length($key2)==2);
    $key2 = "00$key2" if (length($key2)==1);
    for( my $n = 3; $n < @field; $n++ ){
        if( $field[$n] !~ /WT|Hetero|Homo|unknown|NA/ ){
            $field[$n] = (length($field[$n]) == 1) ? "$field[$n]$field[$n]" : $field[$n];
            if( $field[$n] eq $wt_info{ $exp_order{$n} }[7] ){ $field[$n] = $wt_info{ $exp_order{$n} }[8]; }
            elsif( $field[$n] eq $wt_info{ $exp_order{$n} }[9] ){ $field[$n] = $wt_info{ $exp_order{$n} }[10]; }
            elsif( $field[$n] eq $wt_info{ $exp_order{$n} }[11] ){ $field[$n] = $wt_info{ $exp_order{$n} }[12]; }
        }
    }
    $exp_info{$key1}{$key2} = join "\t", @field;
}
close EXP;

open OT1,">$ARGV[3].file_for_next.xls";
print OT1 "Barcode\tIndex\tSample\tInfo\t".(join "\t",@header_list)."\n";

foreach my $barcode( `ls $ARGV[0] | grep barcode | grep -v xls` ){
    chomp $barcode;
    $barcode_id = (split /\_/, $barcode)[1];chomp $barcode_id;
    $barcode_id =~ s/^0+//; $barcode_id = ($barcode_id < 10) ? "0$barcode_id" : $barcode_id;
    $barcode_id{$barcode_id} = 1;

    for( my $i = 1; $i <= 128; $i++ ){
        @out = ();
        $index_id = ($i < 10) ? "00$i" : "0$i";
        $index_id = "$i" if (length($i)==3);
        $index_id = "0$i" if (length($i)==2);
        $index_id = "00$i" if (length($i)==1);
        next unless -e "$ARGV[0]/$barcode/index_$index_id.5.final.bam.vcf";
        open VCF,"$ARGV[0]/$barcode/index_$index_id.5.final.bam.vcf";
        my $bam = "$ARGV[0]/$barcode/index_$index_id.5.final.bam";
        my $tpr_file = "$ARGV[0]/$barcode/index_$index_id.5.final.bam.tpr";
        open TF, $tpr_file;
        while (<TF>)    {
            chomp;
            my ($site, $nid, $nt, $rate) = split;
            $tpr{$site} = "$nid;$nt;$rate";
        }
        close TF;
        while(<VCF>){
            chomp;
            next if /^#/;
            my @field = split /\s+/;
            next unless exists $pos_info{$field[0]}{$field[1]};
            my ($chr, $pos, $ref, $alt) = @field[0,1,3,4];
            my $raw_dp = 0;
            if ($field[7] =~ /DP=(\d+)/)    { $raw_dp = $1; }
            my @result_fld = split /:/, $field[9];
            my ($gt_word, $gq_word) = @result_fld[0,3];
            my $ad_word = $result_fld[1] ? $result_fld[1] : '0,0';
            my ($ad_ref, $ad_alt) = split /,/, $ad_word;
            my $sum_dp = $gt_word eq './.' ? 0 : $ad_ref+$ad_alt;
            my ($gt, $ad, $gq);
            my $site = $field[2];
            $cutoff2 = ($sum_dp == 0) ? 0 : $ad_alt / $sum_dp;
            if ($site !~ /del/) {
                if ($sum_dp > $cutoff)  {
                    $gt = $gt_word eq '0/0' ? 'WT' :
                        $gt_word eq '0/1' ? 'Hetero' :
                        $gt_word eq '1/1' ? 'Homo' : 'Err';
                    ($ad, $gq) = ("$ad_ref;$ad_alt", $gq_word);
                }
                elsif($raw_dp > $cutoff)  { # use mpileup
                    print STDERR "cjeck: $chr:$pos-$pos\n";
                    my $mpileup = `$samtools mpileup -d 20000 -f $ref_file -r $chr:$pos-$pos $bam`;
                    if ($mpileup eq '') {
                        ($gt, $ad, $gq) = ('NA', '0;0', '0');
                    }else   {
                        my $pileup = (split /\t/, $mpileup)[4];
                        my $copy = $pileup;
                        $copy =~ s/\^.//g;
                        $copy =~ s/\$//g;
                        while ($copy =~ /([+-])(\d+)/)  {
                            my ($sign, $num) = ($1, $2);
                            if ($sign eq '+')   {
                                $copy =~ s/.\+$num([ACGTNacgtn]{$num})/+/;
                            }else   {
                                $copy =~ s/.-$num([ACGTNacgtn]{$num})/-/;
                            }
                        }
                        my @char = split //, $copy;
                        my @pos = grep {/\./} @char;
                        my @neg = grep {/,/} @char;
                        my $alt = $site;
                        $alt =~ s/.+(.)$/$1/;
                        my $alt_n = $alt;
                        $alt_n =~ tr/ACGT/acgt/;
                        my @alt_p = grep {/$alt/} @char;
                        my @alt_n = grep {/$alt_n/} @char;
                        my $total_dep = @char;
                        my $ref_dep = @pos+@neg;
                        my $alt_dep = @alt_p+@alt_n;
                        print "$barcode_id\t$index_id\t$chr\t$pos\t$raw_dp\t$ref_dep\t$alt_dep\n";
                        my $alt_r = $alt_dep/($ref_dep+$alt_dep);
                        $gt = &judge( $alt_r, $site );
                        $ad = "$ref_dep;$alt_dep";
                        $gq = 99;
                        $cutoff2 = $alt_r;
                    }
                }
                else{
                    $gt = $ad = "NA";
                    $gq = 0;
                }
            }
            else{ # TPR for del, don`t understand!
                if ($tpr{$site}) {
                    my ($d1, $dt, $rate) = split /;/, $tpr{$site};
                    $cutoff2 = $rate;
                    if ($rate eq 'NA')  {
                        ($gt, $ad, $gq) = ('NA', 'NA', '0');
                    }elsif ($dt < $del_cutoff)  {
                        my $d0 = $dt-$d1;
                        ($gt, $ad, $gq) = ('NA', "$d0;$d1", '0');
                    }else   {
                        my $d0 = $dt-$d1;
                        $gt = &judge( $rate, $site );
                        $ad = "$d0;$d1";
                        $gq = ($dt < $del_cutoff) ? 0 :
                            ($gt eq 'NA')       ? 0 : 99;
                    }
                }else   {
                    ($gt, $ad, $gq) = ("NA", "NA", "0");
                }
            }

            if( $gt eq "NA"){
                $info_vcf{$barcode_id}{$index_id}{ $pos_info{$field[0]}{$field[1]}[4] }[0] = $gt;
                $info_vcf{$barcode_id}{$index_id}{ $pos_info{$field[0]}{$field[1]}[4] }[1] = $gt;
                $info_vcf{$barcode_id}{$index_id}{ $pos_info{$field[0]}{$field[1]}[4] }[2] = $ad;
                $info_vcf{$barcode_id}{$index_id}{ $pos_info{$field[0]}{$field[1]}[4] }[3] = $gq;
            }
            else{
                my $gt_type;
                if( $gt eq "WT" ){
                    $gt_type = 7;
                }
                elsif( $gt eq "Hetero" ){
                    $gt_type = 9;
                }
                elsif( $gt eq "Homo" ){
                    $gt_type = 11;
                }
                $info_vcf{$barcode_id}{$index_id}{ $pos_info{$field[0]}{$field[1]}[4] }[0] = $pos_info{$field[0]}{$field[1]}[ $gt_type + 1 ];
                $info_vcf{$barcode_id}{$index_id}{ $pos_info{$field[0]}{$field[1]}[4] }[1] = $pos_info{$field[0]}{$field[1]}[ $gt_type ];
                $info_vcf{$barcode_id}{$index_id}{ $pos_info{$field[0]}{$field[1]}[4] }[2] = $ad;
                $info_vcf{$barcode_id}{$index_id}{ $pos_info{$field[0]}{$field[1]}[4] }[3] = $gq;
            }
            $roc{$barcode_id}{$index_id}{ $pos_info{$field[0]}{$field[1]}[4] } = $cutoff2;

        }
        close VCF;

        foreach my $pos( @header_list ){
            if( defined $info_vcf{$barcode_id}{$index_id}{ $pos } ){
                for(my $j = 0; $j < 4; $j++){
                    push @{$out[$j]}, $info_vcf{$barcode_id}{$index_id}{$pos}[$j];
                }
            }
            else{
                $info_vcf{$barcode_id}{$index_id}{ $pos }[0] = "WT";
                for(my $j = 0; $j < 4; $j++){
                    if( $j == 1 ){
                        push @{$out[$j]}, $wt_info{$pos}[7];
                    }
                    else{
                        push @{$out[$j]}, "WT";
                    }
                }
            }
        }
        my $sample_name = (defined $exp_info{$barcode_id}{$index_id}) ? (split /\s+/, $exp_info{$barcode_id}{$index_id})[2] : "-";
        for(my $j = 0;  $j < 4; $j++){
            print OT1 "barcode_${barcode_id}\tID_${index_id}\t$sample_name\t$flag_id[$j]\t".(join "\t", @{$out[$j]})."\n";
        }

    }
}

close OT1;

open OT2,">$ARGV[3].detect_sample_stat.xls";
print OT2 "Barcode\tIndex\tSample\tTP\tTN\tFP\tFN\tMISS\tTPR(%)\n";
open ROC,">$ARGV[3].roc_stat.xls";
print ROC "Barcode\tIndex\tSite\tExp\tRate\n";

foreach my $key1 ( sort keys %barcode_id ){
    for(my $key2 = 1; $key2 <= 128; $key2++){
        ($tp, $fp, $fn, $tn, $sum, $miss, $tpr) = (0) x 7;
        $key1 =~ s/^0+//; $key1 = ($key1 < 10) ? "0$key1" : $key1;
        $key2 =~ s/^0+//; #$key2 = ($key2 < 10) ? "00$key2" : "0$key2";
        $key2 = "$key2" if (length($key2)==3);
        $key2 = "0$key2" if (length($key2)==2);
        $key2 = "00$key2" if (length($key2)==1);
        my $error_flag = 0;
        my (@roc_op1, @roc_op2) = () x 2;my $roc_type = "N";
        if( defined $exp_info{$key1}{$key2} ){
            my @field = split /\t/, $exp_info{$key1}{$key2};
            for(my $k = 3; $k < @field; $k++){
                if( !defined $info_vcf{$key1}{$key2}{$exp_order{$k}}[0] ){
                    $error_flag++;
                    $info_vcf{$key1}{$key2}{$exp_order{$k}}[0] = "NULL";
                }
                $miss++ if $info_vcf{$key1}{$key2}{$exp_order{$k}}[0] eq "NA";
                if( $field[$k] eq "unknown" ){
                    if( $info_vcf{$key1}{$key2}{$exp_order{$k}}[0] eq "NA" || $info_vcf{$key1}{$key2}{$exp_order{$k}}[0] eq "WT" ){
                        $tn++;
                    }
                    else{
                        $tp++;
                    }
#next;
                }
                if( $field[$k] eq $info_vcf{$key1}{$key2}{$exp_order{$k}}[0] ){
                    if( $field[$k] eq "WT" || $field[$k] eq "NA" ){ $tn++; }
                    else{ $roc_type = "Y";$tp++; }
                }
                else{
                    if( $field[$k] eq "WT" ){
                        if( $info_vcf{$key1}{$key2}{$exp_order{$k}}[0] eq "NA" ){
                            $fn++;
                        }
                        else{
                            $fp++;
                        }
                    }
                    elsif( $field[$k] eq "NA" ){
                        if( $info_vcf{$key1}{$key2}{$exp_order{$k}}[0] eq "WT" || $info_vcf{$key1}{$key2}{$exp_order{$k}}[0] eq "NULL" ){
                            $fn++;
                        }
                        else{
                            $fp++;
                        }
                    }
                    else{ $roc_type = "Y";$fn++; }
                }
                print ROC "barcode_$key1\tID_$key2\t$exp_order{$k}\t$roc_type\t$roc{$key1}{$key2}{$exp_order{$k}}\n";
            }
            $sum = $tp+$tn+$fp+$fn;
            $tpr = ($sum == 0) ? "-" : int( ($tp+$tn) / $sum * 10000 + 0.5 ) / 100;
            print OT2  "barcode_$key1\tID_$key2\t$field[2]\t$tp\t$tn\t$fp\t$fn\t$miss\t$tpr\n";
            print "ERROR SAMPLE: barcode $key1 index $key2: $ARGV[0]\n" if $error_flag > 0;
        }
        else{
            print OT2 "barcode_$key1\tID_$key2\t-\t-\t-\t-\t-\t-\t-\n";
        }
    }
}
close OT2;
close ROC;

sub judge   {
    my ($ratio, $site) = @_;
    my $thres_line = $th{$site} ? $th{$site} : $th{'default'};
    my @thres = split /\s+/, $thres_line;
    return my $type = $ratio > $thres[4] ? 'Homo' :
        $ratio < $thres[1] ? 'WT' :
        $ratio < $thres[3] && $ratio > $thres[2] ? 'Hetero' : 'NA';
}


