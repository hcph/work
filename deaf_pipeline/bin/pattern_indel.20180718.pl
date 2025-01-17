#!/usr/bin/perl -w
use strict;
use FindBin '$Bin';

die "perl $0 bam (hotspot_indel_file wing_len)\n" if @ARGV < 1;

my ($bam, $indel_info, $wing_len) = @ARGV;
$indel_info ||= "$Bin/../database/del_sites.info";
$wing_len ||= 50;
my $min_half_len = 10;
my $samtools = -f "$Bin/samtools" ? "$Bin/samtools" : "samtools";

open INFO, $indel_info;
while (<INFO>)  {
    chomp;
    my ($chr, $target, $name, $query) = (split)[0..3];
    my $target_length = length($query) - 1;
    open SAM, "$samtools view $bam $chr:$target-$target | ";
    my (%check, $total_num, $effect_num);
    $total_num = 0;     ###########
    while (<SAM>)   {
        my @info = split;
        my ($pos, $cigar) = @info[3,5];
        my $copy = $cigar;
        $copy =~ s/\d+S//g;
        if ($copy =~ /[^\dMID]/)    {
            print STDERR "CIGAR error:\t$info[0]\t$chr\t$pos\t$cigar\n";
            next;
        }
        my ($p_pos, $pattern) = &substr_cigar( $copy, $pos, $target, $wing_len );
        next if $p_pos == -1;
        $total_num ++;
        $check{"$p_pos\t$pattern"} ++;
    }
    close SAM;

    my ($total ,$num_d) = (0, 0);
    for my $k (keys %check) {
        next if $check{$k} < 2;
        my ($pos, $pattern) = split /\t/, $k;
        my @p_len = split /[MID]/, $pattern;
        next if @p_len > 3;
        $total += $check{$k};
        $num_d += $check{$k} if &judge_cigar( $pattern, $pos, $target, $target_length ) eq 'DEL';
    }
    my $ratio = $total > 0 ? sprintf( "%.2f", $num_d/$total ) : 'NA';
    print "$name\t$num_d\t$total\t$ratio\n";
}
close INFO;

sub substr_cigar    {
    my ($cigar, $pos, $target, $wing_len) = @_;
    my $left_half_len = $target-$pos+1;
    return -1 if $left_half_len < $min_half_len;
    my $pat_len = $wing_len*2;
    my $start = $target-$wing_len;
    my $end = $target+$wing_len-1;
    my @c_len = split /[MID]/, $cigar;
    my @c_char = split /\d+/, $cigar;
    shift @c_char;
    my ($temp, $pattern) = $pos-1;
    for my $i (0..$#c_len)  {
        $temp += $c_len[$i] if $c_char[$i] ne 'I';
        next if $temp < $start;
        if ($temp > $end)   {
            my $len_c = $c_len[$i];
            my $len_t = $temp-$end;
            if ($len_t < $len_c)    {
                my $len = $len_c-$len_t;
                die "$len = $len_c-$len_t" if $len < 0;
                $len = $pat_len if $len > $pat_len;
                $pattern .= $len.$c_char[$i];
            }
            last;
        }else   {
            my $len_c = $c_len[$i];
            my $len_t = $temp-$start+1;
            my $len = $len_t < $len_c ? $len_t : $len_c;
            die "$len = $len_t < $len_c ?" if $len < 0;
            $pattern .= $len.$c_char[$i];
        }
    }
    my $right_half_len = $temp-$target+1;
    return -1 if $right_half_len < $min_half_len;
    my $p_pos = $pos < $start ? $start : $pos;
    return( $p_pos, $pattern );
}

sub judge_cigar {
    my ($cigar, $pos, $target, $target_length) = @_;
    my @c_len = split /[MID]/, $cigar;
    my @c_char = split /\d+/, $cigar;
    $c_char[0] = 'M';
    unshift @c_len, 0;
    my ($temp, $judge) = $pos-1;
    for my $i (0..$#c_len-1)    {
        $temp += $c_len[$i] if $c_char[$i] ne 'I';
        my ($j, $temp2) = ($i+1, $temp);
        $temp2 += $c_len[$j];
        if ($temp == $target && $c_len[$j] == $target_length)   {
            $judge = $c_char[$j] eq 'M' ? "REF" :
                     $c_char[$j] eq 'D' ? "DEL" :
                     $c_char[$j] eq 'I' ? "INS" : "ERR";
            last;
        }
    }
    $judge ||= $pos < $target && $temp > $target ? "REF" : "ERR";
    return $judge;
}
