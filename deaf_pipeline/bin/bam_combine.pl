#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);

die "perl $0 <hgbam> <segbam> <bed> > [output]"unless @ARGV==3;

my (%seg_c, %seg_i, %seg_l);
open IN,$ARGV[2];
while(<IN>){
    chomp;
    my @a = split /\t/;
    $seg_i{$a[3]} = "$a[0] $a[1] $a[2]";
    $seg_l{$a[3]} = $a[2] - $a[1];
    $seg_c{$a[0]}{$a[3]} = "$a[1]-$a[2]";
}
close IN;

my $samtools = "$Bin/samtools";

my %header_tag;
open HEAD, "$samtools view -H $ARGV[0] | ";
while(<HEAD>){
    chomp;
    my @tmp = split;
    my $type = $tmp[0];
    if ($type eq '@SQ') {
        my $chrn = (split /chr/, $tmp[1])[1];
        next unless $chrn eq 1 || $chrn eq 7 || $chrn eq 13 || $chrn eq 'M';
    }elsif (/^\@PG\tID:tmap/){
    }else {
        next if $header_tag{$type};
        $header_tag{$type} = 1;
    }
    print "$_\n";
}
close HEAD;

open HGB, "$samtools view $ARGV[0] | ";
open SGB, "$samtools view $ARGV[1] | ";
my $hg_line = <HGB>;
my $sg_line = <SGB>;
my $last_name = (split /\t/, $hg_line)[0];
while(<HGB>){
    my $name = (split)[0];
    if($name eq $last_name){
        $hg_line .= $_;
        while(<HGB>){
            $name = (split)[0];
            if($name eq $last_name){
                $hg_line .= $_;
            }
            else{
                last;
            }
        }
    }
    my ($sgl, $sg_name);
    while($sgl = <SGB>){
        $sg_name = (split /\t/, $sgl)[0];
        if($sg_name eq $last_name){
            $sg_line .= $sgl;
        }
        else{
            last;
        }
    }
    die "$name $sg_name" if $sg_name && $name ne $sg_name;
    combine_alignment( $hg_line, $sg_line );
    ($last_name, $hg_line, $sg_line) = ($name, $_, $sgl);
}
close HGB;
close SGB;

sub combine_alignment {
    my ($hgl, $sgls) = @_;
    chomp( $hgl );
    my ($flag, $out_line) = (0, $hgl);
    my @hfld = split /\t/, $out_line;
    my ($chr, $pos, $cigar1) = @hfld[2,3,5];
    my $seg_h = my $seg_s = '*';
    for my $s (sort keys %{$seg_c{$chr}})    {
        my ($ss, $es) = split /-/, $seg_c{$chr}{$s};
        my $end = &calc_end( $pos, $cigar1 );
        next if $pos > $es || $end < $ss;
        $seg_h = $s;
        last;
    }

    for my $sgl (split /\n/, $sgls) {
        my @sfld = split /\t/, $sgl;
        next if $sfld[2] ne $seg_h;

        $seg_s = $sfld[2];
        my $mq_hg = $hfld[4];
        my $mq_sg = $sfld[4];

        if ($mq_sg > $mq_hg)    {
            my ($as_sg, $xs_sg);
            for my $f (@sfld)   {
                $as_sg = $f if $f =~ /^AS:i/;
                $xs_sg = $f if $f =~ /^XS:i/;
            }
            for my $f (@hfld)   {
                $f = $as_sg if $f =~ /^AS:i/;
                $f = $xs_sg if $f =~ /^XS:i/;
            }
            $hfld[4] = $mq_sg;
        }
        $out_line = join( "\t", @hfld );
        $flag ++;
    }
    print "$out_line\n";
}

sub calc_end    {
    my ($pos, $cigar) = @_;
    $pos --;
    $cigar =~ s/\d+[INSHP]//g;
    my @pl = split /[MD]/, $cigar;
    for my $p (@pl) {
        $pos += $p;
    }
    return $pos;
}


