#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;
use List::Util qw[min max];

warn "\n$0\n";
warn "\n","*********************","\n";


my $numArgs = $#ARGV + 1;

unless( $numArgs == 1) {
    
    print "\n";    
    print "===============================================================================\n";
    print "Usage: perl $0 \n";
    print "\n";
    print "===============================================================================\n";
    print "   Author: Frederic Masclaux\n";
    print "   Date: 2016-07-13\n";
    print "   Version: 1\n";
    print "   Name: $0\n";
    print "   Purpose: \n";
    print "   Convert .vcf file from GATK to .vcf file nearly as formatted by Freebayes\n";
    print "   The input file should be a quality filtered file\n";      
    print "   Comments :\n";
    print "   \n";    
    print "===============================================================================\n";
    exit;
}

my ($VCF_file) = @ARGV;
open (IN1, $VCF_file) || die "Can't open $VCF_file: $!";

#output
my @interm = split(/\./, $VCF_file);
my $NameOUT1 = $interm[0] . ".NearlyFreebayesVCF";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";

print OUT1 "## Nearly-Freebayes-VCF, corresponding to VCFv4.1\n";
print OUT1 "##\n";
print OUT1 "## Field 'INFO' is not changed compared to the field INFO' in GATK VCF file. Only the new sub-field 'TYPE=' is added at the beginning of this field\n";
print OUT1 "## subfields for quality, QR and QA, are not representative of the real quality. A value of 1000 is used as default\n";
print OUT1 "## This file should be generated from a GATK *filtered VCF* (filtered for a minimum coverage and a miminum calling quality)\n";

#.....Variables.....#

my $CHROM;
my $POS;
my $ID;
my $REF;
my $ALT;
my $QUAL;
my $FILTER;
my $INFO;
my $GENOFORMAT;
my $GENO;

my $gatkGT;
my $gatkAD;
my $gatkDP;
my $gatkGQ;
my $gatkPL;

my $freebGT;
my $freebDP;
my $freebRO;
my $freebQR;
my $freebAO;
my $freebQA;

my $freebALT;
my $freebTYPE;

my @interm2;
my @interm3;

my $cursor;
my $test;

my $temp_ratio;

my $counter;

print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                  Parsing the file $VCF_file\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

foreach my $line (<IN1>) {
    chomp($line);
   
    if 	($line =~ /^\s*#/) {
        next;
       
    } elsif ($line =~ /^\s*$/) {
        next;
       
    } else {
        @interm = split(/\t/,$line);
             
        $CHROM = $interm[0];
        $POS = $interm[1];
        $ID = $interm[2];
        $REF = $interm[3];
        $ALT = $interm[4];
        $QUAL = $interm[5];
        $FILTER = $interm[6];
        $INFO = $interm[7];
        
        $GENOFORMAT = $interm[8];
        $GENO = $interm[9];

        @interm = split(/\,/,$ALT);
        
        @interm2 = split(/\:/,$GENO);
        if (@interm2 == 5 and $interm2[2] > 0) {    # GATK variant files have problems : 1/ sometimes not all the fields in $GENO are reported, 2/ Total DP can be 0)
            $counter ++;
            $gatkGT = $interm2[0];
            $gatkAD = $interm2[1];
            $gatkGQ = $interm2[3];
            $gatkPL = $interm2[4];
            
            @interm3 = split(/\,/,$gatkAD);        
                   
            #define genotype (like 0/1/1/1)
            $freebGT = "";
            $gatkDP = 0;
            foreach my $occurrence (@interm3) {
                $gatkDP += $occurrence;
            }

            $cursor = 0;
            foreach my $occurrence (@interm3) {
                $temp_ratio = sprintf("%.2f", int($occurrence/$gatkDP*10 + 0.5));
                if ($temp_ratio > 0) {
                    $freebGT .= ($cursor . "/") x $temp_ratio;
                }
                $cursor ++;
            }
            chop($freebGT);
            
                
            # reference allele
            $freebRO = shift @interm3;
            
            # alleles
            $cursor = 0;
            $freebALT = "";
            $freebAO = "";
            $freebTYPE  = "";
            
            foreach my $allele (@interm) {
                $freebALT .= $allele . ",";
                $freebAO .= $interm3[$cursor] . ",";
                if (length($allele) == 1 and length($REF) == 1) {
                    $freebTYPE .= "snp,";
                } elsif (length($REF) > length($allele)) {
                    $freebTYPE .= "del,";
                } elsif (length($REF) eq length($allele)) {
                    $test = compare_alleles($REF,$allele);
                    if ($test eq "mnp" ) {
                        $freebTYPE .= "mnp,";
                    } else {
                        $freebTYPE .= "snp,";
                    }
                    
                } elsif (length($REF) < length($allele)) {
                    $freebTYPE .= "ins,";
                }
                $cursor ++;
            }
            chop($freebALT);
            chop($freebAO);
            chop($freebTYPE);
           
            print OUT1 $CHROM, "\t", $POS, "\t", $ID, "\t", $REF, "\t", $freebALT, "\t", $QUAL, "\t", $FILTER, "\t",
            "TYPE=", $freebTYPE, ";", $INFO,  "\t", "GT:DP:RO:QR:AO:QA", "\t", $freebGT, ":", $gatkDP, ":", $freebRO, ":1000:", $freebAO, ":1000", "\n";
            
        } else {
            print "This position is problematic:\n", $line, "\n", "=> Not processed\n", "-" x 30, "\n";
        }
    }
}

print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                            Complete \n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

print "$counter positions were processed\n\n";

    
sub compare_alleles {
    my ($ref_sub, $allele_sub) = @_;
    
    my @interm_REF_sub = split(//,$ref_sub);
    my @interm_allele_sub = split(//,$allele_sub);
    my $length_sub = @interm_REF_sub;
    my $cursor_sub = 0;
    
    while ($cursor_sub < $length_sub) {
        if ($cursor_sub != 0 and $interm_REF_sub[$cursor_sub] ne $interm_allele_sub[$cursor_sub]) {
            $test = "mnp";
        }
        $cursor_sub ++;
    }  
}
 
 
 __END__
 
 
 