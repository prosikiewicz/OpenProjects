#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;
use List::Util qw[min max];

warn "\n$0\n"; 				# nom du fichier
warn "\n","*********************","\n";


my $numArgs = $#ARGV + 1;

unless( $numArgs == 0) {
    
    print "\n";    
    print "===============================================================================\n";
    print "Usage: perl $0 \n";
    print "\n";
    print "===============================================================================\n";
    print "   Author: Frederic Masclaux\n";
    print "   Date: 2016-06-08\n";
    print "   Version: 1\n";
    print "   Name: $0\n";
    print "   Purpose: \n";
    print "   Summarize all positions from VCF files\n";
    print "   \n";      
    print "   Comments :\n";
    print "   Automatically reads all *Q30.vcf files in the directory :\n";    
    print "===============================================================================\n";
    exit;
}


my $some_dir = ".";
opendir(DIR, $some_dir) || die "can't opendir $some_dir: $!";
my @files = grep { /Q30.vcf/ } readdir(DIR);
closedir DIR;

warn "VCF files\n";
for my $i (@files) {
   warn $i,"\n";
    
}
warn "\n";




#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);
my $NameOUT1 = "All_positions_with_variants" . $date . ".combinedVCF";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";



#.....Variables.....#

my $current_rad ="";

my @int1;
my @int2;

my (@interm, @interm2,@interm3, @interm4);
   
my $scaffold; 
my $position;
my $type;
my $coverage;
my @ALLallelesID;
my $value;

my $cursor_type;

my @test_type;
my $REFallele;
my @ALTalleleList;

my %AlleleDist;
my %record;

my @test_ref;
my $test_cond;

my @ordered_isolates;
my $str_ordered_isolates;

warn "\n","*********************","\n";


#####  Parsing VCF file  ##########################

warn "---Parsing vcf file :\n\n";

foreach my $f ( @files) {
   open (IN2, "<$f");
   warn " treatment of $f\n";
   

   #sample name
   @int1 = split(/\./,$f);
   @int2 = split(/\_/,$int1[0]);
   my $SampleName = $int2[2];

   

    
   foreach my $line (<IN2>) {
      chomp($line);
       
      if 	($line =~ /^\s*#/) {
         next;
           
      } elsif ($line =~ /^\s*$/) {
         next;
           
      } else {
         $line =~ s/:/;/g;		# all : used by freebayes in its specific fields are replaced with ; as in the standard VCF
         $line =~ s/,/:/g;		# all , are replaced by : . They are used to separate 2 ALT alleles. , is a problem for table (@...)
                 
         @interm = split(/\s/,$line);
                 
         #Scaffold and position
         $scaffold = $interm[0];
         $position = $interm[1];
 
                 
         #	# extract TYPE=
         @interm2 = split(/;/,$interm[7]);
         $cursor_type = 0;
         foreach my $t (@interm2){
            if ($t =~ m/TYPE/) {
               last;
            }
         $cursor_type ++
         }
         @interm3 = split(/=/,$interm2[$cursor_type]);
         if ($interm3[0] eq "TYPE"){
             $type = $interm3[1];
         } else {
            warn "error in recording TYPE"; die;
         }
             
         @test_type = split(/:/ , $type);
   
          
           
       #	# extract allele value-ID
         @interm2 = split(/;/,$interm[9]);
         @ALLallelesID = split (/\//,$interm2[0]);

         my @unique = ();
         my %seen =();
 
         foreach my $t (@ALLallelesID) {		# Report unique values of 'type'
            if (! $seen{$t}++ ) {
               push @unique, $t;
            }
         }
         @ALLallelesID = @unique;
 
         
         
     #	# extract alleles
         $REFallele = $interm[3];
         @ALTalleleList = split(/:/,$interm[4]);
 
   ############# Control if REF allele is not present in ALT alleles
         # # # #
         my $validation = 1;
           
         foreach my $t (@ALTalleleList) {		# to eliminate an error of a single ALT allele which is the same than REF allele (at a position)
            if (@ALTalleleList == 1) {
               if ($REFallele eq $t) {
                  $validation = 0;
               }
            }
         }
   #############
   
   
         if ($validation == 1) {
            if (exists $record{$scaffold}{$position}{isolates}) {
               $record{$scaffold}{$position}{isolates} = $record{$scaffold}{$position}{isolates} . " " . $SampleName;
            } else {
               $record{$scaffold}{$position}{isolates} = $SampleName;
               $record{$scaffold}{$position}{Ref} = substr($REFallele,0,1);
            }
         }
      }
   }
   close(IN2);
}


warn "\n","*********************","\n";

#####  Write file  ##########################

warn "---Write file :\n\n";

print OUT1 "##scaffold\tRADname\tPosition\tSeqID\tCluster\tRepeat\tCoding\tRef\tIsolates\n";

foreach my $ScaffoldRec ( nsort (keys %record)) {
   foreach my $position (nsort (keys %{$record{$ScaffoldRec}})) {
      if ($record{$ScaffoldRec}{$position}{isolates} ne "") {
         
         print OUT1 "$ScaffoldRec\tna\t$position\t", "na\tna\tna\tna", "\t", $record{$ScaffoldRec}{$position}{Ref}, "\t";
         
         $str_ordered_isolates = "";
         @ordered_isolates = split(" ", $record{$ScaffoldRec}{$position}{isolates});
         foreach my $t (nsort @ordered_isolates) {
            $str_ordered_isolates .= $t . " ";
         }
         chop $str_ordered_isolates;
         print OUT1 $str_ordered_isolates , "\n";
      }
   }

}


