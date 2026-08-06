#!/usr/bin/perl -w
use strict;

my $FN = $ARGV[0];
#my %ccsHASH = getCCSIDs($FN);
my %ccsHASH;
while(defined(my $line = <>)) { 
   chomp($line);   
   print "-->$line\n";
   my @wds = split(/\t/, $line);
   my $currCCS = $wds[0];
   print "currCCS $currCCS\n";

   if(defined($ccsHASH{$currCCS})) { 
      print "LINE: $line\n";
   }
}


sub getCCSIDs {
   my($inFN) = @_;
   open(INFILE, $inFN) || die;
   my %h;
   while(defined(my $line = <INFILE>)) {
      chomp($line);
      $h{$line} = $line;
   }
   close(INFILE);
   return(%h);
}


