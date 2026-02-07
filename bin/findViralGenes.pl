#!/usr/bin/perl -w
use strict;
use Getopt::Long;

############################################
## AUTHOR:       Eric C. Rouchka          ##
##               University of Louisville ##
## LAST UPDATED: 6/5/2024                 ##
## V8: Added code to save gene sequences  ##
## V6:                                    ##
##                                        ##
## V5: Added 0,1,2,3 designation          ##
##     Added length histograms            ##
## V4: Added new integrity descriptors    ##
##     Added IPDA probes from PMC9525950  ##
############################################

########################################################
## GLOBAL VARIALBLES -- SHOULD CHANGE ON INSTALLATION ##
########################################################
my $DB_BASE   = "/home/kyinbre/data/HIV_GENES/BLAST_DB/";
my $BLAST_CMD = "/home/kyinbre/bin/ncbi-blast-2.10.0+/bin/blastn";
my $TMPDIR    = "/home/kyinbre/tmp/";
my $IPDA_BLASTDB = "/home/kyinbre/data/HIV_GENES/BLAST_DB/HIV/IPDA/IPDA";
my $IPDA_V2_BLASTDB = "/home/kyinbre/data/HIV_GENES/BLAST_DB/HIV/IPDA/IPDA_PMC9525950";
my $PCT_CUTOFF = 0.7;
my $PCT_CUTOFF_2 = 0.7;
my $MIN_MATCH_LEN = 30; ## Min Len of gene match
my $inFN;
my $OUT_FN;
my $reference;
my $virusType;

GetOptions("reference=s" => \$reference,
           "virus=s" => \$virusType,
           "in=s" => \$inFN, "out=s" => \$OUT_FN) || printUsage();
my $blastDB = validateParameters($inFN, $OUT_FN, $reference, $virusType);
my $blastFA = $blastDB . ".fa";
my %refGeneSeqs = getReferenceSequences($blastFA);
my %geneLengthHASH = getGeneLengths($reference, $virusType);
my $LTR_LEN = $geneLengthHASH{"LTR"};
my $numSeqs = getSequenceCount($inFN);
my @LTRPosArr;
my @LTR5PosArr;
my @LTR3PosArr;

for(my $i = 0; $i < $LTR_LEN; $i++) {
   push(@LTRPosArr, 0);
   push(@LTR5PosArr, 0);
   push(@LTR3PosArr, 0);
}
## print "$numSeqs Sequences to Process\n";

open(INFILE, $inFN) || die("Error opening $inFN for reading");
open(OVERALLOUTFILE, ">$OUT_FN") || die("Error opening $OUT_FN for writing");
print OVERALLOUTFILE "CCS_READ_ID\tSTRAND\tGENE_MATCH_STRING\tMATCH_TYPE";
if($virusType eq "HIV") { 
   print OVERALLOUTFILE "\tIPDA_PSI\tIPDA_RRE\tIPDA_INTACT\tIPDA_LTR_GAG\tIPDA_POL\tIPDA_ENV\tIPDA_V2_INTACT";
}
print OVERALLOUTFILE "\tEPISOME_FLAG\n";

my $seqNum = 0;
my $numFull = 0;
my $num5end = 0;
my $num3end = 0;
my $numTrunc = 0;
my $numOther = 0;
my $numFullMinus5LTR = 0;
my $numFullMinus3LTR = 0;
my $numFullMinusBOTHLTR = 0;
my $numInternalDel = 0;
my %matchHASHCNT;
my $numFull5Trunc = 0;
my $numFull3Trunc = 0;
my $numFullBothTrunc = 0;
my $numPutative = 0;
my $numIndeterminate = 0;
my $numTruncated = 0;
my $numHeavilyTruncated = 0;
my $minLTRLen = 100000;
my %LTR5Lens;
my %LTR3Lens;
my $numBins = int($LTR_LEN / 25) + 1;
for(my $i = 0; $i < $numBins; $i++) {
   $LTR5Lens{$i} = 0;
   $LTR3Lens{$i} = 0;
}

my %geneMatchHASH;

while(defined(my $hdr = <INFILE>)) { 
   $seqNum++;
   chomp($hdr);
   my $seq = <INFILE>;
   chomp($seq);
   $seq =~ s/\s+//g;

   my $firstMatchLoc = 1000000;
   my $lastMatchLoc = -1;

   my $tmpFN = getRandomFileName(".fa");
   open(OUTFILE, ">$tmpFN") || die;
   print OUTFILE "$hdr\n$seq\n";
   close(OUTFILE);

   ## Blast agains HIV/SIV genes
   my $tmpBlastFN = getRandomFileName(".blast");
   my $cmd = "$BLAST_CMD -db $blastDB -query $tmpFN -word_size 12 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 -outfmt \"6 qseqid qstart qend sseqid sstrand sstart send evalue bitscore score length pident gaps\" > $tmpBlastFN";
   system($cmd);

   ## Blast against IPDA regions
   my $PSI_Value = 0;
   my $RRE_Value = 0;
   my $LTR_GAG_Value = 0;
   my $POL_Value = 0;
   my $ENV_Value = 0;
 
   if($virusType eq "HIV") { 
      ($PSI_Value, $RRE_Value) = findIPDA($tmpFN);
      ($LTR_GAG_Value, $POL_Value, $ENV_Value) = findIPDA_V2($tmpFN);
   }
 
   $cmd = "rm $tmpFN";
   system($cmd);
   open(BLASTFILE, "$tmpBlastFN") || die;
   my %matchHASH;
   my %matchHASHLTR;
   my %matchGeneLenHASH;

   my $currSeqID = $hdr;
   $currSeqID =~ s/\>//g;
   my $numPlus=0;
   my $numMinus=0;

   while(defined(my $line = <BLASTFILE>)) { 
      chomp($line);
      my @wds = split(/\t/, $line);
      my @wds2 = split(/\_/, $wds[3]);
      my $currGene = $wds2[0];
      my $matchBeg = $wds[1];
      my $matchEnd = $wds[2];
      my $matchLen = $matchBeg - $matchEnd;
      if($matchLen < 0) { $matchLen *= -1; }
      $matchLen += 1;
      my $strand = $wds[4];
      if(defined($matchHASHLTR{$wds[3]})) { 
      }
      if($matchLen >= $MIN_MATCH_LEN) { 
         my $matchSeq = "";
         if($matchBeg < $matchEnd) { $matchSeq = substr($seq, ($matchBeg - 1), $matchLen); }
         else                      { $matchSeq = substr($seq, ($matchEnd - 1), $matchLen); }

         if($strand eq "plus") { 
            $numPlus++;
         }
         else {
            $numMinus++;
            $matchSeq = reverseComplement($matchSeq);
         }


         if(defined($matchHASH{$currGene})) { 
            if((!($matchHASH{$currGene} =~ /\/\/\//)) && ($currGene eq "LTR")) { 
               if(!defined($geneMatchHASH{$currGene})) {
                  $geneMatchHASH{$currGene} = "$hdr\n$matchSeq\n";
               }
               else {
                  $geneMatchHASH{$currGene} .= "$hdr\n$matchSeq\n";
               } 
            }
            my $modFLAG = 0;
            my @mARR = split(/\/\/\//, $matchHASH{$currGene});
            my @mLenARR = split(/\/\/\//, $matchGeneLenHASH{$currGene});
            my $numM = @mARR;
            for(my $mm = 0; $mm < $numM; $mm++) { 
               my($oMBeg, $oMEnd) = split(/\.\./, $mARR[$mm]);
               if(($matchBeg < $oMBeg) && ($matchEnd < $oMBeg)) { 
                  ## Matches are distinct
                  my $mDiff = $oMBeg - $matchEnd; 
                  if($mDiff < 100) { 
                     ## Merge Blast hits within 100 bp
                     my $newMatchStr = $matchBeg . ".." . $oMEnd;
                     my $newLen = $oMEnd - $matchBeg + 1;
                     $mARR[$mm] = $newMatchStr;
                     $mLenARR[$mm] = $newLen;
                     $modFLAG = 1;
                  }
               }
               else {
                  if(($matchBeg > $oMEnd) && ($matchEnd > $oMEnd)) {
                     ## Matches are distinct
                     my $mDiff = $matchBeg - $oMEnd;
                     if($mDiff < 100) { 
                        ## merge blast hits within 100 bp
                        my $newMatchStr = $oMBeg . ".." . $matchEnd;
                        my $newLen = $matchEnd - $oMBeg;
                        $mARR[$mm] = $newMatchStr;
                        $mLenARR[$mm] = $newLen;
                        $modFLAG = 1;
                     } 
                  }
                  else {
                     ## Match is internal -- do not merge
                 }
               }
            }
            if($modFLAG) { 
              if($numM == 1) { 
                 $matchHASH{$currGene} = $mARR[0];
                 $matchGeneLenHASH{$currGene} = $mLenARR[0];
              }
              else {
                 $matchHASH{$currGene} = join(@mARR, "///");
                 $matchGeneLenHASH{$currGene} = join(@mLenARR, "///");
              }
            }
            else {
               $matchHASH{$currGene} .= "///" . $matchBeg . ".." . $matchEnd;
               $matchGeneLenHASH{$currGene} .= "///" . $matchLen;
            }
         }
         else {
            $matchHASH{$currGene} = $matchBeg . ".." . $matchEnd;
            $matchGeneLenHASH{$currGene} = $matchLen;
            if(!defined($geneMatchHASH{$currGene})) {
               $geneMatchHASH{$currGene} = "$hdr\n$matchSeq\n";
            }
            else {
               $geneMatchHASH{$currGene} .= "$hdr\n$matchSeq\n";
            } 
         }
         if($currGene eq "LTR") { 
            my $LTR_Beg = $wds[5] - 1;
            my $LTR_End = $wds[6] - 1;
            if($LTR_Beg > $LTR_End) { 
               my $tmp = $LTR_Beg;
               $LTR_Beg = $LTR_End;
               $LTR_End = $tmp;
            }
            my $currLTRLen = $LTR_End- $LTR_Beg + 1;
            if($currLTRLen > 12) { 
               if($currLTRLen < $minLTRLen) { 
                  $minLTRLen = $currLTRLen;
               }
               for(my $ltrPos = $LTR_Beg; $ltrPos <= $LTR_End; $ltrPos++) { 
                  $LTRPosArr[$ltrPos]++;
               }
            }
         }
         if($matchBeg < $firstMatchLoc) { $firstMatchLoc = $matchBeg; }
         if($matchBeg > $lastMatchLoc)  { $lastMatchLoc = $matchBeg; }
         if($matchEnd < $firstMatchLoc) { $firstMatchLoc = $matchEnd; }
         if($matchEnd > $lastMatchLoc)  { $lastMatchLoc = $matchEnd; }
      }
   }
   close(BLASTFILE);

   ###########################################################
   ## Now test for consistency in gene positions based on   ##
   ## Matches and remove inconsistent genes (typically NEF) ##
   ###########################################################
   my @geneList = ("LTR", "GAG", "POL", "VIF", "VPR", "VPU", "ENV", "NEF", "LTR");
   if($virusType eq "SIV") { 
      @geneList = ("LTR", "GAG", "POL", "VIF", "VPR", "VPX", "ENV", "NEF", "LTR");
   }
   my %genePos;

   my $overallStrand = getOverallStrand($numMinus, $numPlus);
   my $numG = @geneList;
   for(my $g = 0; $g < $numG; $g++) { 
      my $cMatch = "NA";
      my $cLen = "NA";
      $genePos{$geneList[$g]} = -1;
      if(defined($matchHASH{$geneList[$g]})) { 
         $cMatch = $matchHASH{$geneList[$g]};
         my @wdsabc = split(/\/\/\//, $cMatch);
         my @wds123 = split(/\.\./, $wdsabc[0]);
         my $firstPos = $wds123[0];
         if($overallStrand eq "plus") {  
            $genePos{$geneList[$g]} = $wds123[0];
         } 
         else {
            $genePos{$geneList[$g]} = length($seq) - $wds123[1];
         }
         
         $cLen = $matchGeneLenHASH{$geneList[$g]};
      }
#       print "gene:\t$geneList[$g]\tmatch: $cMatch\tlen: $cLen\n";
   }
#    print "POSITIONS\n";
#    for(my $g = 0; $g < $numG; $g++) { 
#       print "$genePos{$geneList[$g]} ";
#    }
#    print "\n";
   my $lastFound = -1;

   # PRINT TO SHOW GENE LOCATIONS
#   for(my $g = 0; $g < ($numG - 1); $g++) {
#      print "$geneList[$g]\t";
#   }
#   print "\n";
#   for(my $g = 0; $g < ($numG - 1); $g++) {
#      print "$genePos{$geneList[$g]}\t";
#   }
#   print "\n";

   my $numOutOfOrder = 0;
   my $geneOut = -1;

   for(my $g = 1; $g < ($numG - 1); $g++) { 
      if($genePos{$geneList[$g]} != -1) { 
         if($genePos{$geneList[$g]} < $lastFound) { 
   #          print "***** ERROR AT $geneList[$g] OUT OF ORDER *****\n";
            if($geneList[$g] eq "NEF") { 
               delete $matchHASH{"NEF"};
               delete $matchGeneLenHASH{"NEF"};
               delete $geneMatchHASH{"NEF"}; 
            }
            else {
               $numOutOfOrder++;
               $geneOut = $g;
               $lastFound = $genePos{$geneList[$g]};
            }
         }
         else {
            if(($lastFound == -1) && ($geneList[$g] eq "NEF") && ($genePos{"LTR"} != -1)) {
               delete $matchHASH{"NEF"};
               delete $matchGeneLenHASH{"NEF"};
               delete $geneMatchHASH{"NEF"}; 
            }
            else {
               $lastFound = $genePos{$geneList[$g]};
            }
         }
      }
   }
   my $episomeFlag = "";
   if($numOutOfOrder == 1) { 
      $episomeFlag = "EPISOME BREAK IN " . $geneList[$geneOut - 1] . "-" . $geneList[$geneOut];
   }
            
   ## NOW ADD IN THE PORTIONS TO SSEARCH FOR THE PSI AND RRE REGIONS AS WELL
   
   removeTmpFile($tmpBlastFN);
   $overallStrand = getOverallStrand($numMinus, $numPlus);
   my($Trunc5FLAG, $Trunc3FLAG) = getTruncationFlags($firstMatchLoc, $lastMatchLoc, $seq, $overallStrand);
   my($minGene, $maxGene, $minLen, $maxLen) = getMinMaxGene(\%matchHASH, $overallStrand);
#   print "MIN: $minGene ($minLen) MAX: $maxGene ($maxLen)\n";

   my($mHashREF, $mglhREF) = assignLTRMatches(\%matchHASH, $minGene, $maxGene, $minLen, $maxLen, \%matchGeneLenHASH, $episomeFlag);
   %matchHASH = %$mHashREF;
   %matchGeneLenHASH = %$mglhREF;

   my($matchSTR, $matchHASHCNTREF, $LTR5Len, $LTR3Len) = getMatchString(\%matchHASH, \%matchHASHCNT, \%matchGeneLenHASH, \%geneLengthHASH, $Trunc5FLAG, $Trunc3FLAG);
   %matchHASHCNT = %$matchHASHCNTREF;
   if($matchSTR =~ /[123]-{8}/) { 
      if(!($Trunc5FLAG && $Trunc3FLAG)) { 
         my $tmpVal = substr($matchSTR, 0, 1);
         $matchHASHCNT{$matchSTR}--;
         if($matchHASHCNT{$matchSTR} == 0) { 
            delete $matchHASHCNT{$matchSTR};
         }
         $matchSTR = $tmpVal . "00000000";
         if(defined($matchHASHCNT{$matchSTR})) {
            $matchHASHCNT{$matchSTR}++;
         }
         else {
            $matchHASHCNT{$matchSTR} = 1;
         }
      }
   }
   if($matchSTR =~ /-{8}[123]/) { 
      if(!($Trunc5FLAG && $Trunc3FLAG)) { 
         my $tmpVal = substr($matchSTR, 8, 1);
         $matchHASHCNT{$matchSTR}--;
         if($matchHASHCNT{$matchSTR} == 0) { 
            delete $matchHASHCNT{$matchSTR};
         }
         $matchSTR = "00000000" . $tmpVal;
         if(defined($matchHASHCNT{$matchSTR})) {
            $matchHASHCNT{$matchSTR}++;
         }
         else {
            $matchHASHCNT{$matchSTR} = 1;
         }
      }
   }


   if($matchSTR =~ /[123]0{7}[123]/) { 
      my $tmpVal = substr($matchSTR, 0, 1);
      $matchHASHCNT{$matchSTR}--;
      if($matchHASHCNT{$matchSTR} == 0) { 
         delete $matchHASHCNT{$matchSTR};
      }
      $matchSTR = $tmpVal . "00000000";
      $Trunc5FLAG = ($Trunc5FLAG && $Trunc3FLAG);
      $Trunc3FLAG = $Trunc5FLAG;
      $matchSTR = updateMatchString($matchSTR, $Trunc5FLAG, $Trunc3FLAG);
      if(defined($matchHASHCNT{$matchSTR})) { 
         $matchHASHCNT{$matchSTR}++;
      }
      else {
         $matchHASHCNT{$matchSTR} = 1;
      }
   }
 
   if($matchSTR =~ /^\-/) { 
   }
   else {
      if($LTR5Len == 0) { 
         $LTR5Lens{0}++;
      }
      else {
         my $currBin = int($LTR5Len / 25) + 1;
         $LTR5Lens{$currBin}++;
      }
   }
   if($matchSTR =~ /\-$/) { 
   }
   else {
      if($LTR3Len == 0) { 
         $LTR3Lens{0}++;
      }
      else {
         my $currBin = int($LTR3Len / 25) + 1;
         $LTR3Lens{$currBin}++;
      }
   }

   my $currType = "";
   $currType = findMatchType($matchSTR);
   if($currType eq "INTACT") { 
      $numFull++;
   }
   else {
      if($currType eq "PUTATIVELY INTACT") { 
         $numPutative++;
      }
      else {
         if($currType eq "INDETERMINATE") { 
            $numIndeterminate++;
         }
         else {
            if($currType eq "HEAVILY TRUNCATED") { 
               $numHeavilyTruncated++;
            }
            else {
               if($currType eq "TRUNCATED") { 
                  $numTruncated++;
               }
               else {
                  if($currType eq "INTERNAL DELETION") { 
                     $numInternalDel++; 
                  }
                  else {
                     print "TYPE: $currType\n";
                  }
               }
            }
         }
      }
   }
   
   if(($seqNum % 100) == 0) { 
  ##    print "$seqNum sequences processed\n";
   }
   print OVERALLOUTFILE "$currSeqID\t$overallStrand\t$matchSTR\t$currType";
   if($virusType eq "HIV") { 
      my $foundFlag = $PSI_Value && $RRE_Value;
      my $found_V2_Flag = 0;
      if(($LTR_GAG_Value >= 3) && ($POL_Value >= 3) && ($ENV_Value >= 3)) { 
         $found_V2_Flag = 1;
      }
      print OVERALLOUTFILE "\t$PSI_Value\t$RRE_Value\t$foundFlag";
      print OVERALLOUTFILE "\t$LTR_GAG_Value\t$POL_Value\t$ENV_Value\t$found_V2_Flag";
   }
   print OVERALLOUTFILE "\t$episomeFlag";
   print OVERALLOUTFILE "\n";
}
print "\n\n";
print "INTACT: $numFull\nPUTATIVELY INTACT: $numPutative\n";
print "INDETERMINATE: $numIndeterminate\nHEAVILY TRUNCATED: $numHeavilyTruncated\n";
print "TRUNCATED: $numTruncated\nINTERNAL DELETION: $numInternalDel\nOTHER: $numOther\n";
close(INFILE);
close(OVERALLOUTFILE);

print "\n\n";

my @matchGeneARR = keys(%geneMatchHASH);
my $numMatchGenes = @matchGeneARR;
for(my $i = 0; $i < $numMatchGenes; $i++) {
   my $matchFN = $OUT_FN . "_" . $matchGeneARR[$i] . "_matches.fa";
   my $alignFN = $OUT_FN . "_" . $matchGeneARR[$i] . "_ALIGNMENT.fa";
   open(MATCHFILE, ">$matchFN") || die("Error opening $matchFN for writing");
   print MATCHFILE "$refGeneSeqs{$matchGeneARR[$i]}";
   print MATCHFILE "$geneMatchHASH{$matchGeneARR[$i]}";
   close(MATCHFILE);
   my $cmd = "muscle -in $matchFN -out $alignFN";
#   system($cmd);
}

printSummary(\%matchHASHCNT, $virusType);
print "MIN LTR LENGTH: $minLTRLen\n";
makeLTRGraph($LTR_LEN, \@LTRPosArr);
makeLTRHistogram(\%LTR5Lens, \%LTR3Lens);

exit(0); 
#-----------------------------------------------------------------------------
sub makeLTRHistogram {
   my($LTR5hREF, $LTR3hREF) = @_;
   my %LTR5Lens = %$LTR5hREF;
   my %LTR3Lens = %$LTR3hREF;

   open(LTRHISTOFILE, ">LTRHistogramData.txt") || die("Error opening LTRHistogramData.txt for writing");

   my @k = sort{$a <=> $b}(keys(%LTR5Lens));
   my $numK = @k;
   my $currVal = 0;
   print LTRHISTOFILE "UTR\tBASES\tCOUNT\n";
   for(my $i = 0; $i < $numK; $i++) { 
      $currVal = $i * 25;
      print LTRHISTOFILE "5LTR\t$currVal\t$LTR5Lens{$k[$i]}\n";
   }
   @k = sort{$a <=> $b}(keys(%LTR3Lens));
   $numK = @k;
   for(my $i = 0; $i < $numK; $i++) { 
      $currVal = $i * 25;
      print LTRHISTOFILE "3LTR\t$currVal\t$LTR3Lens{$k[$i]}\n";
   }
   close(LTRHISTOFILE);
   open(RFILE, ">makeLTRHistogram.R") || die("Error opening makeLTRHistogram.R for writing");
   print RFILE "";
   print RFILE "library(\"ggplot2\")\n";
   print RFILE "library(\"cowplot\")\n";
   print RFILE "LTRHistoData <- read.csv(file=\"LTRHistogramData.txt\", header=TRUE, sep=\"\\t\")\n";
   print RFILE "LTR3 <- LTRHistoData[LTRHistoData\$UTR==\"3LTR\", ]\n";
   print RFILE "LTR5 <- LTRHistoData[LTRHistoData\$UTR==\"5LTR\", ]\n";
   print RFILE "p3 <- ggplot(data=LTR3, aes(x=BASES, y=COUNT)) +\n";
   print RFILE "             geom_bar(color=\"blue\", fill=\"dodgerblue\", stat=\"identity\") +\n";
   print RFILE "             xlab(\"Number of Bases\")+\n";
   print RFILE "             ylab(\"Count\")+\n";
   print RFILE "             ggtitle(\"3\\' LTR Lengths\")\n";
   print RFILE "             theme(plot.title=element_text(hjust=0.5))\n";
   print RFILE "ggsave(plot=p3, file=\"LTR3_Histogram.png\", width=10, height=5, units=\"in\")\n";
   print RFILE "\n\n";
   print RFILE "p5 <- ggplot(data=LTR5, aes(x=BASES, y=COUNT)) +\n";
   print RFILE "             geom_bar(color=\"red\", fill=\"orange\", stat=\"identity\") +\n";
   print RFILE "             xlab(\"Number of Bases\")+\n";
   print RFILE "             ylab(\"Count\")+\n";
   print RFILE "             ggtitle(\"5\\' LTR Lengths\")+\n";
   print RFILE "             theme(plot.title=element_text(hjust=0.5))\n";
   print RFILE "ggsave(plot=p5, file=\"LTR5_Histogram.png\", width=10, height=5, units=\"in\")\n";
   print RFILE "p <- plot_grid(p5, p3, nrow=2, ncol=1)\n";
   print RFILE "ggsave(plot=p, file=\"LTRBOTH_Histogram.png\", width=10, height=10, units=\"in\")\n";

   close(RFILE);
   my $cmd = "/home/kyinbre/anaconda3/envs/R/bin/R CMD BATCH makeLTRHistogram.R";
   system($cmd);
   $cmd = "chmod g+w LTRHistogramData.txt";
   system($cmd);
   $cmd = "chmod g+w makeLTRHistogram.R";
   system($cmd);

}
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
sub makeLTRGraph {
   my($LTR_LEN, $aREF) = @_;

   my @LTRPosArr = @$aREF;

   open(LTRCNTFILE, ">LTRPosCnts.txt") || die("Error opening LTRPosCnts.txt for writing");
   for(my $i = 0; $i < $LTR_LEN; $i++) { 
      print LTRCNTFILE "$i\t$LTRPosArr[$i]\n";
   }
   close(LTRCNTFILE);

   open(RFILE, ">makeLTRGraphic.R") || die("Error opening makeLTRGraphic.R for writing");
   print RFILE "library(\"ggplot2\")\n";
   print RFILE "LTRData <- read.csv(file=\"LTRPosCnts.txt\", header=FALSE, sep=\"\\t\")\n";
   print RFILE "colnames(LTRData) <- c(\"X\", \"Y\")\n";
   print RFILE "p <- ggplot(data=LTRData, aes(x=X, y=Y)) +\n";
   print RFILE "            geom_line() +\n";
   print RFILE "            xlab(\"LTR Position\") +\n";
   print RFILE "            ylab(\"Count\")\n";
   print RFILE "ggsave(plot=p, file=\"LTRUsage.png\", width=10, height=5, units=\"in\")\n";
   print RFILE "sessionInfo()\n";
   close(RFILE);
   my $cmd = "/home/kyinbre/anaconda3/envs/R/bin/R CMD BATCH makeLTRGraphic.R";
   system($cmd);
   $cmd = "chmod g+w LTRPosCnts.txt";
   system($cmd);
   $cmd = "chmod g+w makeLTRGraphic.R";
   system($cmd);
}

#-----------------------------------------------------------------------------
sub printSummary {
    my($hREF, $vType) = @_;
    my %matchHASHCNT = %$hREF;

   my @k = sort(keys(%matchHASHCNT));
   my $numK = @k;
   for(my $i = 0; $i < $numK; $i++) { 
      print "$k[$i]\t$matchHASHCNT{$k[$i]}\n";
   }

   print "\n\n";
   print "5               3\n";
   print "'               '\n";
   print "L G P V V V E N L\n";
   print "T A O I P P N E T\n";
   if($vType eq "HIV") { 
      print "R G L F R U V F R\n";
   }
   else {
      print "R G L F R X V F R\n";
   }
   print "=========================\n";
   for(my $i = 0; $i < $numK; $i++) { 
      my $currK = $k[$i];
      my $currCnt = $matchHASHCNT{$k[$i]};
      my $l = length($currK);

      for(my $j = 0; $j < $l; $j++) { 
         my $currCh = substr($currK, $j, 1);
         print "$currCh ";
      }
      print "\t$currCnt\n";
   }
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getRandomFileName {
   my($extension) = @_;

   my $randFN = $TMPDIR . "/" . time() . "_" . rand() . "$extension";
   while(-e $randFN) { 
      $randFN = $TMPDIR . "/" . time() . "_" . rand() . "$extension";
   }
   return($randFN);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub reverseComplement {
   my($s) = @_;

   $s = uc($s);
   $s =~ s/A/t/g;
   $s =~ s/T/a/g;
   $s =~ s/C/g/g;
   $s =~ s/G/c/g;
   $s = uc($s);
   $s = reverse($s);
   return($s);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub printUsage {
   print "USAGE: -in <IN_FILE> -out <OUT_FILE>\n";
   exit(0);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getSequenceCount {
   my($in_FN) = @_;

   my $tmpCntFN = getRandomFileName(".txt");

   my $cmd = "grep -c \">\" $in_FN > $tmpCntFN";
   system($cmd);
   open(INFILE, "$tmpCntFN") || die;
   my $numSeqs = <INFILE>;
   chomp($numSeqs);
   close(INFILE);
   $cmd = "rm $tmpCntFN";
   system($cmd);
   return($numSeqs);
}
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
sub updateMatchString {
   my($s, $begFLAG, $endFLAG) = @_;
   my $l = length($s);
   my $origS = $s;

   my $begZeros = 0;
   my $endZeros = 0;
   if($begFLAG) {
      my $begReplace = "";

      for(my $i = 0; $i < $l; $i++) {
         my $currCH = substr($s, $i, 1);
         if($currCH eq "0") {
            $begReplace .= "-";
         }
         else {
            $i=$l;
         }
      }
      $s =~ s/^0+/$begReplace/;
   }
   if($endFLAG) {
      my $endReplace = "";
      for(my $i = ($l - 1); $i >= 0; $i--) {
         my $currCH = substr($s, $i, 1);
         if($currCH eq "0") {
            $endReplace .= "-";
         }
         else {
            $i=-1;
         }
      }
      $s =~ s/0+$/$endReplace/;
   }
   return($s);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getMatchString {
   my($mhREF, $mhcREF, $mglhREF, $glhREF, $Trunc5FLAG, $Trunc3FLAG) = @_;

   my %matchHASH = %$mhREF;
   my %matchHASHCNT = %$mhcREF;
   my %matchGeneLenHASH = %$mglhREF;
   my %geneLenghHASH = %$glhREF;

   my @testARR = ("5LTR", "GAG", "POL", "VIF", "VPR", "VPU", "ENV", "NEF", "3LTR");
   if($virusType eq "SIV") {
      $testARR[5] = "VPX";
   }
   my $matchSTR = "";
   my $numToTest = @testARR;
   for(my $i = 0; $i < $numToTest; $i++) { 
      my $currGene = $testARR[$i];
      my $val = 0; 
      if(defined($matchHASH{$currGene})) { 
         $val = 1;
         my $gLen = 1;
         my $testGene = $currGene;

         if($currGene =~ /LTR/) { 
            $testGene = "LTR";
         }
         my $testVal = $geneLengthHASH{$testGene};
         if(defined($testVal)) { 
            my @wds123 = split(/\/\/\//, $testVal);
            $gLen = $wds123[0];
         }
         else {
          #  print "UNDEFINED LENGTH FOR $currGene\n"; 
         }
         my $matchGeneLen = 0;
         $testVal = $matchGeneLenHASH{$currGene};
         if(defined($testVal)) { 
            my @wds123 = split(/\/\/\//, $testVal);
            $matchGeneLen = $wds123[0];
         }
        
         my $pctFound = $matchGeneLen / $gLen;
         # print "$currGene PCT FOUND: $pctFound\n";
         if($pctFound > $PCT_CUTOFF_2) { 
            $val = 2;
         }
         if($pctFound > $PCT_CUTOFF) { 
            $val = 3;
         }
         if($currGene =~ /LTR/) { 
         }
      }
      $matchSTR .= $val; 
   }
   $matchSTR = updateMatchString($matchSTR, $Trunc5FLAG, $Trunc3FLAG);
   if(!defined($matchHASHCNT{$matchSTR})) { 
      $matchHASHCNT{$matchSTR} = 1;
   }
   else {
      $matchHASHCNT{$matchSTR}++;
   }
   my $LTR5Len = $matchGeneLenHASH{"5LTR"};
   my $LTR3Len = $matchGeneLenHASH{"3LTR"};
   if(!defined($LTR5Len)) { $LTR5Len = 0; }
   if(!defined($LTR3Len)) { $LTR3Len = 0; }

   return($matchSTR, \%matchHASHCNT, $LTR5Len, $LTR3Len);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getGeneLengths {
   my($reference, $virusType) = @_;
   my @testGeneARR = ("LTR", "GAG", "POL", "VIF", "VPR", "VPU", "ENV", "NEF");
   my $numGenes = @testGeneARR;
   my %lHASH;
   for(my $i = 0; $i < $numGenes; $i++) { 
      my $currGene = $testGeneARR[$i];
      my $currFN = $DB_BASE . "/$virusType/$reference/$currGene" . ".fa";
      open(GENEFILE, $currFN) || die("Cannot open $currFN for reading");
      my $hdr = <GENEFILE>;
      my $seq = "";
      while(defined(my $line = <GENEFILE>)) { 
         chomp($line);
         $seq .= $line;
      }
      close(GENEFILE);
      $seq =~ s/\s+//g;
      my $geneLen = length($seq);
      $lHASH{$currGene} = $geneLen; 
   }
   return(%lHASH);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

sub validateParameters {
   my($inFN, $OUT_FN, $reference, $virusType) = @_;

   if((!defined($inFN)) || (!defined($OUT_FN)) || (!defined($reference)) || (!defined($virusType))){
      printUsage();
   }
   if(!(-e $inFN)) {
      die("Input file $inFN does not exist");
   }

   if(!(($virusType eq "HIV") || ($virusType eq "SIV"))) {
      die("virus must be either HIV or SIV");
   }
   my $blastDB = $DB_BASE . "/$virusType/$reference/allGenes";
   my $blastTest = $blastDB . ".fa";
   if(!(-e $blastTest)) {
      die("Database $blastDB does not exist");
   }
   return($blastDB);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getTruncationFlags {
   my($firstMatchLoc, $lastMatchLoc, $seq, $strand) = @_;

   my $l = length($seq);

   my $Trunc5FLAG = 0;
   my $Trunc3FLAG = 0;

   if($firstMatchLoc < 100) {
      if($strand eq "plus") {
         $Trunc5FLAG = 1;
      }
      else {
         $Trunc3FLAG = 1;
      }
   }
   if($lastMatchLoc > ((length($seq) - 100))) {
      if($strand eq "plus") {
         $Trunc3FLAG = 1;
      }
      else {
         $Trunc5FLAG = 1;
      }
   }
   return($Trunc5FLAG, $Trunc3FLAG);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getMinMaxGene {
   my($hREF, $overallStrand) = @_;

   my %matchHASH = %$hREF;
   my $min = 100000;
   my $max = -100000;
   my $minGene = "";
   my $maxGene = "";
   my $minLen = 0;
   my $maxLen = 0;
   my @k = sort(keys(%matchHASH));
   if(defined($matchHASH{"LTR"})) { 
      push(@k, "LTR");  ## Make LTR last one added 
   }
   my $numK = @k;
   for(my $i = 0; $i < $numK; $i++) {
      my @wds22 = split(/\/\/\//, $matchHASH{$k[$i]});
      my $numWds2 = @wds22;
      if($k[$i] ne "LTR") { 
         $numWds2 = 1; 
         ## Only look at best hit for genes besides LTR
      }
      for(my $j = 0; $j < $numWds2; $j++) { 
         my $cVal = $wds22[$j]; 
         #print "GENE: $k[$i] CVAL: $cVal\n";
         my @wds = split(/\.\./, $cVal);
         my $pos = $wds[0];
         if($pos > $max) {
            $max = $pos;
            $maxGene = $k[$i];
            $maxLen = $wds[1] - $wds[0];
            if($maxLen < 0) { $maxLen *= -1; }
            $maxLen += 1;
         }
         if($pos <= $min) {
           # if($pos == $min) {
               if(($minGene ne "LTR") || ($k[$i] eq "LTR")) { 
                  $min = $pos;
                  $minGene = $k[$i];
                  $minLen = $wds[1] - $wds[0];
                  if($minLen < 0) { $minLen *= -1; }
                  $minLen += 1;
               }
           # }
         }
      }
   }
   if($overallStrand eq "minus") {
      my $tmp = $minGene;
      $minGene = $maxGene;
      $maxGene = $tmp;
      $tmp = $minLen;
      $minLen = $maxLen;
      $maxLen = $tmp;
      $tmp = $min;
      $min = $max;
      $max = $tmp;
   }
   return($minGene, $maxGene, $minLen, $maxLen);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub assignLTRMatches {
   my($hREF, $minGene, $maxGene, $minLen, $maxLen, $mglhREF, $eFlag) = @_;

   my %matchHASH = %$hREF;
   my %matchGeneLenHASH = %$mglhREF;

   if(defined($matchHASH{"LTR"})) {
      ## print "LTR MATCH FOUND!!! (MIN $minGene; MAX $maxGene)\n";
      my $val = $matchHASH{"LTR"};
      if($matchHASH{"LTR"} =~ /\/\/\//) {
         if($minGene eq "LTR") { 
            $matchHASH{"5LTR"} = 1;
            $matchGeneLenHASH{"5LTR"} = $minLen;
         }
         if($maxGene eq "LTR") { 
            $matchHASH{"3LTR"} = 1;
            $matchGeneLenHASH{"3LTR"} = $maxLen;
         }
         if($eFlag ne "") {
            ## EPISOMAL WITH TWO LTR MATCHES
            my $tmpM = $matchHASH{"LTR"};
            my @tmpA = split(/\/\/\//, $matchHASH{"LTR"});
            my @tmpB = split(/\/\/\//, $matchGeneLenHASH{"LTR"});
            $matchHASH{"5LTR"} = 1;
            $matchHASH{"3LTR"} = 1;
            if($tmpA[0] < $tmpA[1]) { 
               $matchGeneLenHASH{"5LTR"} = $tmpB[0];
               $matchGeneLenHASH{"3LTR"} = $tmpB[1];
            }
            else {
               $matchGeneLenHASH{"5LTR"} = $tmpB[1];
               $matchGeneLenHASH{"3LTR"} = $tmpB[0];
            }
         }
        # $matchHASH{"3LTR"} = 1;
        # $matchHASH{"5LTR"} = 1;
      }
      else {
         if($eFlag ne "") {
            ## EPISOME -- ASSIGN TO 5LTR
            $matchHASH{"5LTR"} = 1;
            $matchGeneLenHASH{"5LTR"} = $matchGeneLenHASH{"LTR"};
         }
         if($minGene eq "LTR") {
            $matchHASH{"5LTR"} = 1;
            $matchGeneLenHASH{"5LTR"} = $minLen;
         }
         else {
            if($maxGene eq "LTR") {
               $matchHASH{"3LTR"} = 1;
               $matchGeneLenHASH{"3LTR"} = $maxLen;
            }
            else {
   #            print "*+*+*+++** LTR MATCHES IN MIDDLE _*_*_*__*_*\n";
            }
         }
      }
   }
   return(\%matchHASH, \%matchGeneLenHASH);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getOverallStrand {
   my($numMinus, $numPlus) = @_;

   my $overallStrand = "plus";
   if($numMinus > $numPlus) { 
      $overallStrand = "minus";
   }
   return($overallStrand);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub removeTmpFile {
   my($fn) = @_;

   my $cmd = "rm $fn";
   system($cmd);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub findIPDA_V2 {
   my($tmpFN) = @_;

   my $LTR_GAG_Value = 0;
   my $POL_Value = 0;
   my $ENV_Value = 0;

   my $tmpIPDABlastFN = getRandomFileName(".IPDA_V2.blast");
   my $cmd = "$BLAST_CMD -db $IPDA_V2_BLASTDB -query $tmpFN -word_size 9 -outfmt \"6 qseqid qstart qend sseqid sstrand sstart send evalue bitscore score length pident gaps\" > $tmpIPDABlastFN";
   system($cmd);
   open(IPDAFILE, $tmpIPDABlastFN) || die("Error opening $tmpIPDABlastFN for reading");
   while(defined(my $IPDAblast = <IPDAFILE>)) { 
      chomp($IPDAblast);
      my @wds = split(/\t/, $IPDAblast);
      my $region = $wds[3];
      if($region eq "LTR_FWD") { 
         if($LTR_GAG_Value == 0) { 
            $LTR_GAG_Value = 1; 
         }
         else {
            if($LTR_GAG_Value == 2) { 
               $LTR_GAG_Value = 3;
            }
         }
      }
      if($region eq "LTR_REV") { 
         if($LTR_GAG_Value == 0) { 
            $LTR_GAG_Value = 2;
         }
         else {
            if($LTR_GAG_Value == 1) { 
               $LTR_GAG_Value = 3;
            }
         }
      }
      if($region eq "POL_FWD") { 
         if($POL_Value == 0) { 
            $POL_Value = 1; 
         }
         else {
            if($POL_Value == 2) { 
               $POL_Value = 3;
            }
         }
      }
      if($region eq "POL_REV") { 
         if($POL_Value == 0) { 
            $POL_Value = 2;
         }
         else {
            if($POL_Value == 1) { 
               $POL_Value = 3;
            }
         }
      }
      if($region eq "ENV_FWD") { 
         if($ENV_Value == 0) { 
            $ENV_Value = 1; 
         }
         else {
            if($ENV_Value == 2) { 
               $ENV_Value = 3;
            }
         }
      }
      if($region eq "ENV_REV") { 
         if($ENV_Value == 0) { 
            $ENV_Value = 2;
         }
         else {
            if($ENV_Value == 1) { 
               $ENV_Value = 3;
            }
         }
      }
   }
   close(IPDAFILE); 
   $cmd = "rm $tmpIPDABlastFN";
   system($cmd);

   return($LTR_GAG_Value, $POL_Value, $ENV_Value);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub findIPDA {
   my($tmpFN) = @_;
   my $PSI_Value = 0;
   my $RRE_Value = 0;

   ## SEE IF PSI AND RRE REGIONS ARE FOUND ##
   my $tmpIPDABlastFN = getRandomFileName(".IPDA.blast");
   my $cmd = "$BLAST_CMD -db $IPDA_BLASTDB -query $tmpFN -word_size 11 -outfmt \"6 qseqid qstart qend sseqid sstrand sstart send evalue bitscore score length pident gaps\" > $tmpIPDABlastFN";
   system($cmd);
   open(IPDAFILE, $tmpIPDABlastFN) || die("Error opening $tmpIPDABlastFN for reading");
   while(defined(my $IPDAblast = <IPDAFILE>)) { 
      chomp($IPDAblast);
      my @wds = split(/\t/, $IPDAblast);
      my $region = $wds[3];
      if($region eq "PSI") { $PSI_Value = 1; }
      if($region eq "RRE") { $RRE_Value = 1; }
   }
   close(IPDAFILE); 
   $cmd = "rm $tmpIPDABlastFN";
   system($cmd);
   return($PSI_Value, $RRE_Value);
}
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
sub findMatchType {
   my($p) = @_;
   #$p =~ s/2/1/g;
   #$p =~ s/3/1/g;

   my $mType = "INTERNAL DELETION";
   # if($p =~ /^0+/) { 
   #    $mType = "INDETERMINATE";
  #  }
   # if($p =~ /0+$/) {
    #   $mType = "INDETERMINATE";
   # }
   if($p eq "333333333") { 
      $mType = "INTACT";
   }
   else {
      if(($p eq "33333333-") || ($p eq "-33333333") || ($p eq "-3333333-") ||
         ($p eq "33333332-") || ($p eq "-23333333") || ($p eq "-2333333-") ||
         ($p eq "33333331-") || ($p eq "-13333333") || ($p eq "-1333333-") ||
         ($p eq "-1333332-") || ($p eq "-1333331-") || ($p eq "-2333332-") ||
         ($p eq "-2333331-") || ($p eq "--3333333") || ($p eq "3333333--") ||
         ($p =~ /[0123]3{7}[0123]/)) {
         $mType = "PUTATIVELY INTACT";
      } 
      else {
        if(($p =~ /^0{4,}/) || ($p =~ /0{4,}$/))  {
           $mType = "HEAVILY TRUNCATED";
           $mType = "TRUNCATED";
        }
        else {
           if(($p =~ /^0{1,}[123]3*[123]-+$/) || ($p =~ /^-+[123]3*[123]0{1,}$/)) { 
              $mType = "TRUNCATED";
           }
           else {
              if(($p =~ /^0{1,}3*$/) || ($p =~ /^3*0{1,}$/)) { 
                 $mType = "TRUNCATED";
              }
           
              else {
                if(($p =~ /^-{2,}[123]3*$/)) { 
                   $mType = "INDETERMINATE";
                }
                else {
                   if(($p =~ /^3*[123]-{2,}$/)) { 
                      $mType = "INDETERMINATE";
                   }
                   else {
                      if(($p =~ /^-{2,}[123]3*[123]-{2,}$/)) {
                         $mType = "INDETERMINATE";
                      }
                      else {
                         if($p =~ /^-+[12]0{2,}$/) { 
                            $mType = "INDETERMINATE";
                            $mType = "TRUNCATED";
                         }
                         else {
                            if($p =~ /^-+[12]{2}-+$/) { 
                               $mType = "INDETERMINATE";
                            }
                            else {
                               if($p =~ /^-+[123]3*[123]-+$/) {
                                  $mType = "INDETERMINATE";
                               } 
                               else {
                                  if(($p =~ /^[12]3*[12]-+$/) || ($p =~ /^-+[12]3*[12]$/)) { 
                                     $mType = "INDETERMINATE";
                                  }
                               }
                            }
                         }
                      }
                   }
                }
             }
          }
        }
      }
   }

   return($mType);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getReferenceSequences {
   my($inFN) = @_;
   my %h;
   open(INFILE, $inFN) || die("Error opening $inFN for reading");
   my $lastHDR = "";
   my $seq = "";
   my $lastGene = "";
   while(defined(my $line = <INFILE>)) { 
      chomp($line);
      if($line =~ /^\>/) { 
         if($lastHDR ne "") { 
            if(!defined($h{$lastGene})) { 
               $h{$lastGene} = "$lastHDR\n$seq\n";
            }
            else {
               $h{$lastGene} .= "$lastHDR\n$seq\n";
            }
         }
         $lastHDR = $line;
         my @wds = split(/\_/, $lastHDR);
         $lastGene = $wds[0];
         $lastGene =~ s/^\>//g;
         $seq = "";
      }
      else {
         $seq .= $line;
      }
   }
   if(!defined($h{$lastGene})) { 
      $h{$lastGene} = "$lastHDR\n$seq\n";
   }
   else {
      $h{$lastGene} .= "$lastHDR\n$seq\n";
   }
   close(INFILE);
   return(%h);
}
#-----------------------------------------------------------------------------