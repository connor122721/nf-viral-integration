#!/usr/bin/perl -w
use strict;

################################
## Authors: Eric C. Rouchka   ##
##          Ghazal Sadri Mehs ##
## Date:    06/09/2026        ##
## Last Update: 06/26/2026    ##
################################

################################
## Get command line arguments ##
################################
my $numArgs = @ARGV;
if($numArgs != 3) { 
   die("USAGE: parseSAM <bamFN> <HIV_label> <output_prefix>\n");
}
print "==== RUNNING FOR $ARGV[2] ====\n";

######################
## GLOBAL VARIABLES ##
######################
my $MIN_FLANK_LEN = 50;
my $MAX_INTEGRATION_DIFF = 10;
my $MAX_HIV_DIFF = 10;
my $MAX_HIV_POS_DIFF = 10;
my $MAX_SHEAR_DIFF = 10;
my $HIV_CHR_NAME = $ARGV[1];
my $SEPARATOR = "xXxXx";  ## Potentially can be found in quality score
my $SEPARATOR2 = "///";	  ## Potentially can be found in quality score

my ($bamFN, $HIV_label, $OUTPUT_PRE) = @ARGV;
print "retrieving CCS reads that hit HIV (may take some time)\n";
my %HIV_ccsHASH = getHIV_CCSIDs($bamFN, $HIV_label);
print "converting bam to sam for processing (may take some time -- likely several minutes)\n";
my $samFN = createSAMFile($bamFN);
print "parsing sam file for only CCS reads that hit HIV (may take some time)\n";
my $HIVparsedSAM_FN = getHIVSAM($samFN, \%HIV_ccsHASH);


open(my $ERRFILE, ">error_log_$OUTPUT_PRE.txt") || die;
open(my $WARNFILE, ">warning_log_$OUTPUT_PRE.txt") || die;

my $outFN = "$OUTPUT_PRE" . "_MasterOfMasterFrame.tsv";
open(my $OUTFILE, ">$outFN") || die("Error opening $outFN for writing");
print $OUTFILE "READ\tFLANK_TYPE\tPOS\tINTEGRATION_SITE\tSHEAR_SITE\tCLONE_ID\tPCR_ID\tHIV_POS\tHIV_LEN\tHOST_SEQ_OVERLAP\n";

my %HIV_NOHOST_HASH;
my %CCS_NOHOST_HASH;

print $ERRFILE "================\n";
print $ERRFILE "= $OUTPUT_PRE   =\n";
print $ERRFILE "================\n";
print $WARNFILE "================\n";
print $WARNFILE "= $OUTPUT_PRE   =\n";
print $WARNFILE "================\n";
my %CHR_HASH;
my %HIV_HASH;
my %CCS_HASH;
my %FLANK_HASH;
my %POS_HASH;
my %INTEGRATION_HASH;
my %SHEAR_HASH;
my %CLONE_HASH;
my %HIV_LEN_HASH;
my %HIV_POS_HASH;
my %PCR_HASH;
my %CLONE_PCR_HASH;
my %OVERLAP_SEQ_HASH;

my $currBAM = $HIVparsedSAM_FN; 
my @clonePos;
if(!(-e $currBAM)) { 
   print "SKIPPING $currBAM (does not exist)\n";
}
else {
   print "processing sam file\n";
   my($cHashREF, $hHashREF, $ccsHashREF) = parseSAMFile($currBAM, \%CHR_HASH, \%HIV_HASH, \%CCS_HASH, $ERRFILE, $WARNFILE);
   print "sam file processed\n";
   %CHR_HASH = %$cHashREF; 
   %HIV_HASH = %$hHashREF;
   %CCS_HASH = %$ccsHashREF;

   my @ccsIDs = keys(%CCS_HASH);
   my $numCCS = @ccsIDs;
   my $currShearSite = -1;
   my $currIntegrationSite = -1;
   my $currFlankType = "NONE";
   my $currHIVLen = 0;
   my $currHIVPos = "UNK";
   my $currPosVal = -1;
   my $currChrNum = "UNK";
   for(my $c = 0; $c < $numCCS; $c++) { 
      my $currSeqOverlap = "";
      my $currCCSID = $ccsIDs[$c];
      my $currCHRmatch = $CHR_HASH{$currCCSID};
      my $currHIVmatch = $HIV_HASH{$currCCSID};
      if(!defined($currCHRmatch)) { 
         $CCS_NOHOST_HASH{$currCCSID} = $currCCSID;
         $HIV_NOHOST_HASH{$currCCSID} = $currHIVmatch;
      }
      else {
         ## CHECK FOR DUPLICATES
         my @chrHits = split(/$SEPARATOR/, $currCHRmatch);
         my $numHits = @chrHits;
         my @wds123 = split(/\t/, $chrHits[0]);
         $currChrNum = $wds123[2];
         if($numHits > 1) {
            ## Multiple host hits
            if($numHits > 2) { 
               print $ERRFILE "TOO MANY HOST HITS FOR $currCCSID\n";
            }
            else {
               my($shear1, $shear2, $int1, $int2, $overlapBases, $HIVLen) = getOverlappingInformation($currCCSID, $chrHits[0], $chrHits[1], $ERRFILE, $WARNFILE);
               if(($shear1 != -1) && ($shear2 != -1)) {
                  $currShearSite = $shear1 . $SEPARATOR2 . $shear2;
               }
               $currSeqOverlap = $overlapBases;
               $currIntegrationSite = $int1;
               $currFlankType = "Dually Flanked";
               $currHIVLen = $HIVLen;
               $currPosVal = $int1;
            }
         }
         else { 
            ## Single Host Hit
            my($shear, $int, $HIVLen, $currPos) = getSingleHostHitInformation($currCCSID, $chrHits[0], $ERRFILE, $WARNFILE);
            $currShearSite = $shear;
            $currIntegrationSite = $int;
            $currPosVal = $currPos;
            $currHIVLen = $HIVLen;
         }
         if(!defined($currHIVmatch)) { 
            print $WARNFILE "No HIV hit for\t$currCCSID\n";
            if(!($currCHRmatch =~ /$SEPARATOR/)) { 
               ## Only calculate for single host hit 
               my @wds = split(/\t/, $currCHRmatch);
               my $currCigar = $wds[5];
               my $largeLeft = "";
               my $largeRight = "";
               my $largestIns = 0;
               while($currCigar =~ /([0-9]+)I/g) {
                  if ($1 > $largestIns) { 
                     $largestIns = $1;
                     $largeLeft = $`;
                     $largeRight = $';
                  }
               }
               my($leftM, $leftI, $leftD) = findMatchInsDel($largeLeft);
               my($rightM, $rightI, $rightD) = findMatchInsDel($largeRight);
               my $leftLen = $leftM + $leftD;
               my $rightLen = $rightM + $rightD;
               my $FlankDesc = "NONE";
               if(($leftLen >= $MIN_FLANK_LEN) && ($rightLen >= $MIN_FLANK_LEN)) { 
                  $FlankDesc = "Dually Flanked";
               }
               else {
                  if($leftLen >= $MIN_FLANK_LEN) { 
                     $FlankDesc = "5\'-Flanked";
                  }
                  else {
                     $FlankDesc = "3\'-Flanked";
                  }
               }
               $currFlankType = $FlankDesc;
               my $leftShear = $wds[3];
               my $rightShear = $leftShear + $leftM + $leftD + $rightM + $rightD - 1;
               $currShearSite = $leftShear . $SEPARATOR2 . $rightShear;
               $currHIVLen = $largestIns;
               $currIntegrationSite = $wds[3] + $leftM + $leftD;
            }
         }
         else {
            my @wds = split(/\t/, $currHIVmatch);
            $currHIVPos = $wds[3];
            my $currCigar = $wds[5]; 
            my $FivePrimeSH = 0;  
            my $ThreePrimeSH = 0;
            my $Flank5Len = 0;
            my $Flank3Len = 0;
            my $FlankDesc = "NONE";
            # print "HIV $currCigar\n";
            if($currCigar =~ /^([0-9]+)[SH]/) { 
               $Flank5Len = $1;
            }
            if($currCigar =~ /([0-9]+)[SH]$/) { 
               $Flank3Len = $1; 
            }
            if(($Flank5Len >= $MIN_FLANK_LEN) && 
               ($Flank3Len >= $MIN_FLANK_LEN)) { 
               my @chrwds = split(/\t/, $currCHRmatch);
               my $chrCigar = $chrwds[5];
               my $leftMatch = -1;
               my $rightMatch = -1;
               if($chrCigar =~ /^([0-9]+)[M]/) {
                  $leftMatch = $1;
               }
               if($chrCigar =~ /([0-9]+)[M]$/) {
                  $rightMatch = $1;
               } 
   
               ## CHECK CHR HIT
               if(($leftMatch > $MIN_FLANK_LEN) && ($rightMatch >= $MIN_FLANK_LEN)) {
                  $FlankDesc = "Dually Flanked";
               }
               else {
                  if($Flank5Len >= $Flank3Len) {
                     $FlankDesc = "5\'-Flanked";
                  }
                  else {
                     $FlankDesc = "3\'-Flanked";
                  }
               }
            }
           else {             
               if($Flank5Len >= $MIN_FLANK_LEN) { 
                  $FlankDesc = "5\'-Flanked";
               }
               if($Flank3Len >= $MIN_FLANK_LEN) { 
                  $FlankDesc = "3\'-Flanked";
               }
            }
            my($numMatch, $numIns, $numDel) = findMatchInsDel($currCigar);
   
            $currHIVLen = $numMatch + $numDel; 
            if(!($currCHRmatch =~ /$SEPARATOR/)) { 
               $currFlankType = $FlankDesc;
            }
         }
         ## Set Everything
         my $cloneNum = getCloneID(\@clonePos, $currIntegrationSite);
         if($cloneNum == -1) {
            my $numClones = @clonePos;
            push(@clonePos, $currIntegrationSite);
            $cloneNum = $numClones + 1;
         }
         my $currCloneNum = $currChrNum .  "_" . $cloneNum; 
         $SHEAR_HASH{$currCCSID} = $currShearSite;
         $INTEGRATION_HASH{$currCCSID} = $currIntegrationSite;
         $CLONE_HASH{$currCCSID} = $currCloneNum;
         $FLANK_HASH{$currCCSID} = $currFlankType;
         $HIV_LEN_HASH{$currCCSID} = $currHIVLen;  
         $HIV_POS_HASH{$currCCSID} = $currHIVPos;
         $POS_HASH{$currCCSID} = $currPosVal;
         $OVERLAP_SEQ_HASH{$currCCSID} = $currSeqOverlap;
      }
      my @ccsK = sort {$CLONE_HASH{$a} cmp $CLONE_HASH{$b}} keys(%CLONE_HASH);
      my $numCCSK = @ccsK;
      for(my $c = 0; $c < $numCCSK; $c++) { 
         #my $currCCSID = $ccsIDs[$c]; 
         my $currCCSID = $ccsK[$c]; 
         my $currFlankDesc = $FLANK_HASH{$currCCSID};
         my $currPos = $POS_HASH{$currCCSID};
         my $currIntegration = $INTEGRATION_HASH{$currCCSID};
         my $currShear = $SHEAR_HASH{$currCCSID};
         my $currClone = $CLONE_HASH{$currCCSID};
         my $currHIVLen = $HIV_LEN_HASH{$currCCSID};
         my $pcrClones = $CLONE_PCR_HASH{$currClone};
         my $currHIVPos = $HIV_POS_HASH{$currCCSID};
         my $currPCRClone = -1;

         if(defined($pcrClones)) { 
            my @wds22 = split(/$SEPARATOR/, $pcrClones);
            my $numWds22 = @wds22;
            for(my $p = 0; $p < $numWds22; $p++) { 
               my($pcrID, $hivLen, $shearPos)  = split(/\|/, $wds22[$p]);
   
               my $hivDIFF = $currHIVLen - $hivLen;
               my @shear1Wds = split(/$SEPARATOR2/, $currShear);
               my @shear2Wds = split(/$SEPARATOR2/, $shearPos);
               my $numShear1 = @shear1Wds;
               my $numShear2 = @shear2Wds;
               if($numShear1 == $numShear2) { ## Both are Dually flanked
                  my $shearOK = 1;
                  for(my $s = 0; $s < $numShear1; $s++) { 
                     my $shearDIFF = $shear1Wds[$s] - $shear2Wds[$s];
                     if($shearDIFF < 0) { $shearDIFF *= -1; }
                     if($shearDIFF > $MAX_SHEAR_DIFF) { 
                        $shearOK = 0;
                     }
                  }
                  if($hivDIFF < 0) { $hivDIFF *= -1; }
                  if(($hivDIFF <= $MAX_HIV_DIFF) && ($shearOK)) { 
                     $currPCRClone = $pcrID;
                  }
               }
            }
            if($currPCRClone eq "-1") {
               my $numPCR = $numWds22 + 1;
               $currPCRClone = $currClone . "_" . $numPCR;
               $CLONE_PCR_HASH{$currClone} .= $SEPARATOR . $currPCRClone . "|" . $currHIVLen . "|" . $currShear;
            } 
         }
         else {
            $currPCRClone = $currClone . "_1";
            $CLONE_PCR_HASH{$currClone} = $currPCRClone . "|" . $currHIVLen . "|" . $currShear;
         }
         $PCR_HASH{$currCCSID} = $currPCRClone;
      }
   }
   printHostHits(\%PCR_HASH, \%CLONE_HASH, \%FLANK_HASH, \%POS_HASH, \%INTEGRATION_HASH, \%SHEAR_HASH, \%HIV_LEN_HASH, \%HIV_POS_HASH, \%OVERLAP_SEQ_HASH, $OUTFILE);
}

processNoHostReads(\%HIV_NOHOST_HASH, $ERRFILE, $WARNFILE, $OUTFILE);
my @noHostCCS = keys(%HIV_NOHOST_HASH);
my $numNoHost = @noHostCCS;
close($OUTFILE);
close($ERRFILE);
close($WARNFILE);

#-----------------------------------------------------------------------
sub processNoHostReads {
   my($hREF, $EFH, $WarnFH, $OFH) = @_;
   my %h = %$hREF;
   my @ccs = keys(%h);
   my $numCCS = @ccs;
   my %flankHASH;
   my %lenHASH;
   my %posHASH;
   my @goodCCS;

   for(my $i = 0; $i < $numCCS; $i++) { 
      my $currCCS = $ccs[$i];
      my $currHIV = $h{$currCCS};
      my @hivARR = split(/$SEPARATOR/, $currHIV);
      my $numHits = @hivARR;
      if($numHits == 1) { 
         push(@goodCCS, $currCCS);
         my @wds = split(/\t/, $hivARR[0]);
         my $currCigar = $wds[5];
         my $currPos = $wds[3];
         my ($numM, $numI, $numD) = findMatchInsDel($currCigar);
         my $currLen = $numM + $numD;
         my $leftMask = 0;
         my $rightMask = 0;
         if($currCigar =~ /^([0-9]+)[HS]/) { 
            $leftMask = $1;
         }
         if($currCigar =~ /([0-9]+)[HS]$/) { 
            $rightMask = $1;
         }
         if($leftMask >= 50) { 
            $flankHASH{$currCCS} = "No Flank Potential Episomal HIV";
         }
         else { 
            $flankHASH{$currCCS} = "No Flank Linear HIV";
         }
         $lenHASH{$currCCS} = $currLen;
         $posHASH{$currCCS} = $currPos;
      }
      else {
         if($numHits > 2) {
            my @posARR;
            my @lenARR;
            for(my $j = 0; $j < $numHits; $j++) { 
               my @tmpWds = split(/\t/, $hivARR[$j]);
               push(@posARR, $tmpWds[3]);
               my ($tmpM, $tmpI, $tmpD) = findMatchInsDel($tmpWds[5]);
               my $tmpL = $tmpM + $tmpD;
               push(@lenARR, $tmpL);
            }
            my $tmpPos = join($SEPARATOR2, @posARR);
            my $tmpLen = join($SEPARATOR2, @lenARR);
            print $EFH "MORE THAN 2 HITS FOR $currCCS (NO HOST -- HIV ONLY) POS: $tmpPos LEN: $tmpLen\n";
         }
         else {
            push(@goodCCS, $currCCS);
            my $hit1 = $hivARR[0];
            my $hit2 = $hivARR[1];
            my @wds1 = split(/\t/, $hivARR[0]);
            my @wds2 = split(/\t/, $hivARR[1]);
            my $cigar1 = $wds1[5];
            my $cigar2 = $wds2[5];
            my $cPos1  = $wds1[3];
            my $cPos2  = $wds2[3];

            my ($numM1, $numI1, $numD1) = findMatchInsDel($cigar1);
            my ($numM2, $numI2, $numD2) = findMatchInsDel($cigar2);
            my $len1 = $numM1 + $numD1;
            my $len2 = $numM2 + $numD2;
            my $leftMask1 = 0;
            my $rightMask1 = 0;
            my $leftMask2 = 0;
            my $rightMask2 = 0;
            if($cigar1 =~ /^([0-9]+)[HS]/) {
               $leftMask1 = $1;
            }
            if($cigar1 =~ /([0-9]+)[HS]$/) { 
               $rightMask1 = $1;
            }
            if($cigar2 =~ /^([0-9]+)[HS]/) {
               $leftMask2 = $1;
            }
            if($cigar2 =~ /([0-9]+)[HS]$/) { 
               $rightMask2 = $1;
            }
            my $currPos = $cPos1 . $SEPARATOR2 . $cPos2;
            my $currLen = $len1 . $SEPARATOR2 . $len2;
            if($cPos1 > $cPos2) { 
               $currPos = $cPos2 . $SEPARATOR2 . $cPos1;
               $currLen = $len2 . $SEPARATOR2 . $len1;
               if(($rightMask1 > 50) && ($leftMask2 > 50)) { 
                  $flankHASH{$currCCS} = "No Flank Potential episomal HIV";
               }
               else {
                  print $WarnFH "$currCCS is HIV only with two hits yet linear\n";
                  $flankHASH{$currCCS} = "No Flank Linear HIV";
               }
              
            } 
            else {
               if(($rightMask2 > 50) && ($leftMask1 > 50)) { 
                  $flankHASH{$currCCS} = "No Flank Potential episomal HIV";
               }
               else {
                  $flankHASH{$currCCS} = "No Flank Linear HIV";
                  print $WarnFH "$currCCS is HIV only with two hits yet linear\n";
               }
            }
            $lenHASH{$currCCS} = $currLen;
            $posHASH{$currCCS} = $currPos;

            
         }
      }
   }
 
   ####################################
   ## ASSIGN PCR CLONES FOR HIV HITS ##
   ####################################
   my $numGood = @goodCCS;
   my @HIV_Clones;
   my %cloneHASH;

   for(my $i = 0; $i < $numGood; $i++) { 
      my $currCCS = $goodCCS[$i];
      my $currLenInfo = $lenHASH{$currCCS};
      my $currPosInfo = $posHASH{$currCCS};
      my @currLenARR = split(/$SEPARATOR2/, $currLenInfo);
      my @currPosARR = split(/$SEPARATOR2/, $currPosInfo);
      my $currNumHits = @currLenARR;

      my $numClones = @HIV_Clones;
      my $clonePos = -1;
      for(my $j = 0; $j < $numClones; $j++) {
         my $currCloneInfo = $HIV_Clones[$j];
         my($cloneLenInfo, $clonePosInfo) = split(/\|/, $currCloneInfo);
         my @cloneLenARR = split(/$SEPARATOR2/, $cloneLenInfo);
         my @clonePosARR = split(/$SEPARATOR2/, $clonePosInfo);
         my $cloneNumHits = @cloneLenARR;

         if($cloneNumHits == $currNumHits) {
            my $maxLenDiff = 0;
            my $maxPosDiff = 0;
            for(my $k = 0; $k < $cloneNumHits; $k++) {
               my $currLenDiff = $currLenARR[$k] - $cloneLenARR[$k];
               my $currPosDiff = $currPosARR[$k] - $clonePosARR[$k];
               if($currLenDiff < 0) { $currLenDiff *= -1; }
               if($currPosDiff < 0) { $currPosDiff *= -1; }
               if($currLenDiff > $maxLenDiff) {
                  $maxLenDiff = $currLenDiff;
               }
               if($currPosDiff > $maxPosDiff) {
                  $maxPosDiff = $currPosDiff;
               }
            }
            if(($maxLenDiff <= $MAX_HIV_DIFF) && ($maxPosDiff <= $MAX_HIV_POS_DIFF)) { 
               $clonePos = $j;
            }   
         }            
      }
      if($clonePos == -1) {
         push(@HIV_Clones, $currLenInfo . "|" . $currPosInfo);
         $clonePos = $numClones;
      }
      $clonePos++;
      $cloneHASH{$currCCS} = "HIV_$clonePos";
   }
   my @k = sort {$cloneHASH{$a} cmp $cloneHASH{$b}} keys(%cloneHASH);
   my $numK = @k;
   for(my $i = 0; $i < $numK; $i++) { 
      my $currCCS = $k[$i];
      my $currFlank = $flankHASH{$currCCS};
      my $currLen = $lenHASH{$currCCS};
      my $currPos = $posHASH{$currCCS};
      my $currClone = $cloneHASH{$currCCS};
      print $OFH "$currCCS\t$currFlank\tNA\tNA\tNA\t$currClone\t$currClone\t$currPos\t$currLen\t\n";
   }

}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub getOverlappingInformation {
   my($ccs, $s1, $s2, $errFH, $warnFH) = @_;

   my $shear1 = -1;
   my $shear2 = -1;
   my $int1 = -1;
   my $int2 = -1;
   my $overlapBases = "";

   my @wds1 = split(/\t/, $s1);
   my @wds2 = split(/\t/, $s2);
   my $cigar1 = $wds1[5];
   my $cigar2 = $wds2[5];
   my $HIVlen1=-1;
   my $HIVlen2=-1; 
   if($cigar1 =~ /^([0-9]+)[HS](([0-9]+[MID])*[0-9]+M)$/) {
     $HIVlen1 = $1;
     my $c1Match = $2;
     if($cigar2 =~ /^([0-9]+M([0-9]+[MID])*)([0-9]+)[HS]$/) {
        my $c2Match = $1;
        $HIVlen2 = $3;
        my ($c1M, $c1I, $c1D) = findMatchInsDel($c1Match);
        my ($c2M, $c2I, $c2D) = findMatchInsDel($c2Match);
        $shear1 = $wds2[3];
        $int1 = $shear1 + $c2M + $c2D;
        $int2 = $wds1[3];
        $shear2 = $int2 + $c1M + $c1D;
        my $nOverlap = $int1 - $int2 - 1;
        if($nOverlap > 0) {
           $overlapBases = substr($wds2[9], length($wds2[9]) - $nOverlap, $nOverlap);
        }
     }
     else {
        print $errFH "UNRESOLVABLE DUPLICATE ENTRY FOR\t$ccs\n";
     }
   }
   else {
      if($cigar2 =~ /^([0-9]+)[HS](([0-9]+[MID])*[0-9]+M)$/) {
         $HIVlen1 = $1;
         my $c2Match = $2;
         if($cigar1 =~ /^([0-9]+M([0-9]+[MID])*)([0-9]+)[HS]$/) {
            $HIVlen2 = $3;
            my $c1Match = $1;
            my ($c1M, $c1I, $c1D) = findMatchInsDel($c1Match);
            my ($c2M, $c2I, $c2D) = findMatchInsDel($c2Match);
            $shear1 = $wds1[3];
            $int1 = $shear1 + $c1M + $c1D;
            $int2 = $wds2[3];
            $shear2 = $int2 + $c2M + $c2D;
            my $nOverlap = $int1 - $int2 - 1;
            if($nOverlap > 0) {
               $overlapBases = substr($wds1[9], length($wds1[9]) - $nOverlap, $nOverlap);
            }
         }
         else {
            print $errFH "UNRESOLVABLE DUPLICATE ENTRY FOR\t$ccs\n";
         }
      }
      else {
        print $errFH "UNRESOLVABLE DUPLICATE ENTRY FOR\t$ccs\n";
      }
   }
   my $HIVlen = $HIVlen1;
   if($HIVlen2 > $HIVlen1) { 
      $HIVlen = $HIVlen2;
   }
   return($shear1, $shear2, $int1, $int2, $overlapBases, $HIVlen);
}


#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub updateHash {
   my($hREF, $lab, $ccs, $l, $errFH, $warnFH) = @_;

   my $SAM_chr = "TEST";
   my %h = %$hREF;
   if(defined($h{$ccs})) { 
      if($lab eq "HIV") { 
         my @wds1 = split(/\t/, $l);
         my @wds2 = split(/\t/, $h{$ccs});
         my $len1 = length($wds1[9]);
         my $len2 = length($wds2[9]);
         if(($len1 > $len2) && ($len2 < 300)) { 
           $h{$ccs} = $l;
         }
         else {
            if(($len1 > 300) && ($len2 > 300)) {
               if(!($SAM_chr eq "NF")) { 
                  print $warnFH "POTENTIAL REARRANGEMENT: Two separate long HIV entries for $ccs: $len1 and $len2\n";
               }
               $h{$ccs} .= $SEPARATOR . $l;
            }
         }            
      }
      else {
         $h{$ccs} .= $SEPARATOR . $l;
      }
   }
   else {
      $h{$ccs} = $l;
   }
   return(%h);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub parseSAMFile {
   my($bamFN, $chREF, $hhREF, $ccshREF, $EFH, $WFH) = @_;

   my %cHash = %$chREF;
   my %hHash = %$hhREF; 
   my %ccsHash = %$ccshREF;

   ##########################################
   ## Remove secondary hits -- only retain ##
   ## Primary and supplementary hits       ##
   ##########################################

   my $cmd = "samtools view -F 256 $bamFN > tmp.sam";
   system($cmd);
   print "parsed sam file created\n";

   open(SAMFILE, "tmp.sam") || die("Error in opening tmp.sam");
   while(defined(my $line = <SAMFILE>)) { 
      chomp($line);
      my @wds = split(/\t/, $line);
      my $numWds = @wds;
      ## KEY SAM FIELDS:
      ## 0: ccs read id
      ## 2: chromosome
      ## 3: chromosome position
      ## 5: cigar string

      my $CCS_ID = $wds[0];
      my $chr = $wds[2];
      if($chr eq "*") { ## Non mapping -- skip
      }
      else {
         $ccsHash{$CCS_ID} = $CCS_ID;

         if($chr eq $HIV_CHR_NAME) { 
            %hHash = updateHash(\%hHash, "HIV", $CCS_ID, $line, $EFH, $WFH);
         }
         else {
            %cHash = updateHash(\%cHash, "--CHR--", $CCS_ID, $line, $EFH, $WFH);
         }
      }
   }
   close(SAMFILE);
   return(\%cHash, \%hHash, \%ccsHash);
}
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
sub findMatchInsDel {
   my($cStr) = @_;

   my $numMatch = 0;
   my $numIns = 0;
   my $numDel = 0;
   while($cStr =~ /([0-9]+)M/g) { 
      $numMatch+=$1;
   } 
   while($cStr =~ /([0-9]+)I/g) { 
      $numIns+=$1;
   } 
   while($cStr =~ /([0-9]+)[DN]/g) { 
      $numDel+=$1;
   }
   return($numMatch, $numIns, $numDel);
}
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
sub getCloneID {
   my($cPosREF, $intSite) = @_;
 
   my @cPos = @$cPosREF; 
   my $cNum = -1;
   my $numClones = @cPos;

   for(my $cl = 0; $cl < $numClones; $cl++) {
      my $currDiff = $cPos[$cl] - $intSite;
      if($currDiff < 0) { 
         $currDiff *= -1;
      }
      if($currDiff <= $MAX_INTEGRATION_DIFF) { 
         $cNum = $cl + 1;
      }
   }
   return($cNum);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub getSingleHostHitInformation {
   my($ccsID, $chrMatch, $errFH, $warnFH) = @_;

   #########################################################
   ## parses the reads that only have a single hit to the ##
   ## host genome.  Typically, this means they are single ##
   ## flanking reads, but they could be dually flanked    ##
   #########################################################

   my @wds = split(/\t/, $chrMatch);
   my $currPos = $wds[3];

   my $currCigar = $wds[5];
   my $currIntegrationSite = -1;
   my $currShearSite = -1;
   my($numMatch, $numIns, $numDel) = findMatchInsDel($currCigar);

   my $currHIVLen = 0;
   my $shear1 = -1;
   my $shear2 = -1;
   my $int1 = -1;
   my $int2 = -1;
   my $len1 = -1;
   my $len2 = -1;
   if($currCigar =~ /^([0-9]+)[SH]/) {   ## Mask at beginning of read
      my $maskedLen = $1;

      if($maskedLen > $MIN_FLANK_LEN) {
         ## LEFT FLANKED
         $int1 = $wds[3];
         $shear1 = $wds[3] + $numMatch + $numDel;
         $len1 = $maskedLen;
      }
   }
   if($currCigar =~ /([0-9]+)[SH]$/) {  ## Mask at end of read
      my $maskedLen = $1;
      if($maskedLen > $MIN_FLANK_LEN) {
         $shear2 = $wds[3];
         $int2 = $wds[3] + $numMatch + $numDel;
         $len2 = $maskedLen;
      }
   }
   if(($len1 > $len2) && ($len1 > 0)) {
      $currShearSite = $shear1;
      $currIntegrationSite = $int1;
      $currHIVLen = $len1;
   }
   else {
      if($len2 > 0) {
         $currShearSite = $shear2;
         $currIntegrationSite = $int2;
         $currHIVLen = $len2;
      }
   }
   return($currShearSite, $currIntegrationSite, $currHIVLen, $currPos);
}
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
sub getHIV_CCSIDs {
   my($bamFN, $HIV_lab) = @_;

   #########################################################
   ## parses the bam file for CCS ids that map to the HIV ##
   ## genome specified by the HIV label between positions ##
   ## 1-1000000 (should only be about 10K bases, but to   ##
   ## expanded to be sure)                                ##
   #########################################################

   if(!(-e $bamFN)) {
      die("BAM file $bamFN does not exist\n");
   }

   my %h;
   my $cmd = "samtools view $bamFN $HIV_lab:1-1000000 > tmp$HIV_lab.sam";
   system($cmd);
   open(INFILE, "tmp$HIV_lab.sam") || die;
   open(OUTFILE, ">CCS_ReadIDs_$HIV_lab.txt") || die;
   while(defined(my $line = <INFILE>)) {
      my($ccsID, @rest) = split(/\t/, $line);
      $h{$ccsID} = $ccsID;
      print OUTFILE "$ccsID\n";
   }
   close(INFILE);
   close(OUTFILE);
   return(%h);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub createSAMFile {
   my ($bamFN) = @_;

   #############################################################
   ## Converts the bam file to a sam file so it can be parsed ##
   ## for specific CCS ids                                    ##
   #############################################################

   my @wds = split(/\//, $bamFN);
   my $samFN = $wds[-1];
   $samFN =~ s/\.bam$/\.sam/;

   if(!(-e $samFN)) {
      my $cmd = "samtools view -h $bamFN > $samFN";
      print "RUNNING $cmd\n";
      system($cmd);
   }
   return($samFN);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub getHIVSAM {
   my($fn, $hREF) = @_;

   ########################################################
   ## Filters SAM file to only retain CCS reads that map ##
   ## to the HIV genome (both HIV and host read hits)    ##
   ########################################################

   my %h = %$hREF;
   my $outFN = $fn;
   $outFN =~ s/\.sam$/\.hivFiltered\.sam/;
   open(INFILE, $fn) || die("Error opening $fn for reading");
   open(OUTFILE, ">$outFN") || die("Error opening $outFN for writing");
   while(defined(my $line = <INFILE>)) {
      if($line =~ /^\@/) {     ## Comment lines
         print OUTFILE $line;
      }
      else {                   ## Check to see if the ccs ID is in the HIV hash
         my @wds = split(/\t/, $line);
         if(defined($h{$wds[0]})) {
            print OUTFILE $line;
         }
      }
   }
   close(INFILE);
   close(OUTFILE);
   return($outFN);
}
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
sub printHostHits {
   my($pcrREF, $cloneREF, $flankREF, $posREF, $intREF, $shearREF, $hivLenREF, $hivPosREF, $overlapSeqREF, $OUTFILE) = @_;

   my %PCR_HASH = %$pcrREF;
   my %CLONE_HASH = %$cloneREF;
   my %FLANK_HASH = %$flankREF;
   my %POS_HASH = %$posREF;
   my %INTEGRATION_HASH = %$intREF;
   my %SHEAR_HASH = %$shearREF;
   my %HIV_LEN_HASH = %$hivLenREF;
   my %HIV_POS_HASH = %$hivPosREF;
   my %OVERLAP_SEQ_HASH = %$overlapSeqREF;
   my @ccsK = sort {$POS_HASH{$a} cmp $POS_HASH{$b}} keys(%POS_HASH);
   my $numCCSK = @ccsK;
   my @validChrARR;
   my %validChrHASH;

   #####################################
   ## START WITH STANDARD CHROMOSOMES ##
   #####################################
   for(my $chr = 1; $chr <= 24; $chr++) {
      my $currChr = "chr$chr";
      if($chr == 23) {
         $currChr = "chrX";
      }
      if($chr == 24) {
         $currChr = "chrY";
      }
      push(@validChrARR, $currChr);
      $validChrHASH{$currChr} = $currChr;
   }

   ##########################################
   ## ADD IN OTHER NONSTANDARD CHROMOSOMES ##
   ##########################################
   for(my $c = 0; $c < $numCCSK; $c++) {
      my $currCCSID = $ccsK[$c];
      my $currClone = $CLONE_HASH{$currCCSID};
      my ($currChr, @rest) = split(/\_/, $currClone);
      if(!defined($validChrHASH{$currChr})) {
         $validChrHASH{$currChr} = $currChr;
         push(@validChrARR, $currChr);
      }
   }
   my $numValidChrs = @validChrARR;
   for(my $chr = 0; $chr < $numValidChrs; $chr++) {
      my $currChr = $validChrARR[$chr];
      my $numCCScurrChr = 0;
      my @toPrintARR;
      my %currCloneHASH;
      my @currCloneARR;
      for(my $c = 0; $c < $numCCSK; $c++) {
         my $currCCSID = $ccsK[$c];
         my $currFlankDesc = $FLANK_HASH{$currCCSID};
         my $currPos = $POS_HASH{$currCCSID};
         my $currIntegration = $INTEGRATION_HASH{$currCCSID};
         my $currShear = $SHEAR_HASH{$currCCSID};
         my $currClone = $CLONE_HASH{$currCCSID};
         my $currHIVLen = $HIV_LEN_HASH{$currCCSID};
         my $currHIVPos = $HIV_POS_HASH{$currCCSID};
         my $currPCR = $PCR_HASH{$currCCSID};
         if(!defined($currHIVLen)) { $currHIVLen = "NA"; }
         my $currSeqOverlap = $OVERLAP_SEQ_HASH{$currCCSID};
         my($currCCSChr, @rest) = split(/\_/, $currClone);
         if($currChr eq $currCCSChr) { 
            $numCCScurrChr++;
            if($currPos != -1) {
               if(!defined($currCloneHASH{$currClone})) { 
                  $currCloneHASH{$currClone} = $currClone;
                  push(@currCloneARR, $currClone);
               }
               push(@toPrintARR, "$currCCSID\t$currFlankDesc\t$currPos\t$currIntegration\t$currShear\t$currClone\t$currPCR\t$currHIVPos\t$currHIVLen\t$currSeqOverlap\n");
            }
         }
      }
      my $numClones = @currCloneARR;
      my %cloneRenameHASH;
      for(my $i = 0; $i < $numClones; $i++) { 
         my $currClone = $currCloneARR[$i];
         my $newCloneNum = $i + 1;
         my $newClone = $currChr . "_" . $newCloneNum;
         $cloneRenameHASH{$currClone} = $newClone;
         print "Renaming $currClone --> $newClone\n";
      }

      my $numToPrint = @toPrintARR;
      for(my $i = 0; $i < $numToPrint; $i++) { 
         my @wds = split(/\t/, $toPrintARR[$i]);
         $wds[5] = $cloneRenameHASH{$wds[5]};
         my @wds2 = split(/\_/, $wds[6]);
         $wds[6] = $wds[5] . "_" . $wds2[-1];
         $toPrintARR[$i] = join("\t", @wds);
         print $OUTFILE $toPrintARR[$i];
      }

      if($numCCScurrChr > 0) { 
         print "$currChr Reads: $numCCScurrChr clones; $numClones\n";
      }
   }
}

#-----------------------------------------------------------------------
