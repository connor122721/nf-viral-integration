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
if(($numArgs != 3) && ($numArgs != 4)) { 
   die("USAGE: parseSAM <bamFN> <HIV_label> <output_prefix> [unmasked_reads.fa]\n");
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

my ($bamFN, $HIV_label, $OUTPUT_PRE, $UNMASKED_FA) = @ARGV;
print "retrieving CCS reads that hit HIV (may take some time)\n";
my %HIV_ccsHASH = getHIV_CCSIDs($bamFN, $HIV_label);

################################################################
## Optional: pull HIV_SEQ from an unmasked FASTA instead of    ##
## the BAM SEQ field.  Lookups are keyed by READ NAME, so the  ##
## FASTA and the BAM do not need to be in the same order.      ##
################################################################
my %READ_SEQ;
my $USE_FASTA = 0;
my $FASTA_MISSING = 0;
my $FASTA_LEN_MISMATCH = 0;
my $FASTA_USED = 0;
my $ORIENT_OK = 0;
my $ORIENT_FLIPPED = 0;
my $ORIENT_CONFLICT = 0;
my $ORIENT_UNVERIFIED = 0;
my $MOSTLY_MASKED = 0;
my $MIN_ORIENT_IDENTITY = 0.90;
if(defined($UNMASKED_FA)) {
   if(!(-e $UNMASKED_FA)) {
      die("Unmasked FASTA $UNMASKED_FA does not exist\n");
   }
   print "loading unmasked read sequences from $UNMASKED_FA (name keyed)\n";
   %READ_SEQ = loadReadFasta($UNMASKED_FA, \%HIV_ccsHASH);
   $USE_FASTA = 1;
   my @loaded = keys(%READ_SEQ);
   my $numLoaded = @loaded;
   my @wanted = keys(%HIV_ccsHASH);
   my $numWanted = @wanted;
   print "  matched $numLoaded of $numWanted HIV-hitting reads by name\n";
   if($numLoaded == 0) {
      print "  WARNING: no read names matched -- check that FASTA headers use the same IDs as the BAM\n";
   }
}

print "converting bam to sam for processing (may take some time -- likely several minutes)\n";
my $samFN = createSAMFile($bamFN);
print "parsing sam file for only CCS reads that hit HIV (may take some time)\n";
my $HIVparsedSAM_FN = getHIVSAM($samFN, \%HIV_ccsHASH);


open(my $ERRFILE, ">error_log_$OUTPUT_PRE.txt") || die;
open(my $WARNFILE, ">warning_log_$OUTPUT_PRE.txt") || die;

my $outFN = "$OUTPUT_PRE" . "_MasterOfMasterFrame.tsv";
open(my $OUTFILE, ">$outFN") || die("Error opening $outFN for writing");
print $OUTFILE "READ\tFLANK_TYPE\tPOS\tINTEGRATION_SITE\tSHEAR_SITE\tCLONE_ID\tPCR_ID\tHIV_POS\tHIV_LEN\tHOST_SEQ_OVERLAP\tHIV_SEQ\n";

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
my %HIV_SEQ_HASH;

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

   my @ccsIDs = sort(keys(%CCS_HASH));
   my $numCCS = @ccsIDs;
   my $currShearSite = -1;
   my $currIntegrationSite = -1;
   my $currFlankType = "NONE";
   my $currHIVLen = 0;
   my $currHIVPos = "UNK";
   my $currPosVal = -1;
   my $currChrNum = "UNK";
   for(my $c = 0; $c < $numCCS; $c++) { 
      ###########################################################
      ## Reset per read state.  Without this a read that skips  ##
      ## a branch inherits the previous read's values, and      ##
      ## because hash key order is randomized per process the   ##
      ## inherited value changes from run to run.               ##
      ###########################################################
      $currShearSite = -1;
      $currIntegrationSite = -1;
      $currFlankType = "NONE";
      $currHIVLen = 0;
      $currHIVPos = "UNK";
      $currPosVal = -1;
      $currChrNum = "UNK";
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

         ################################################################
         ## Retrieve the portion of the CCS read that is HIV derived.  ##
         ## If the read has an explicit HIV alignment, pull the query  ##
         ## bases covered by that alignment.  Otherwise (HIV called    ##
         ## from a large insertion in a single host hit), pull the     ##
         ## inserted bases out of the host alignment instead.          ##
         ################################################################
         if(defined($currHIVmatch)) {
            $HIV_SEQ_HASH{$currCCSID} = getHIVReadSeq($currHIVmatch);
         }
         else {
            if($currCHRmatch =~ /$SEPARATOR/) {
               $HIV_SEQ_HASH{$currCCSID} = "NA";
            }
            else {
               $HIV_SEQ_HASH{$currCCSID} = getLargestInsertionSeq($currCHRmatch);
            }
         }
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
   printHostHits(\%PCR_HASH, \%CLONE_HASH, \%FLANK_HASH, \%POS_HASH, \%INTEGRATION_HASH, \%SHEAR_HASH, \%HIV_LEN_HASH, \%HIV_POS_HASH, \%OVERLAP_SEQ_HASH, \%HIV_SEQ_HASH, $OUTFILE);
}

processNoHostReads(\%HIV_NOHOST_HASH, $ERRFILE, $WARNFILE, $OUTFILE);
my @noHostCCS = keys(%HIV_NOHOST_HASH);
my $numNoHost = @noHostCCS;
close($OUTFILE);
if($USE_FASTA) {
   print "\n--- HIV_SEQ source summary ---\n";
   print "resolved from unmasked FASTA by name : $FASTA_USED\n";
   print "read name not found in FASTA         : $FASTA_MISSING\n";
   print "CIGAR/FASTA length conflict          : $FASTA_LEN_MISMATCH\n";
   print "orientation verified as written      : $ORIENT_OK\n";
   print "orientation flipped to match BAM     : $ORIENT_FLIPPED\n";
   print "orientation could not be verified    : $ORIENT_UNVERIFIED\n";
   print "FASTA/BAM disagree both orientations : $ORIENT_CONFLICT\n";
   print "sequences that are >50% N            : $MOSTLY_MASKED\n";
   if($ORIENT_FLIPPED > 0) {
      print "NOTE: the FASTA was written in ALIGNMENT orientation for some reads.\n";
      print "      parseSAM corrected these, but the masking/unmasking step should\n";
      print "      be fixed to use one coordinate system (get_forward_sequence()).\n";
   }
   if($MOSTLY_MASKED > 0) {
      print "NOTE: some HIV_SEQ values are mostly N. unmask.py emits the RESTORED\n";
      print "      region with everything else set to N -- if the restored region is\n";
      print "      the host flank then the HIV window is masked out. Supply the\n";
      print "      pre-masking read FASTA instead.\n";
   }
   if($ORIENT_CONFLICT > 0) {
      print "NOTE: some FASTA records match the BAM in neither orientation, which\n";
      print "      indicates a read identity desync. Run check_bam_fasta_sync.pl.\n";
      print $WARNFILE "HIV_SEQ orientation conflicts: $ORIENT_CONFLICT\n";
   }
   if(($FASTA_MISSING > 0) || ($FASTA_LEN_MISMATCH > 0)) {
      print "NOTE: unresolved reads fell back to the BAM SEQ field.\n";
      print "      A nonzero length conflict count means the FASTA record and the\n";
      print "      BAM record sharing a name describe different reads -- run\n";
      print "      check_bam_fasta_sync.pl before trusting those rows.\n";
      print $WARNFILE "HIV_SEQ FASTA lookup: $FASTA_MISSING missing, $FASTA_LEN_MISMATCH length conflicts\n";
   }
   print "------------------------------\n";
}
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
   my %seqHASH;
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
         $seqHASH{$currCCS} = getHIVSeqFromLine($hivARR[0]);
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
            my $seq1 = getHIVSeqFromLine($hivARR[0]);
            my $seq2 = getHIVSeqFromLine($hivARR[1]);
            my $currPos = $cPos1 . $SEPARATOR2 . $cPos2;
            my $currLen = $len1 . $SEPARATOR2 . $len2;
            my $currSeq = $seq1 . $SEPARATOR2 . $seq2;
            if($cPos1 > $cPos2) { 
               $currPos = $cPos2 . $SEPARATOR2 . $cPos1;
               $currLen = $len2 . $SEPARATOR2 . $len1;
               $currSeq = $seq2 . $SEPARATOR2 . $seq1;
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
            $seqHASH{$currCCS} = $currSeq;

            
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
      my $currSeq = $seqHASH{$currCCS};
      if(!defined($currSeq)) { $currSeq = "NA"; }
      print $OFH "$currCCS\t$currFlank\tNA\tNA\tNA\t$currClone\t$currClone\t$currPos\t$currLen\t\t$currSeq\n";
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
   my($pcrREF, $cloneREF, $flankREF, $posREF, $intREF, $shearREF, $hivLenREF, $hivPosREF, $overlapSeqREF, $hivSeqREF, $OUTFILE) = @_;

   my %PCR_HASH = %$pcrREF;
   my %CLONE_HASH = %$cloneREF;
   my %FLANK_HASH = %$flankREF;
   my %POS_HASH = %$posREF;
   my %INTEGRATION_HASH = %$intREF;
   my %SHEAR_HASH = %$shearREF;
   my %HIV_LEN_HASH = %$hivLenREF;
   my %HIV_POS_HASH = %$hivPosREF;
   my %OVERLAP_SEQ_HASH = %$overlapSeqREF;
   my %HIV_SEQ_HASH = %$hivSeqREF;
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
         my $currHIVSeq = $HIV_SEQ_HASH{$currCCSID};
         if(!defined($currHIVSeq)) { $currHIVSeq = "NA"; }
         my($currCCSChr, @rest) = split(/\_/, $currClone);
         if($currChr eq $currCCSChr) { 
            $numCCScurrChr++;
            if($currPos != -1) {
               if(!defined($currCloneHASH{$currClone})) { 
                  $currCloneHASH{$currClone} = $currClone;
                  push(@currCloneARR, $currClone);
               }
               push(@toPrintARR, "$currCCSID\t$currFlankDesc\t$currPos\t$currIntegration\t$currShear\t$currClone\t$currPCR\t$currHIVPos\t$currHIVLen\t$currSeqOverlap\t$currHIVSeq\n");
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

#-----------------------------------------------------------------------
sub loadReadFasta {
   my($fn, $keepREF) = @_;

   ############################################################
   ## Loads read sequences into a name keyed hash.  Headers   ##
   ## may be "<qname>" or "<qname>/<0|1>" (the strand         ##
   ## qualified form written by unmask.py); both are indexed  ##
   ## so lookup does not depend on the upstream convention.   ##
   ## Only reads present in the keep hash are retained.       ##
   ############################################################

   my %keep = %$keepREF;
   my %h;
   my $key = "";
   my $seq = "";
   open(my $FA, $fn) || die("Error opening $fn for reading\n");
   while(defined(my $line = <$FA>)) {
      chomp($line);
      $line =~ s/\r$//;
      if($line =~ /^>(\S+)/) {
         storeReadRecord(\%h, \%keep, $key, $seq);
         $key = $1;
         $seq = "";
      }
      else {
         $seq .= $line;
      }
   }
   storeReadRecord(\%h, \%keep, $key, $seq);
   close($FA);
   return(%h);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub storeReadRecord {
   my($hREF, $keepREF, $key, $seq) = @_;

   if($key eq "") { return; }
   my $bare = $key;
   if($key =~ /^(.+)\/([01])$/) { $bare = $1; }
   if(!defined($$keepREF{$bare})) { return; }
   $$hREF{$key} = $seq;
   ## also index the bare name when it is unambiguous
   if(!defined($$hREF{$bare})) {
      $$hREF{$bare} = $seq;
   }
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub lookupReadSeq {
   my($qname, $flag) = @_;

   ############################################################
   ## Resolves a read's FASTA record, preferring the record   ##
   ## qualified with this alignment's own strand.             ##
   ############################################################

   my $rev = (($flag & 0x10) ? "1" : "0");
   my $other = (($rev eq "1") ? "0" : "1");
   foreach my $k ($qname . "/" . $rev, $qname . "/" . $other, $qname) {
      if(defined($READ_SEQ{$k})) { return($READ_SEQ{$k}); }
   }
   return(undef);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub seqIdentity {
   my($a, $b) = @_;

   ############################################################
   ## Fraction of comparable (non-N on both sides) positions  ##
   ## that agree.  Returns (identity, numCompared); identity  ##
   ## is -1 when nothing could be compared.                   ##
   ############################################################

   my $len = length($a);
   if($len != length($b)) { return(-1, 0); }
   my $n = 0;
   my $ok = 0;
   for(my $i = 0; $i < $len; $i++) {
      my $x = uc(substr($a, $i, 1));
      my $y = uc(substr($b, $i, 1));
      if(($x eq "N") || ($y eq "N")) { next; }
      $n++;
      if($x eq $y) { $ok++; }
   }
   if($n == 0) { return(-1, 0); }
   return($ok / $n, $n);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub revComp {
   my($s) = @_;
   $s = reverse($s);
   $s =~ tr/ACGTUacgtuRYKMrykmBVDHbvdh/TGCAAtgcaaYRMKyrmkVBHDvbhd/;
   return($s);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub getQuerySpan {
   my($cigar) = @_;

   ############################################################
   ## Returns (leadingClip, alignedQueryLen, totalQueryLen)   ##
   ## in ALIGNMENT orientation.  Hard clips are counted here  ##
   ## because the full read is being recovered from FASTA,    ##
   ## where the hard clipped bases are still present.         ##
   ############################################################

   if(!defined($cigar) || ($cigar eq "*")) { return(0, 0, 0); }
   my $lead = 0;
   my $aln = 0;
   my $total = 0;
   my $seenAln = 0;
   while($cigar =~ /([0-9]+)([MIDNSHP=X])/g) {
      my $opLen = $1;
      my $op = $2;
      if(($op eq "S") || ($op eq "H")) {
         $total += $opLen;
         if(!$seenAln) { $lead += $opLen; }
      }
      else {
         if(($op eq "M") || ($op eq "I") || ($op eq "=") || ($op eq "X")) {
            $aln += $opLen;
            $total += $opLen;
            $seenAln = 1;
         }
      }
   }
   return($lead, $aln, $total);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub getHIVSeqFromLine {
   my($line) = @_;

   ############################################################
   ## Single entry point for HIV sequence retrieval.  When an ##
   ## unmasked FASTA was supplied the sequence is recovered   ##
   ## by read name, which removes any dependence on BAM and   ##
   ## FASTA record ordering.  Falls back to the BAM SEQ field ##
   ## whenever the name cannot be resolved safely.            ##
   ############################################################

   if(!defined($line) || ($line eq "")) { return("NA"); }
   my @wds = split(/\t/, $line);
   my $numWds = @wds;
   if($numWds < 10) { return("NA"); }

   my $name = $wds[0];
   my $flag = 0 + $wds[1];
   my $cigar = $wds[5];
   my $seq = $wds[9];

   if($USE_FASTA) {
      my $full = lookupReadSeq($name, $flag);
      if(!defined($full)) {
         $FASTA_MISSING++;
      }
      else {
         my($lead, $alnLen, $totalLen) = getQuerySpan($cigar);
         if($totalLen != length($full)) {
            ###########################################################
            ## The CIGAR says this read is a different length than    ##
            ## the FASTA record of the same name.  That is a genuine  ##
            ## identity conflict, not a formatting quirk, so refuse   ##
            ## to emit a sequence that may belong to another read.    ##
            ###########################################################
            $FASTA_LEN_MISMATCH++;
         }
         else {
            if(($alnLen > 0) && (($lead + $alnLen) <= $totalLen)) {
               ###########################################################
               ## Do not assume which orientation the FASTA was written  ##
               ## in.  Slice it both ways and let the alignment decide:  ##
               ##   A = record is in original read orientation           ##
               ##   B = record is already in alignment orientation       ##
               ## Score each against the aligned bases carried by the    ##
               ## BAM record, ignoring masked (N) positions.             ##
               ###########################################################
               my $rcFull = revComp($full);
               my $candA = substr((($flag & 0x10) ? $rcFull : $full), $lead, $alnLen);
               my $candB = substr((($flag & 0x10) ? $full : $rcFull), $lead, $alnLen);

               my $bamWindow = getAlignedQuerySeq($cigar, $seq);
               my $chosen = $candA;
               if(($bamWindow ne "NA") && (length($bamWindow) == $alnLen)) {
                  my($idA, $nA) = seqIdentity($candA, $bamWindow);
                  my($idB, $nB) = seqIdentity($candB, $bamWindow);
                  if(($idA < 0) && ($idB < 0)) {
                     $ORIENT_UNVERIFIED++;
                  }
                  else {
                     if($idB > $idA) {
                        if($idB >= $MIN_ORIENT_IDENTITY) {
                           $chosen = $candB;
                           $ORIENT_FLIPPED++;
                           if($ORIENT_FLIPPED <= 20) {
                              print $WARNFILE "ORIENTATION: $name FASTA record was in alignment orientation (fwd=" . sprintf("%.3f", $idA) . " rev=" . sprintf("%.3f", $idB) . "); flipped\n";
                           }
                        }
                        else {
                           $ORIENT_CONFLICT++;
                        }
                     }
                     else {
                        if($idA >= $MIN_ORIENT_IDENTITY) {
                           $ORIENT_OK++;
                        }
                        else {
                           $ORIENT_CONFLICT++;
                           if($ORIENT_CONFLICT <= 20) {
                              print $WARNFILE "ORIENTATION: $name FASTA and BAM disagree in both orientations (fwd=" . sprintf("%.3f", $idA) . " rev=" . sprintf("%.3f", $idB) . ") -- possible read identity desync\n";
                           }
                        }
                     }
                  }
               }
               else {
                  $ORIENT_UNVERIFIED++;
               }

               my $nCount = ($chosen =~ tr/Nn//);
               if(($nCount * 2) > length($chosen)) { $MOSTLY_MASKED++; }

               $FASTA_USED++;
               return($chosen);
            }
         }
      }
   }
   return(getAlignedQuerySeq($cigar, $seq));
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub getAlignedQuerySeq {
   my($cigar, $seq) = @_;

   ############################################################
   ## Walks a CIGAR string and returns only the query (read)  ##
   ## bases that are part of the alignment.  M/=/X/I consume  ##
   ## the read; S is clipped away; D/N/H/P consume no read    ##
   ## bases.  For hard clipped supplementary alignments SEQ   ##
   ## already contains only the aligned bases, so this is a   ##
   ## no-op in that case.  Reverse strand alignments are      ##
   ## stored reverse complemented in the SAM, so the returned ##
   ## sequence is always in reference (HIV) orientation.      ##
   ############################################################

   if(!defined($cigar) || !defined($seq)) { return("NA"); }
   if(($cigar eq "*") || ($seq eq "*") || ($seq eq "")) { return("NA"); }

   my $seqLen = length($seq);
   my $qPos = 0;
   my $aligned = "";
   while($cigar =~ /([0-9]+)([MIDNSHP=X])/g) {
      my $opLen = $1;
      my $op = $2;
      if(($op eq "M") || ($op eq "=") || ($op eq "X") || ($op eq "I")) {
         if(($qPos + $opLen) <= $seqLen) {
            $aligned .= substr($seq, $qPos, $opLen);
         }
         $qPos += $opLen;
      }
      else {
         if($op eq "S") {
            $qPos += $opLen;
         }
      }
   }
   if($aligned eq "") { return("NA"); }
   return($aligned);
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub getHIVReadSeq {
   my($hivMatch) = @_;

   ############################################################
   ## Returns the HIV derived portion of a CCS read given the ##
   ## stored HIV SAM line(s).  Multiple HIV hits (potential   ##
   ## rearrangements) are returned separated by $SEPARATOR2   ##
   ## to match the formatting of the other multi-hit columns. ##
   ############################################################

   if(!defined($hivMatch) || ($hivMatch eq "")) { return("NA"); }

   my @hits = split(/$SEPARATOR/, $hivMatch);
   my $numHits = @hits;
   my @seqARR;
   for(my $i = 0; $i < $numHits; $i++) {
      my @wds = split(/\t/, $hits[$i]);
      my $numWds = @wds;
      if($numWds >= 10) {
         push(@seqARR, getHIVSeqFromLine($hits[$i]));
      }
   }
   my $numSeq = @seqARR;
   if($numSeq == 0) { return("NA"); }
   return(join($SEPARATOR2, @seqARR));
}
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
sub getLargestInsertionSeq {
   my($chrMatch) = @_;

   ############################################################
   ## For reads with a single host hit but no separate HIV    ##
   ## alignment, the HIV segment shows up as a large I in the ##
   ## host CIGAR.  This pulls those inserted bases out of the ##
   ## host SEQ field so HIV_SEQ is still populated.           ##
   ############################################################

   if(!defined($chrMatch) || ($chrMatch eq "")) { return("NA"); }

   my @wds = split(/\t/, $chrMatch);
   my $numWds = @wds;
   if($numWds < 10) { return("NA"); }

   my $cigar = $wds[5];
   my $seq = $wds[9];
   if(!defined($cigar) || !defined($seq)) { return("NA"); }
   if(($cigar eq "*") || ($seq eq "*") || ($seq eq "")) { return("NA"); }

   my $seqLen = length($seq);
   my $qPos = 0;
   my $bestLen = 0;
   my $bestStart = -1;
   while($cigar =~ /([0-9]+)([MIDNSHP=X])/g) {
      my $opLen = $1;
      my $op = $2;
      if($op eq "I") {
         if($opLen > $bestLen) {
            $bestLen = $opLen;
            $bestStart = $qPos;
         }
         $qPos += $opLen;
      }
      else {
         if(($op eq "M") || ($op eq "=") || ($op eq "X") || ($op eq "S")) {
            $qPos += $opLen;
         }
      }
   }
   if(($bestStart >= 0) && ($bestLen > 0)) {
      if($USE_FASTA) {
         my $full = $READ_SEQ{$wds[0]};
         if(defined($full)) {
            my($lead, $alnLen, $totalLen) = getQuerySpan($cigar);
            if($totalLen == length($full)) {
               ## bestStart was measured in SEQ coordinates, which exclude
               ## hard clipped bases; shift it into full read coordinates
               my $hLead = 0;
               if($cigar =~ /^([0-9]+)H/) { $hLead = $1; }
               my $oriented = $full;
               if((0 + $wds[1]) & 0x10) { $oriented = revComp($oriented); }
               my $st = $bestStart + $hLead;
               if(($st + $bestLen) <= length($oriented)) {
                  $FASTA_USED++;
                  return(substr($oriented, $st, $bestLen));
               }
            }
            else {
               $FASTA_LEN_MISMATCH++;
            }
         }
         else {
            $FASTA_MISSING++;
         }
      }
      if(($bestStart + $bestLen) <= $seqLen) {
         return(substr($seq, $bestStart, $bestLen));
      }
   }
   return("NA");
}
#-----------------------------------------------------------------------