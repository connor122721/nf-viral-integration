import pysam, sys
import re
from os import listdir
from os.path import isfile, join
from RTFclass import RTFoutput
from OUTPUTclass import OUTPUTfile

class IntegrationData:
   hivs = []
   flanks = ""
   human = ""
   flanknum = 0
   flankname = ""
   done = False
   max_count = 0
   max_count_name = ""
   printing = False
   n_gt_20 = 0
   count = 0
   flank_tmp = ""
   hiv_reads = []
   human_read = ""
   allnames = ""
   got_flanks = False
   flank_reads = []
   seqstr = []
   seqformat = []
   flanklocs = []
   seq = ""
   f_alts = ""
   flank_ref = ""
   hiv_ref_coords = ""
   hiv_query_coords = ""
   row_data = ""
   human_group_num = 0
   last_human_group_num = 0
   last_hiv_group_num = 0
   human_groups = []
   files = ""
   hiv_files = ""
   flanks_file = ""
   human_file = ""
   max_hiv_file = ""
   unmapped = 0

   #-------------------------------------------------------------------------------------------
   def setHostInformation(self):
      # NOTHING
      test = True
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def setOrientationFormat(self, RTF):

      ################################################################################  
      # This function determines the orientation of the HIV integration and          #
      # resets the HIV integration and flanking host integration locations if needed #
      # and then sets the RTF format string to either black (fwd) or black bold (rev)#
      ################################################################################

      self.hiv_query_coords = [self.hiv_reads[0].query_alignment_start, self.hiv_reads[0].query_alignment_end]
      self.hiv_ref_coords   = [self.hiv_reads[0].reference_start, self.hiv_reads[0].reference_end]   
 
      hiv_orient = self.hiv_reads[0].is_reverse
   
      if hiv_orient:  ## Reverse orientation
         self.hiv_query_coords = [len(self.seqformat)-self.hiv_query_coords[1], len(self.seqformat)-self.hiv_query_coords[0]]
         for h in self.hiv_reads:
            if not h.is_unmapped:
               hiv_query_coords_tmp = [h.query_alignment_start, h.query_alignment_end]
               if h.is_reverse:
                  # Reassign HIV coordinates
                  hiv_query_coords_tmp = [len(self.seqformat)-hiv_query_coords_tmp[1], len(self.seqformat)-hiv_query_coords_tmp[0]]              
               if hiv_query_coords_tmp[0] < self.hiv_query_coords[0]:
                  self.hiv_query_coords[0] = hiv_query_coords_tmp[0]
               if hiv_query_coords_tmp[1] > self.hiv_query_coords[1]:
                  self.hiv_query_coords[1] = hiv_query_coords_tmp[1]
               
               # Reassign reference coordinates
               hiv_ref_coords_tmp = [h.reference_start, h.reference_end]
               if hiv_ref_coords_tmp[0] < self.hiv_ref_coords[0]:
                  self.hiv_ref_coords[0] = hiv_ref_coords_tmp[0]
               if hiv_ref_coords_tmp[1] > self.hiv_ref_coords[1]:
                  self.hiv_ref_coords[1] = hiv_ref_coords_tmp[1]
               if hiv_orient != h.is_reverse:
                  self.row_data[RTF.columns.index("HIV_DIR_ERR")] = 1
               pairs = h.get_aligned_pairs()
               for p in pairs:
                   if (p[0]!=None) and (p[1]!=None):
                      if h.is_reverse:
                         self.seqformat[len(self.seqformat)-p[0]-1] = 11  ## black bold -- reverse
                      else:
                         self.seqformat[p[0]] = 1                    ## black  -- forward
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def setIntegrationCols(self, RTF):

      ##################################################################
      # This function sets the integration location and length columns #
      ##################################################################
   
      self.row_data[RTF.columns.index("RTF_NUM")] = RTF.rtfnum
      self.row_data[RTF.columns.index("INSERT")] = self.hiv_reads[0].reference_name + ":" + str(self.hiv_ref_coords[0]) + "-" +str(self.hiv_ref_coords[1])
      self.row_data[RTF.columns.index("INSERT_LEN")] = self.hiv_query_coords[1]-self.hiv_query_coords[0]
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def setViralInformation(self, RTF):

      #######################################################
      ## Sets the information for the integration sequence ##
      #######################################################

      self.seqstr = []
      self.seq = self.hiv_reads[0].get_forward_sequence()
      self.seqformat = [0]*len(self.seq)
      self.flanklocs = [0]*len(self.seq)
      self.f_alts = ""
        
      self.row_data[RTF.columns.index("READ")] = self.hiv_reads[0].qname
      self.setOrientationFormat(RTF)
      self.setIntegrationCols(RTF)

      self.flank_ref = ""
      if len(self.flank_reads)>0:
         if not self.flank_reads[0].is_unmapped:
            self.flank_ref = self.flank_reads[0].reference_name
            self.flank_dir = self.flank_reads[0].is_reverse
            self.flank_coords = [self.flank_reads[0].reference_start, self.flank_reads[0].reference_end]
            for f in self.flank_reads:
               f_coords = [f.reference_start, f.reference_end]
               if f_coords[0]<self.flank_coords[0]:
                  self.flank_coords[0] = f_coords[0]
               if f_coords[1]>self.flank_coords[1]:
                  self.flank_coords[1] = f_coords[1]
                        
               min_dist = min(abs(f_coords[0]-self.flank_coords[0]), abs(f_coords[0]-self.flank_coords[1]), abs(f_coords[1]-self.flank_coords[0]), abs(f_coords[1]-self.flank_coords[1]))
               if f.reference_name==self.flank_ref and min_dist<20000:
                  if f.is_reverse != self.flank_dir:
                     self.row_data[RTF.columns.index("FLANK_DIR_ERR")] = 1
                  pairs = f.get_aligned_pairs(matches_only=True)
                  for p in pairs:
                      if(p[0] != None) and (p[1] != None):
                         if(p[0] < len(self.seqformat)):
                            if f.is_reverse:
                               if (self.seqformat[len(self.seqformat)-p[0]-1])%10 == 1:
                                  self.row_data[RTF.columns.index("OVERLAP_ERR")] = 1
                               self.seqformat[len(self.seqformat)-p[0]-1] = 12
                               self.flanklocs[len(self.seqformat)-p[0]-1] = p[1]
                            else:
                               if (self.seqformat[p[0]-1])%10 == 1:
                                  self.row_data[RTF.columns.index("OVERLAP_ERR")] = 1
                               self.seqformat[p[0]-1] = 2
                               self.flanklocs[p[0]-1] = p[1]
               else:
                  f_alt = f.reference_name+":"+str(f_coords[0])+"-"+str(f_coords[1])
                  if self.f_alts!="":
                     self.f_alts = self.f_alts + ";"
                  self.f_alts = self.f_alts + f_alt
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def setHostInformation(self, RTF):

      ### BEGIN setHostInformation ###
      self.human_group_num = 0
      if not self.human_read.is_unmapped:
         self.row_data[RTF.columns.index("HUMAN_CHECK")] = self.human_read.reference_name + ":" + str(self.human_read.reference_start) + "-" + str(self.human_read.reference_end)
         self.human_group_num = self.last_human_group_num + 1
         grp = [self.human_group_num, [self.human_read.reference_name, self.human_read.reference_start, self.human_read.reference_end]] 
         for g in self.human_groups:
            min_dist = min(abs(self.human_read.reference_start-g[1][1]), abs(self.human_read.reference_start-g[1][2]), abs(self.human_read.reference_end-g[1][1]), abs(self.human_read.reference_end-g[1][2]))
            if (g[1][0]==grp[1][0]) and (min_dist<10000):
               self.human_group_num = g[0]
         if self.human_group_num == self.last_human_group_num + 1:
            self.human_groups.append(grp)
            self.last_human_group_num = self.human_group_num

         pairs = self.human_read.get_aligned_pairs()
         for p in pairs:
            if (p[0]!=None) and (p[1]!=None):
               if self.human_read.is_reverse:
                  self.seqformat[len(self.seqformat)-p[0]-1] += 100
               else:
                  if p[0]>len(self.seqformat):
                     print("UH OH")
                     print(self.seqformat)
                     print(p[0])
                     print(len(self.seqformat))
                     print(len(self.human_read.seq))
                     print(len(self.hiv_reads[0].seq))
                     print(self.human_read.qname)
                     print(self.hiv_reads[0].qname)
                  #print(seqformat[p[0]])
                  self.seqformat[p[0]] += 100
   ### END SetHostInformation

   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def setHIVReads(self):
      self.hiv_reads = []
      self.allnames = "HIV:"
      for h in self.hivs:
         self.hiv_reads.append(next(h))
         self.allnames = self.allnames + " " + self.hiv_reads[-1].qname
      self.human_read  = next(self.human)
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def setFlanksInformation(self, RTF):
      ## BEGIN setFlanksInformation ##
      self.row_data[RTF.columns.index("HUMAN_GROUP")] = self.human_group_num
      self.row_data[RTF.columns.index("HUMAN_ALTS")] = self.f_alts
      fmt = self.seqformat[0]
      self.seqstr = [RTF.fmts[fmt]]
      previ = 0
      got_flank_1 = False
      got_flank_2 = False
      got_hiv = False
      min_hiv = -1
      max_hiv = 0
      self.unmapped = 0
      flank_1_start = 0
      flank_1_end = 0
      flank_2_start = len(self.seq)
      flank_2_end = len(self.seq)
      for i in range(len(self.seqformat)):
         if (self.seqformat[i]%10==0):
            self.unmapped += 1
         if (self.seqformat[i]%10==1):
            max_hiv = i
            if min_hiv==-1:
               min_hiv = i
         if (self.seqformat[i]%10==2) and (not got_hiv):
            if not got_flank_1:
               flank_1_start = i
            flank_1_end = i
            got_flank_1 = True

         if (self.seqformat[i]%10==1):
             got_hiv = True
         if (self.seqformat[i]%10==2) and (got_hiv):
            if not got_flank_2:
               flank_2_start = i
            got_flank_2 = True
            flank_2_end = i
         if self.seqformat[i]!=fmt:
            fmt = self.seqformat[i]
            self.seqstr.append(self.seq[previ:i])
            self.seqstr.append(RTF.fmts[fmt])
            previ = i
         if i==(len(self.seqformat)-1):
            self.seqstr.append(self.seq[previ:])

      ## Now test all fnks to make sure they are in the correct order ##
      ## Flank1, then Flank2                                          ##
      if got_flank_1:
         flankstr = self.flank_ref + ":" + str(self.flanklocs[flank_1_start])+"-"+str(self.flanklocs[flank_1_end])
         self.row_data[RTF.columns.index("LEFT_FLANK")] = flankstr
         if flank_1_end < flank_1_start:
            print("flank 1 order")
            sys.exit()
      if got_flank_2:
         flankstr = self.flank_ref + ":" + str(self.flanklocs[flank_2_start])+"-"+str(self.flanklocs[flank_2_end])
         self.row_data[RTF.columns.index("RIGHT_FLANK")] = flankstr
         if flank_2_end < flank_2_start:
            print("flank 2 order")
            sys.exit()
      if got_flank_1 and got_flank_2:
         print("Both flanks! " + str(self.count))
         print(self.hiv_reads[0].qname)
      seqtest = ""
      for i in self.seqstr:
         if isinstance(i, str):
            seqtest = seqtest + i
      if seqtest!=self.seq:
         print("Final test: " + str(seqtest==seq))
         print(seqtest)
         print()
         print(seq)
         sys.exit()
           
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def writeRTF(self, RTF):
      ## Writes the results to the RTF file ##
      RTF.rtffile.write("\\b0 \\ul0 \\cf1 ")
      RTF.rtffile.write(">"+self.hiv_reads[0].qname+"\\line\n")
      fmt = 0
      for i in range(len(self.seqstr)):
         if not type(self.seqstr[i]) is str:
            curfmt = RTF.fmts.index(self.seqstr[i])
            RTF.rtffile.write(RTF.rtffmts[curfmt])
         else:
            RTF.rtffile.write(self.seqstr[i])
            RTF.rtffile.write("\\b0 \\ul0 \\cf1 ")
      RTF.rtffile.write("\\b0 \\ul0 \\cf1 ")
      RTF.rtffile.write("\\line\n")
      RTF.rtfseqs = RTF.rtfseqs + 1
      if RTF.rtfseqs>25000:
         print("More than 25000")
         #print("Too many rtfseqs")
         #RTF.rtfseqs = 0
         #RTF.rtfnum = RTF.rtfnum + 1
         #RTF.closertf()
         #RTF.openrtf()
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def setFlankNames(self):
      self.got_flanks = False
      self.flank_reads = []

      while ((self.hiv_reads[0].qname.split('/')[:2]==self.flankname.split('/')[:2]) and (not self.done)):
         self.flank_reads.append(self.flank_tmp)
         try:
            self.flank_tmp = next(self.flanks)
         except:
            self.done = True
         self.flanknum += 1
         self.flankname = self.flank_tmp.qname
      self.allnames = self.allnames + " FLANKS:" 
      for f in self.flank_reads:
         self.allnames = self.allnames + " " + f.qname
      self.allnames = self.allnames + " HUMAN:"
      self.allnames = self.allnames + " " + self.human_read.qname
      self.allnames = self.allnames.split(" ")
      self.allnames = [x for x in self.allnames if x not in ['HIV:', 'FLANKS:', 'HUMAN:']]
      self.allnames = ["/".join(x.split("/")[0:2]) for x in self.allnames]
      if not all(x==self.allnames[0] for x in self.allnames):
         print("NAMES NOT EQUAL")
         sys.exit() 
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def checkHostMappings(self, RTF):

      n_count = 0
      self.row_data[RTF.columns.index("UNMAPPED")] = self.unmapped
      ## ECR ADDED ##
      if len(self.flank_reads)>0:
         flank_seq = self.flank_reads[0].get_forward_sequence()

      ## If host read does not map, substitute N's ##
      human_seq = self.human_read.get_forward_sequence()
      if self.human_read.is_unmapped:
         human_seq = 'N'*len(human_seq)

      ## Reverse coordinates if necessary ##
      human_coords = [self.human_read.qstart, self.human_read.qend]
      if self.human_read.is_reverse:
         human_coords = [(len(human_seq)-self.human_read.qend), len(human_seq)-self.human_read.qstart]
   
      ## Host read maps -- make sure there are fewer than 20 N's ## 
      if not self.human_read.is_unmapped:
         ## ECR ADDED ##
         n_count = 0
         if len(self.flank_reads)>0:
            n_count = flank_seq[human_coords[0]:human_coords[1]].count('N')
         if n_count > 20:
            self.n_gt_20 += 1
         self.row_data[RTF.columns.index("HUMAN_MAP_ERR")]= n_count
  
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def getHIVArr(self):

      for i in range(self.max_hiv_file):
          hiv_tmp = self.hiv_files[0].split('.')
          hiv_tmp[-2] = str(i+1)
          hiv_tmp = ".".join(hiv_tmp)
          print("HIV TMP "+hiv_tmp)
          self.hivs.append(pysam.Samfile(hiv_tmp, 'r'))  ## hivs is an array of all hiv sam files
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def writeFiles(self, RTF, txtFile, csvFile):
      for col_num, data in enumerate(self.row_data):
         RTF.ws.write(self.count, col_num, data)
         writeStringTXT = ""
         writeStringCSV = ""
         if col_num > 0:
            writeStringTXT = '\t'
            writeStringCSV = ','
         writeStringTXT = writeStringTXT + str(self.row_data[col_num])
         writeStringCSV = writeStringCSV + str(self.row_data[col_num])
         txtFile.write(writeStringTXT)
         csvFile.write(writeStringCSV)
      txtFile.write('\n')
      csvFile.write('\n') 
   #-------------------------------------------------------------------------------------------

   #-------------------------------------------------------------------------------------------
   def testZMW(self):

      ###################################################################
      ## THIS FUNCTION IS CODE THAT DOES NOT GET RUN BUT IS ADDED HERE ##
      ## AS A STUB FOR FUTURE CHECKS                                   ##
      ###################################################################
      left_read = ""
      right_read = ""
      prev_hiv_zmw = ""
      human_read = ""

      left_read = next(left)
      right_read = next(right)
      human_read = next(human)
      print("HIV_LEN,LEFT_CHR,LEFT_POS,HIV_REF,RIGHT_CHR,RIGHT_POS,LEFT_SEQ,HIV_SEQ,RIGHT_SEQ,LEFT_READ,HIV_READ,RIGHT_READ")
      for hiv_read in hiv:
          hiv_zmw = hiv_read.qname.split('/')[1]

          if hiv_zmw!=prev_hiv_zmw:
             got_left = False
             left_read_final = ""
             right_read_final = ""
             got_right = False
             got_human = False
             human_read_final = ""
             left_zmw = left_read.qname.split('/')[1]
             right_zmw = right_read.qname.split('/')[1]
             human_zmw = human_read.qname.split('/')[1]
             while left_zmw == hiv_zmw:
                if '/'.join(left_read.qname.split('/')[:-2]) == hiv_read.qname:
                   got_left = True
                   left_read_final = left_read
                try:
                   left_read = next(left)
                except:
                   break
                left_zmw = left_read.qname.split('/')[1]
             while right_zmw == hiv_zmw:
                if '/'.join(right_read.qname.split('/')[:-2]) == hiv_read.qname:
                   got_right = True
                   right_read_final = right_read
                try:
                   right_read = next(right)
                except:
                   break
                right_zmw = right_read.qname.split('/')[1]
             while human_zmw == hiv_zmw:
                if '/'.join(human_read.qname.split('/')[:-1]) == hiv_read.qname:
                   got_human = True
                   human_read_final = human_read
                try:
                   human_read = next(human)
                except:
                   break
                human_zmw = human_read.qname.split('/')[1]
             hiv_read_pairs = [pair for pair in hiv_read.aligned_pairs if pair[0]!=None and pair[1]!=None] 
             match_len_hiv = len(hiv_read_pairs) 
             match_hiv = hiv_read_pairs[-1][0]
             hiv_read_pairs = [(pair[0]+hiv_read.qstart, pair[1]) for pair in hiv_read_pairs]
             got_human = False
             if (match_len_hiv>1500 and abs(match_hiv - match_len_hiv)<(match_len_hiv*0.1)):
                human_overlap = False
                if got_human:
                   human_read_pairs = [pair for pair in human_read.aligned_pairs if pair[0]!=None and pair[1]!=None]
                   human_read_pairs = [(pair[0]+human_read.qstart, pair[1]) for pair in human_read_pairs]
                   min_hiv = hiv_read_pairs[0][0]
                   max_hiv = hiv_read_pairs[-1][0]
                   overlap_pairs = [pair for pair in human_read_pairs if pair[0]>=min_hiv and pair[0]<=max_hiv]
                   match_len_human = len(overlap_pairs)
                   match_human = overlap_pairs[-1][0] - overlap_pairs[0][0]
                   if abs(match_len_hiv - match_len_human)<(match_len_hiv*0.2) and abs(match_human - match_len_human)<(match_len_human*0.1):
                      human_overlap = True
                human_overlap = False
                if not human_overlap:
                   outstr = str(len(hiv_read.seq[hiv_read.qstart:hiv_read.qend])) + ","
                   if got_left:
                      outstr += left.getrname(left_read_final.rname)
                      outstr += ","+str(left_read_final.pos)
                   else:
                      outstr += ","
                   outstr += ","+hiv.getrname(hiv_read.rname)
                   if got_right:
                      outstr += ","+right.getrname(right_read_final.rname)
                      outstr += ","+str(right_read_final.pos)
                   else:
                      outstr += ",," 
         
                   outstr += "," + hiv_read.seq[:hiv_read.qstart]
                   outstr += "," + hiv_read.seq[hiv_read.qstart:hiv_read.qend]
                   outstr += "," + hiv_read.seq[hiv_read.qend:]
                   if got_left:
                      outstr += ","+left_read_final.qname
                   else:
                      outstr += ","
                   outstr += ","+hiv_read.qname
                   if got_right:
                      outstr += ","+right_read_final.qname
                   else:
                      outstr += ","
                   print(outstr)
          prev_hiv_zmw = hiv_zmw
   #-------------------------------------------------------------------------------------------
   
   #-------------------------------------------------------------------------------------------
   def __init__(self, d, pre):
      self.files = [join(d, f) for f in listdir(d) if isfile(join(d, f))]
      self.files = [f for f in self.files if pre==f.split('/')[-1][:len(pre)]]

      self.hiv_files    = [x for x in self.files if re.search("trim.[0-9]*.sam", x)]
      self.flanks_file  = [x for x in self.files if re.search("flanks.sam", x)][0]
      self.human_file   = [x for x in self.files if re.search("human.filtered.sam", x)][0]
      self.max_hiv_file =  max([int(x.split('.')[-2]) for x in self.hiv_files])
 
      self.getHIVArr()
      self.flanks = pysam.Samfile(self.flanks_file, 'r')      ## only one flank and 
      self.human = pysam.Samfile(self.human_file, 'r')        ## one human samfile for each sample
      try:
         self.flank_tmp = next(self.flanks)              ## get first flank
      except:
         print("No HIV inserts found")
         sys.exit()
      self.flanknum += 1
      self.flankname = self.flank_tmp.qname      
   #-------------------------------------------------------------------------------------------

