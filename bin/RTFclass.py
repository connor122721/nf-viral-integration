import xlsxwriter

#---------------------------------------------------------------------
# Define a class for holding the formatting and results for RTF output
class RTFoutput:
   fmts = [""] * 200
   rtffmts = [""] * 200
   black = 0                   # forward DNA
   red = 0                     # unmapped (SNP or insertion)
   green = 0                   # flank forward
   black_bold = 0
   red_bold = 0
   green_bold = 0
   black_underline = 0
   red_underline = 0
   green_underline = 0
   black_bold_underline = 0
   red_bold_underline = 0
   green_bold_underline = 0
   columns = ["RTF_NUM","HUMAN_GROUP", "INSERT", "INSERT_LEN", "LEFT_FLANK", "RIGHT_FLANK", "HUMAN_CHECK", "HUMAN_ALTS", "READ", "HIV_DIR_ERR", "FLANK_DIR_ERR", "HUMAN_MAP_ERR", "OVERLAP_ERR", "UNMAPPED"] 
   rtfnum = 0
   rtfseqs = 0
   rtffile = ""
   
   def setColors(self):
      #DNA: black forward, gray backwards
      #Flank: green forwards, cyan backwards
      #unmapped: red (snp or insertion)
      #deletion: <10bp, underline, >10bp, insert it and strikethrough

      ## Formats for excel file
      black = self.wb.add_format({'underline': False,'bold': False, 'color': 'black'})
      red = self.wb.add_format({'underline': False,'bold': False, 'color': 'red'})
      green = self.wb.add_format({'underline': False,'bold': False, 'color': 'green'})
    
      black_bold = self.wb.add_format({'underline': False, 'bold': True, 'color': 'black'})
      red_bold = self.wb.add_format({'underline': False,'bold': True, 'color': 'red'})
      green_bold = self.wb.add_format({'underline': False,'bold': True, 'color': 'green'})

      black_underline = self.wb.add_format({'underline': True,'bold': False, 'color': 'black'})
      red_underline = self.wb.add_format({'underline': True,'bold': False, 'color': 'red'})
      green_underline = self.wb.add_format({'underline': True,'bold': False, 'color': 'green'})

      black_bold_underline = self.wb.add_format({'underline': True,'bold': True, 'color': 'black'})
      red_bold_underline = self.wb.add_format({'underline': True,'bold': True, 'color': 'red'})
      green_bold_underline = self.wb.add_format({'underline': True,'bold': True, 'color': 'green'})


   def setFormats(self):
      self.fmts[0] = self.red
      self.fmts[1] = self.black
      self.fmts[2] = self.green

      self.fmts[10] = self.red_bold
      self.fmts[11] = self.black_bold
      self.fmts[12] = self.green_bold

      self.fmts[100] = self.red_underline
      self.fmts[101] = self.black_underline
      self.fmts[102] = self.green_underline

      self.fmts[110] = self.red_bold_underline
      self.fmts[111] = self.black_bold_underline
      self.fmts[112] = self.green_bold_underline


   def setRTFFormats(self):
      ## Formats for RTF file
      self.rtffmts[0] = "\\cf2 "
      self.rtffmts[1] = "\\cf1 "
      self.rtffmts[2] = "\\cf3 "

      self.rtffmts[10] = "\\cf2 \\b "
      self.rtffmts[11] = "\\cf1 \\b "
      self.rtffmts[12] = "\\cf3 \\b "

      self.rtffmts[100] = "\\cf2 \\ul "
      self.rtffmts[101] = "\\cf1 \\ul "
      self.rtffmts[102] = "\\cf3 \\ul "

      self.rtffmts[110] = "\\cf2 \\b \\ul "
      self.rtffmts[111] = "\\cf1 \\b \\ul "
      self.rtffmts[112] = "\\cf3 \\b \\ul "


   def setColumnInformation(self):

      ## Set the coumn widths for excel
      self.text_wrap = self.wb.add_format({'text_wrap': True})
      for col_num, data in enumerate(self.columns):
         self.ws.write(0, col_num, data)
      sizes = [(len(x)+2) for x in self.columns]
      sizes[self.columns.index("INSERT")] = 40
      sizes[self.columns.index("READ")] = 35
      sizes[self.columns.index("LEFT_FLANK")] = 40
      sizes[self.columns.index("RIGHT_FLANK")] = 40
      sizes[self.columns.index("HUMAN_CHECK")] = 40
      for i in range(len(sizes)):
         self.ws.set_column(i,i,sizes[i])

   def openrtf(self, pre, num):
      self.rtffile = open(pre+"."+str(num)+".rtf","w")
      self.rtffile.write("{\\rtf1{\\colortbl;")
      self.rtffile.write("\\red0\\green0\\blue0;")
      self.rtffile.write("\\red255\\green0\\blue0;")
      self.rtffile.write("\\red0\\green255\\blue0;")
      self.rtffile.write("}")

   def closertf(self):
      self.rtffile.write("}")
      self.rtffile.close()
      self.wb.close()
  
   def __init__(self, pre):
      self.wb = xlsxwriter.Workbook(pre+'.xlsx')
      self.wb.use_zip64()
      self.ws = self.wb.add_worksheet()
      self.setColors()
      self.setFormats()
      self.setRTFFormats()
      self.setColumnInformation()
 
#---------------------------------------------------------------------
