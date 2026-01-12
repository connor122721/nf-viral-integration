import xlsxwriter

#===== BEGIN OF OUTPUTfile class =====
class OUTPUTfile:
   fh = ""
   fn = ""

#---------------------------------------------------------------------
   def write(self, s):
      self.fh.write(s)
#---------------------------------------------------------------------

#---------------------------------------------------------------------
   def close(self):
      self.fh.close()
#---------------------------------------------------------------------

#---------------------------------------------------------------------
   def __init__(self, pre, ext, delim, col):
      self.fn = pre + ext
      self.fh = open(self.fn, "w")

      self.fh.write(col[0])
      for i in range(1,len(col)):
         self.fh.write(delim + col[i])
      self.fh.write('\n')
#---------------------------------------------------------------------

#=====  END OF OUTPUTfile class =====
