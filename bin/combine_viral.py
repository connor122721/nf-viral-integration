import pysam, sys
from os import listdir
from os.path import isfile, join
import re
import xlsxwriter
import warnings

###################################################################
## PROGRAM: combine_viral.py                                     ##
## AUTHOR:  David Sachs, Ichan School of Medicine at Mount Sinai ##
##          Modified by Eric Rouchka, University of Louisville   ##
##          Copyright, University of Louisville                  ##
## DATE:    4/23/2025                                            ##
## LAST MODIFIIED: 4/23/2025                                     ##
## GITHUB:                                                       ##
## LICENSE:                                                      ##
###################################################################

##################################################################################
## USAGE:                                                                       ##
##    combine_viral.py <PREFIX>                                                 ##
##                                                                              ##
## NEED TO REFORMAT AND CREATE SUBROUTINES                                      ##
##################################################################################

warnings.filterwarnings("error")

directory_prefix = sys.argv[1]
directory = "./output/".join(directory_prefix.split('/')[:-1])
prefix = directory_prefix.split('/')[-1]
files = [join(directory, f) for f in listdir(directory) if isfile(join(directory, f))]
files = [f for f in files if prefix==f.split('/')[-1][:len(prefix)]]
hiv_files = [x for x in files if re.search("hiv.[0-9]*.sam", x)]
flanks_file = [x for x in files if re.search("flanks.sam", x)][0]
human_file = [x for x in files if re.search("human.filtered.sam", x)][0]
max_hiv_file =  max([int(x.split('.')[-2]) for x in hiv_files])

hivs = []

for i in range(max_hiv_file):
    hiv_tmp = hiv_files[0].split('.')
    hiv_tmp[-2] = str(i+1)
    hiv_tmp = ".".join(hiv_tmp)
    hivs.append(pysam.Samfile(hiv_tmp, 'r'))  ## hivs is an array of all hiv sam files


flanks = pysam.Samfile(flanks_file, 'r')      ## only one flank and 
human = pysam.Samfile(human_file, 'r')        ## one human samfile for each sample

done = False

flanknum = 0

try:
   flank_tmp = next(flanks)       ## get first flank
except:
   print ("No HIV inserts found")
   sys.exit()

flanknum += 1
flankname = flank_tmp.qname
done = False
max_count = 0
max_count_name = ""
printing = False
n_gt_20 = 0
count = 0

wb = xlsxwriter.Workbook(prefix+'.xlsx')
wb.use_zip64()
ws = wb.add_worksheet()
txtFile = open(prefix + '.tab', "w")   ## ECR Added
csvFile = open(prefix + '.csv', "w")   ## ECR Added

#DNA: black forward, gray backwards
#Flank: green forwards, cyan backwards
#unmapped: red (snp or insertion)
#deletion: <10bp, underline, >10bp, insert it and strikethrough


## Formats for excel file
black_bold = wb.add_format({'underline': False, 'bold': True, 'color': 'black'})
red_bold = wb.add_format({'underline': False,'bold': True, 'color': 'red'})
green_bold = wb.add_format({'underline': False,'bold': True, 'color': 'green'})

black = wb.add_format({'underline': False,'bold': False, 'color': 'black'})
red = wb.add_format({'underline': False,'bold': False, 'color': 'red'})
green = wb.add_format({'underline': False,'bold': False, 'color': 'green'})

black_bold_underline = wb.add_format({'underline': True,'bold': True, 'color': 'black'})
red_bold_underline = wb.add_format({'underline': True,'bold': True, 'color': 'red'})
green_bold_underline = wb.add_format({'underline': True,'bold': True, 'color': 'green'})

black_underline = wb.add_format({'underline': True,'bold': False, 'color': 'black'})
red_underline = wb.add_format({'underline': True,'bold': False, 'color': 'red'})
green_underline = wb.add_format({'underline': True,'bold': False, 'color': 'green'})

fmts = [""]*200
rtffmts =[""]*200

fmts[0] = red
fmts[1] = black
fmts[2] = green

fmts[10] = red_bold
fmts[11] = black_bold
fmts[12] = green_bold

fmts[100] = red_underline
fmts[101] = black_underline
fmts[102] = green_underline

fmts[110] = red_bold_underline
fmts[111] = black_bold_underline
fmts[112] = green_bold_underline

## Formats for RTF file
rtffmts[0] = "\\cf2 "
rtffmts[1] = "\\cf1 "
rtffmts[2] = "\\cf3 "

rtffmts[10] = "\\cf2 \\b "
rtffmts[11] = "\\cf1 \\b "
rtffmts[12] = "\\cf3 \\b "

rtffmts[100] = "\\cf2 \\ul "
rtffmts[101] = "\\cf1 \\ul "
rtffmts[102] = "\\cf3 \\ul "

rtffmts[110] = "\\cf2 \\b \\ul "
rtffmts[111] = "\\cf1 \\b \\ul "
rtffmts[112] = "\\cf3 \\b \\ul "

text_wrap = wb.add_format({'text_wrap': True})
columns = ["RTF_NUM","HUMAN_GROUP", "INSERT", "INSERT_LEN", "LEFT_FLANK", "RIGHT_FLANK", "HUMAN_CHECK", "HUMAN_ALTS", "READ", "HIV_DIR_ERR", "FLANK_DIR_ERR", "HUMAN_MAP_ERR", "OVERLAP_ERR", "UNMAPPED"] 
for col_num, data in enumerate(columns):
    ws.write(0, col_num, data)


## ECR ADDED
#txtFile.write(columns[0])
#csvFile.write(columns[0])

#for i in range(1,len(columns)):
#   txtFile.write('\t' + columns[i])
#   csvFile.write('\,' + columns[i])

txtFile.write(columns[0])
csvFile.write(columns[0])

for i in range(1,len(columns)):
   txtFile.write('\t' + columns[i])
   csvFile.write(',' + columns[i])
txtFile.write('\n')
csvFile.write('\n')

## End ECR add


## Set the coumn widts for excel
sizes = [(len(x)+2) for x in columns]
sizes[columns.index("INSERT")] = 40
sizes[columns.index("READ")] = 35
sizes[columns.index("LEFT_FLANK")] = 40
sizes[columns.index("RIGHT_FLANK")] = 40
sizes[columns.index("HUMAN_CHECK")] = 40
for i in range(len(sizes)):
	ws.set_column(i,i,sizes[i])
rtfnum = 0
rtfseqs = 0
rtffile = ""
def openrtf():
    global rtffile
    rtffile = open(prefix+"."+str(rtfnum)+".rtf","w")
    rtffile.write("{\\rtf1{\\colortbl;")
    rtffile.write("\\red0\\green0\\blue0;")
    rtffile.write("\\red255\\green0\\blue0;")
    rtffile.write("\\red0\\green255\\blue0;")
    rtffile.write("}")
def closertf():
    global rtffile
    rtffile.write("}")
    rtffile.close()
openrtf()
human_groups = []
hiv_groups = []
last_human_group_num = 0
last_hiv_group_num = 0
while (not done):
    row_data = [""]*(len(columns))
    count += 1
    hiv_reads = []
    allnames = "HIV:"
    for h in hivs:
        hiv_reads.append(next(h))
        allnames = allnames + " " + hiv_reads[-1].qname
    human_read  = next(human)
    
    got_flanks = False
    flank_reads = []
    while ((hiv_reads[0].qname.split('/')[:2]==flankname.split('/')[:2]) and (not done)):
        flank_reads.append(flank_tmp)
        try:
            flank_tmp = next(flanks)
        except:
            done = True
        flanknum += 1
        flankname = flank_tmp.qname
    allnames = allnames + " FLANKS:" 
    for f in flank_reads:
        allnames = allnames + " " + f.qname
    allnames = allnames + " HUMAN:"
    allnames = allnames + " " + human_read.qname
    allnames = allnames.split(" ")
    allnames = [x for x in allnames if x not in ['HIV:', 'FLANKS:', 'HUMAN:']]
    allnames = ["/".join(x.split("/")[0:2]) for x in allnames]
    if not all(x==allnames[0] for x in allnames):
        print("NAMES NOT EQUAL")
        sys.exit()
    if 1:#0:#human_read.qname=="m54016_200320_224145/4653955/ccs" or printing: 
        printing = True
        seqlist = hiv_reads[0].get_aligned_pairs()
        seqstr = []
        seq = hiv_reads[0].get_forward_sequence()
        seqformat = [0]*len(seq)
        flanklocs = [0]*len(seq)
        f_alts = ""
        hiv_dir = hiv_reads[0].is_reverse
        row_data[columns.index("READ")] = hiv_reads[0].qname
        hiv_query_coords = [hiv_reads[0].query_alignment_start, hiv_reads[0].query_alignment_end]
        hiv_ref_coords = [hiv_reads[0].reference_start, hiv_reads[0].reference_end]
        if hiv_dir:
            hiv_query_coords = [len(seqformat)-hiv_query_coords[1], len(seqformat)-hiv_query_coords[0]]
        for h in hiv_reads:
            if not h.is_unmapped:
                hiv_query_coords_tmp = [h.query_alignment_start, h.query_alignment_end]
                if h.is_reverse:
                    hiv_query_coords_tmp = [len(seqformat)-hiv_query_coords_tmp[1], len(seqformat)-hiv_query_coords_tmp[0]]
                    
                if hiv_query_coords_tmp[0] < hiv_query_coords[0]:
                    hiv_query_coords[0] = hiv_query_coords_tmp[0]
                if hiv_query_coords_tmp[1] > hiv_query_coords[1]:
                    hiv_query_coords[1] = hiv_query_coords_tmp[1]
                hiv_ref_coords_tmp = [h.reference_start, h.reference_end]
                if hiv_ref_coords_tmp[0] < hiv_ref_coords[0]:
                    hiv_ref_coords[0] = hiv_ref_coords_tmp[0]
                if hiv_ref_coords_tmp[1] > hiv_ref_coords_tmp[1]:
                    hiv_ref_coords[1] = hiv_ref_coords_tmp[1]
                if hiv_dir != h.is_reverse:
                    row_data[columns.index("HIV_DIR_ERR")] = 1
                pairs = h.get_aligned_pairs()
                for p in pairs:
                    if (p[0]!=None) and (p[1]!=None):
                        if h.is_reverse:
                            seqformat[len(seqformat)-p[0]-1] = 11
                        else:
                            seqformat[p[0]] = 1
        row_data[columns.index("RTF_NUM")] =rtfnum
        row_data[columns.index("INSERT")] = hiv_reads[0].reference_name + ":" + str(hiv_ref_coords[0]) + "-" +str(hiv_ref_coords[1])
        row_data[columns.index("INSERT_LEN")] = hiv_query_coords[1]-hiv_query_coords[0]
        flank_ref = ""
        if len(flank_reads)>0:
            if not flank_reads[0].is_unmapped:
                flank_ref = flank_reads[0].reference_name
                flank_dir = flank_reads[0].is_reverse
                flank_coords = [flank_reads[0].reference_start, flank_reads[0].reference_end]
                for f in flank_reads:
                    f_coords = [f.reference_start, f.reference_end]
                    if f_coords[0]<flank_coords[0]:
                        flank_coords[0] = f_coords[0]
                    if f_coords[1]>flank_coords[1]:
                        flank_coords[1] = f_coords[1]
                        
                    min_dist = min(abs(f_coords[0]-flank_coords[0]), abs(f_coords[0]-flank_coords[1]), abs(f_coords[1]-flank_coords[0]), abs(f_coords[1]-flank_coords[1]))
                    if f.reference_name==flank_ref and min_dist<20000:
                        if f.is_reverse != flank_dir:
                            row_data[columns.index("FLANK_DIR_ERR")] = 1
                        pairs = f.get_aligned_pairs(matches_only=True)
                        for p in pairs:
                           if(p[0] != None) and (p[1] != None):
                              if(p[0] < len(seqformat)):
                                 if f.is_reverse:
                                    if (seqformat[len(seqformat)-p[0]-1])%10 == 1:
                                        # print("OVERLAP!R")
                                        # print(p)
                                        # print(count)   
                                        # print(seqformat)
                                        row_data[columns.index("OVERLAP_ERR")] = 1
                                    seqformat[len(seqformat)-p[0]-1] = 12
                                    flanklocs[len(seqformat)-p[0]-1] = p[1]
                                 else:
                                    if (seqformat[p[0]-1])%10 == 1:
                                        # print("OVERLAP!")
                                        # print(p)
                                        # print(count)
                                        # print(seqformat)
                                        row_data[columns.index("OVERLAP_ERR")] = 1
                                    seqformat[p[0]-1] = 2
                                    flanklocs[p[0]-1] = p[1]
                    else:
                        f_alt = f.reference_name+":"+str(f_coords[0])+"-"+str(f_coords[1])
                        if f_alts!="":
                            f_alts = f_alts + ";"
                        f_alts = f_alts + f_alt
        human_group_num = 0
        if not human_read.is_unmapped:
            row_data[columns.index("HUMAN_CHECK")] = human_read.reference_name + ":" + str(human_read.reference_start) + "-" + str(human_read.reference_end)
            human_group_num = last_human_group_num + 1
            grp = [human_group_num, [human_read.reference_name, human_read.reference_start, human_read.reference_end]] 
            for g in human_groups:
                min_dist = min(abs(human_read.reference_start-g[1][1]), abs(human_read.reference_start-g[1][2]), abs(human_read.reference_end-g[1][1]), abs(human_read.reference_end-g[1][2]))
                if (g[1][0]==grp[1][0]) and (min_dist<10000):
                    human_group_num = g[0]
            if human_group_num == last_human_group_num + 1:
                human_groups.append(grp)
                last_human_group_num = human_group_num

            pairs = human_read.get_aligned_pairs()
            for p in pairs:
                if (p[0]!=None) and (p[1]!=None):
                    if human_read.is_reverse:
                        seqformat[len(seqformat)-p[0]-1] += 100
                    else:
                        if p[0]>len(seqformat):
                            print("UH OH")
                            print(seqformat)
                            print(p[0])
                            print(len(seqformat))
                            print(len(human_read.seq))
                            print(len(hiv_reads[0].seq))
                            print(human_read.qname)
                            print(hiv_reads[0].qname)
                        #print(seqformat[p[0]])
                        seqformat[p[0]] += 100
        row_data[columns.index("HUMAN_GROUP")] = human_group_num
        row_data[columns.index("HUMAN_ALTS")] = f_alts
        fmt = seqformat[0]
        seqstr = [fmts[fmt]]
        previ = 0
        got_flank_1 = False
        got_flank_2 = False
        got_hiv = False
        min_hiv = -1
        max_hiv = 0
        unmapped = 0
        flank_1_start = 0
        flank_1_end = 0
        flank_2_start = len(seq)
        flank_2_end = len(seq)
        for i in range(len(seqformat)):
            if (seqformat[i]%10==0):
                unmapped += 1
            if (seqformat[i]%10==1):
                max_hiv = i
                if min_hiv==-1:
                    min_hiv = i
            if (seqformat[i]%10==2) and (not got_hiv):
                if not got_flank_1:
                    flank_1_start = i
                flank_1_end = i
                got_flank_1 = True

            if (seqformat[i]%10==1):
                got_hiv = True
            if (seqformat[i]%10==2) and (got_hiv):
                if not got_flank_2:
                    flank_2_start = i
                got_flank_2 = True
                flank_2_end = i
            if seqformat[i]!=fmt:
                fmt = seqformat[i]
                seqstr.append(seq[previ:i])
                seqstr.append(fmts[fmt])
                previ = i
            if i==(len(seqformat)-1):
                seqstr.append(seq[previ:])
        if got_flank_1:
            flankstr = flank_ref + ":" + str(flanklocs[flank_1_start])+"-"+str(flanklocs[flank_1_end])
            row_data[columns.index("LEFT_FLANK")] = flankstr
            if flank_1_end < flank_1_start:
                print("flank 1 order")
                sys.exit()
        if got_flank_2:
            flankstr = flank_ref + ":" + str(flanklocs[flank_2_start])+"-"+str(flanklocs[flank_2_end])
            row_data[columns.index("RIGHT_FLANK")] = flankstr
            if flank_2_end < flank_2_start:
                print("flank 2 order")
                sys.exit()
        if got_flank_1 and got_flank_2:
            print("Both flanks! " + str(count))
            print(hiv_reads[0].qname)
        seqtest = ""
        for i in seqstr:
            if isinstance(i, str):
                seqtest = seqtest + i
        if seqtest!=seq:
            print("Final test: " + str(seqtest==seq))
            print(seqtest)
            print()
            print(seq)
            sys.exit()
        rtffile.write("\\b0 \\ul0 \\cf1 ")
        rtffile.write(">"+hiv_reads[0].qname+"\\line\n")
        fmt = 0
        for i in range(len(seqstr)):
            if not type(seqstr[i]) is str:
                curfmt = fmts.index(seqstr[i])
                rtffile.write(rtffmts[curfmt])
            else:
                rtffile.write(seqstr[i])
                rtffile.write("\\b0 \\ul0 \\cf1 ")
        rtffile.write("\\b0 \\ul0 \\cf1 ")
        rtffile.write("\\line\n")
        rtfseqs = rtfseqs + 1
        if rtfseqs>25000:
            print("Too many rtfseqs")
            rtfseqs = 0
            rtfnum = rtfnum + 1
            closertf()
            openrtf()
    n_count = 0
    row_data[columns.index("UNMAPPED")] = unmapped
    ## ECR ADDED ##
    if len(flank_reads)>0:
       flank_seq = flank_reads[0].get_forward_sequence()
    human_seq = human_read.get_forward_sequence()
    if human_read.is_unmapped:
        human_seq = 'N'*len(human_seq)
    human_coords = [human_read.qstart, human_read.qend]
    if human_read.is_reverse:
        human_coords = [(len(human_seq)-human_read.qend), len(human_seq)-human_read.qstart]
    
    if not human_read.is_unmapped:
        ## ECR ADDED ##
        n_count = 0
        if len(flank_reads)>0:
           n_count = flank_seq[human_coords[0]:human_coords[1]].count('N')
        if n_count > 20:
            n_gt_20 += 1
        row_data[columns.index("HUMAN_MAP_ERR")]= n_count
    for col_num, data in enumerate(row_data):
        ws.write(count, col_num, data)
        writeStringTXT = ""
        writeStringCSV = ""
        if col_num > 0:
            writeStringTXT = '\t'
            writeStringCSV = ','
        writeStringTXT = writeStringTXT + str(row_data[col_num])
        writeStringCSV = writeStringCSV + str(row_data[col_num])
        txtFile.write(writeStringTXT)
        csvFile.write(writeStringCSV)
    txtFile.write('\n')
    csvFile.write('\n')

closertf()
wb.close()
txtFile.close()
csvFile.close()
sys.exit()






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