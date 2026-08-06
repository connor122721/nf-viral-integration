import pysam, sys

###################################################################
## PROGRAM: unmask.py                                            ##
## AUTHOR:  David Sachs, Ichan School of Medicine at Mount Sinai ##
##          Modified by Eric Rouchka, University of Louisville   ##
##          Copyright, University of Louisville                  ##
## DATE:    4/23/2025                                            ##
## LAST MODIFIED: 07/30/2026                                     ##
###################################################################

##################################################################################
## USAGE:                                                                       ##
##    unmask.py <bamfile> <maskfile.fa> [--strict]                              ##
##                                                                              ##
## This program reads in a sam/bam file and adds back in the masked portions    ##
## of the reads, which could either be the inserted viral sequence or the       ##
## flanking host, determined by the presence of an "N" in the alignment that is ##
## a non-N in the actual fasta sequence.                                        ##
##                                                                              ##
## ORIENTATION                                                                  ##
## ------------                                                                 ##
## read.get_forward_sequence() returns the read in ORIGINAL (as sequenced)      ##
## orientation.  A mask record built from read.query_sequence instead is in     ##
## ALIGNMENT orientation, which for a reverse strand read is the reverse        ##
## complement.  Indexing one with the other silently restores the mirrored      ##
## window: the record is the right length, so no length check catches it, and   ##
## only reverse strand reads are affected.                                      ##
##                                                                              ##
## At every non-N position the mask holds the read's real base, so orientation  ##
## can be settled by base identity rather than by convention.  Each read is     ##
## scored in both orientations and the better one is used.  Reads that agree    ##
## in neither orientation are reported and skipped rather than emitted.         ##
##                                                                              ##
## --strict : exit nonzero if any read fails to resolve                         ##
##################################################################################

MIN_AGREEMENT = 0.90   ## fraction of non-N mask positions that must match

args = [a for a in sys.argv[1:] if not a.startswith("--")]
flags = set(a for a in sys.argv[1:] if a.startswith("--"))
if len(args) != 2:
    sys.stderr.write("USAGE: unmask.py <bamfile> <maskfile.fa> [--strict]\n")
    sys.exit(2)
BAMFILE, MASKFILE = args
STRICT = "--strict" in flags

COMP = str.maketrans("ACGTUacgtuRYKMrykmBVDHbvdhNn", "TGCAAtgcaaYRMKyrmkVBHDvbhdNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


#----------------------------------------------------------------------------------
# Load the maskfile. Records may be keyed as "<qname>/<reversed>" or as a bare
# "<qname>". Both are indexed so that lookup never depends on which convention
# the upstream masking step used.
masks = {}
with open(MASKFILE) as mf:
    key = None
    for line in mf:
        line = line.rstrip("\n").rstrip("\r")
        if not line:
            continue
        if line.startswith(">"):
            key = line[1:].split()[0]
            masks[key] = []
        elif key is not None:
            masks[key].append(line)
masks = {k: "".join(v) for k, v in masks.items()}

# secondary index: bare qname -> list of full keys
by_qname = {}
for k in masks:
    bare = k.rsplit("/", 1)[0] if k.rsplit("/", 1)[-1] in ("0", "1") else k
    by_qname.setdefault(bare, []).append(k)
#----------------------------------------------------------------------------------

stats = {
    "emitted": 0,
    "as_written": 0,
    "orientation_corrected": 0,
    "undetermined_all_masked": 0,
    "no_mask_record": 0,
    "length_mismatch": 0,
    "agrees_in_neither_orientation": 0,
    "skipped_not_primary": 0,
    "duplicate_suppressed": 0,
}


def lookup_mask(qname, reversed_flag):
    """Find this read's mask record regardless of the key convention used."""
    for k in (qname + "/" + reversed_flag,
              qname + "/" + ("0" if reversed_flag == "1" else "1"),
              qname):
        if k in masks:
            return masks[k]
    cands = by_qname.get(qname)
    if cands and len(cands) == 1:
        return masks[cands[0]]
    return None


def agreement(mask, seq):
    """Fraction of non-N mask positions whose base matches seq. Returns
    (score, n_compared); n_compared == 0 means the mask carries no evidence."""
    n = 0
    ok = 0
    for i, m in enumerate(mask):
        if m != "N" and m != "n":
            n += 1
            if m.upper() == seq[i].upper():
                ok += 1
    if n == 0:
        return (0.0, 0)
    return (ok / float(n), n)


emitted_keys = set()

#----------------------------------------------------------------------------------
def process_read(read):
    ###############################################################
    ## Function process_read                                     ##
    ## INPUT: read: current sam alignment read                   ##
    ##-----------------------------------------------------------##
    ## Resolves the read's mask record and its orientation, then ##
    ## restores the masked positions and writes a FASTA record   ##
    ## in ORIGINAL read orientation.                             ##
    ###############################################################

    if read.is_secondary or read.is_supplementary:
        stats["skipped_not_primary"] += 1
        return

    reversed_flag = "1" if read.is_reverse else "0"
    qname = read.query_name

    mask = lookup_mask(qname, reversed_flag)
    if mask is None:
        stats["no_mask_record"] += 1
        return

    fullseq = read.get_forward_sequence()
    if fullseq is None or len(fullseq) != len(mask):
        sys.stderr.write("WARN: length mismatch or missing SEQ for %s "
                         "(seq=%s, mask=%d); skipping\n"
                         % (qname,
                            "None" if fullseq is None else len(fullseq),
                            len(mask)))
        stats["length_mismatch"] += 1
        return

    ###########################################################
    ## Settle orientation by base identity instead of trusting ##
    ## the header suffix.  fullseq is original orientation, so ##
    ## the mask must be too; if it scores better reversed, the ##
    ## upstream step wrote it in alignment orientation.        ##
    ###########################################################
    fwd_score, n_cmp = agreement(mask, fullseq)
    if n_cmp == 0:
        # Mask is entirely N: every position gets restored, so the output is
        # just the read itself and orientation cannot change the result.
        stats["undetermined_all_masked"] += 1
        chosen = mask
    else:
        rev_mask = revcomp(mask)
        rev_score, _ = agreement(rev_mask, fullseq)
        if fwd_score >= rev_score and fwd_score >= MIN_AGREEMENT:
            chosen = mask
            stats["as_written"] += 1
        elif rev_score > fwd_score and rev_score >= MIN_AGREEMENT:
            chosen = rev_mask
            stats["orientation_corrected"] += 1
            sys.stderr.write("ORIENT: mask for %s was in alignment orientation "
                             "(fwd=%.3f rev=%.3f, is_reverse=%s); reverse "
                             "complemented before unmasking\n"
                             % (qname, fwd_score, rev_score, read.is_reverse))
        else:
            stats["agrees_in_neither_orientation"] += 1
            sys.stderr.write("MISMATCH: mask for %s matches the read in neither "
                             "orientation (fwd=%.3f rev=%.3f over %d bases); "
                             "skipping -- this record likely belongs to a "
                             "different read\n"
                             % (qname, fwd_score, rev_score, n_cmp))
            return

    out_key = qname + "/" + reversed_flag
    if out_key in emitted_keys:
        stats["duplicate_suppressed"] += 1
        return
    emitted_keys.add(out_key)

    unmasked = []
    for i in range(len(chosen)):
        if chosen[i] == "N" or chosen[i] == "n":
            unmasked.append(fullseq[i])
        else:
            unmasked.append("N")

    stats["emitted"] += 1
    print(">" + out_key)
    print("".join(unmasked))
#----------------------------------------------------------------------------------


mode = "rb" if BAMFILE.lower().endswith(".bam") else "r"
sam = pysam.AlignmentFile(BAMFILE, mode)
for read in sam:
    process_read(read)
sam.close()

sys.stderr.write("\n=========== UNMASK SUMMARY ===========\n")
for k in ("emitted", "as_written", "orientation_corrected",
          "undetermined_all_masked", "duplicate_suppressed",
          "skipped_not_primary", "no_mask_record", "length_mismatch",
          "agrees_in_neither_orientation"):
    sys.stderr.write("%-32s : %d\n" % (k, stats[k]))
sys.stderr.write("======================================\n")

failed = stats["agrees_in_neither_orientation"] + stats["length_mismatch"]
if stats["orientation_corrected"]:
    sys.stderr.write("NOTE: %d read(s) had a mask in alignment orientation. Fix the\n"
                     "      masking step to use get_forward_sequence() so masks and\n"
                     "      reads share one coordinate system.\n"
                     % stats["orientation_corrected"])
if failed and STRICT:
    sys.exit(1)