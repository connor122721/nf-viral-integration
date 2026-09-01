#!/usr/bin/env python

import argparse
import re
import sys
import warnings
from itertools import product
from math import prod
from scipy.stats import fisher_exact


def check_width(regexpstring, culprit):
    try:
        tmp = re.compile("(?<=" + regexpstring + ")")
    except Exception:
        raise ValueError(
            f"{culprit} pattern must correspond to a fixed length expression (example: T|GC is not allowed) Please try again using a {culprit} pattern that has only one possible width."
        )


def check_chars(chars, good_chars, error_message):
    bad_chars = [x for x in set(chars) if x not in good_chars]
    if len(bad_chars):
        raise ValueError(f"{error_message}. Yours contains: {bad_chars}")


def check_input_patterns(
    mutfrom, mutto, primaryupstream, primarydownstream, iupac_dict
):
    # check that all patterns consist of IUPAC characters only
    check_chars(
        primaryupstream,
        list(iupac_dict.keys()) + ["|"],
        "The upstream pattern must include only IUPAC characters and '|'.",
    )
    check_chars(
        primarydownstream,
        list(iupac_dict.keys()) + ["|"],
        "The downstream pattern must include only IUPAC characters and '|'.",
    )
    # stop with error if any pattern isn't fixed width
    check_width(primaryupstream, "upstream")
    check_width(primarydownstream, "downstream")
    # require mutation to be only one base
    try:
        mutfrom = iupac_dict[mutfrom]
        mutto = iupac_dict[mutto]
    except Exception:
        raise ValueError(
            "Mutation from and to must each be a single IUPAC character.")
    # check for (undesired) overlap in mutation from and to
    if sum([1 for x in mutfrom if x in mutto]) != 0:
        warnings.warn("Mutation from and to have overlapping bases.")
    # check for (undesired) overlap in context
    primary_pattern = re.sub(
        "\\.",
        "",
        "|".join(
            [
                y + x
                for x in primarydownstream.split("|")
                for y in primaryupstream.split("|")
            ]
        ),
    )
    base_info_primary = [
        [iupac_dict[y] for y in list(x)] for x in primary_pattern.split("|")
    ]
    contexts_primary = [
        "".join(list(y)) for x in base_info_primary for y in product(*x)
    ]
    if len(contexts_primary) != len(set(contexts_primary)):
        raise ValueError(
            "Context is redundant. Please provide non-redundant patterns.")


def check_partial_enforce(match_type, enforce_type):
    if match_type == "partial" and enforce_type != "D":
        raise ValueError("When match is partial, enforce must be D.")


def compute_context_prop(refseq, seq, context, enforce, iupac_dict):
    prop = 1
    if context[0] != "":
        prop = 0
        if enforce == "D":
            for u in context:
                prop += prod(
                    [
                        len([x for x in iupac_dict[s] if x in iupac_dict[c]])
                        / len(iupac_dict[s])
                        for s, c in zip(seq, list(u))
                    ]
                )
        elif enforce == "A":
            for u in context:
                prop += prod(
                    [
                        len([x for x in iupac_dict[s] if x in iupac_dict[c]])
                        / len(iupac_dict[s])
                        for s, c in zip(refseq, list(u))
                    ]
                )
        elif enforce == "B":
            for u in context:
                prop += prod(
                    [
                        len([x for x in iupac_dict[s] if x in iupac_dict[c]])
                        / len(iupac_dict[s])
                        for s, c in zip(seq, list(u))
                    ]
                )
                prop += prod(
                    [
                        len([x for x in iupac_dict[s] if x in iupac_dict[c]])
                        / len(iupac_dict[s])
                        for s, c in zip(refseq, list(u))
                    ]
                )
            prop = prop / 2
    return prop


def slice_seq(sequence, start, context_len, context_type, keep_gaps):
    if context_type == "upstream":
        if keep_gaps:
            seq_sliced = list(sequence[:start][-context_len:])
        else:
            seq_sliced = list(re.sub("-", "", sequence[:start])[-context_len:])
    if context_type == "downstream":
        if keep_gaps:
            seq_sliced = list(sequence[start + 1:][:context_len])
        else:
            seq_sliced = list(
                re.sub("-", "", sequence[start + 1:])[:context_len])
    return seq_sliced


def find_match_weight(
    refseq,
    sequence,
    start,
    end,
    mutto,
    upstream_context,
    downstream_context,
    enforce,
    iupac_dict,
    match,
    keep_gaps,
):
    base = sequence[start:end]
    site_primary = matchval_primary = site_control = matchval_control = 0
    if base != "-":
        upstream_ref = upstream_seq = downstream_ref = downstream_seq = []
        up_same_len = down_same_len = True
        if enforce == "D":
            if upstream_context[0] != "":
                upstream_seq = slice_seq(sequence, start, len(
                    upstream_context[0]), "upstream", keep_gaps)
            if downstream_context[0] != "":
                downstream_seq = slice_seq(sequence, start, len(
                    downstream_context[0]), "downstream", keep_gaps)
                down_same_len = len(
                    downstream_context[0]) == len(downstream_seq)
        elif enforce == "A":
            if upstream_context[0] != "":
                upstream_ref = slice_seq(refseq, start, len(
                    upstream_context[0]), "upstream", keep_gaps)
                up_same_len = len(upstream_context[0]) == len(upstream_ref)
            if downstream_context[0] != "":
                downstream_ref = slice_seq(refseq, start, len(
                    downstream_context[0]), "downstream", keep_gaps)
                down_same_len = len(
                    downstream_context[0]) == len(downstream_ref)
        elif enforce == "B":
            if upstream_context[0] != "":
                upstream_ref = slice_seq(refseq, start, len(
                    upstream_context[0]), "upstream", keep_gaps)
                upstream_seq = slice_seq(sequence, start, len(
                    upstream_context[0]), "upstream", keep_gaps)
                up_same_len = len(
                    upstream_context[0]) == len(upstream_ref) and len(
                    upstream_context[0]) == len(upstream_seq)
            if downstream_context[0] != "":
                downstream_ref = slice_seq(refseq, start, len(
                    downstream_context[0]), "downstream", keep_gaps)
                downstream_seq = slice_seq(sequence, start, len(
                    downstream_context[0]), "downstream", keep_gaps)
                down_same_len = len(downstream_context[0]) == len(
                    downstream_ref
                ) and len(downstream_context[0]) == len(downstream_seq)
        if up_same_len and down_same_len:
            upstream_seq_prop_primary = compute_context_prop(
                upstream_ref, upstream_seq, upstream_context, enforce, iupac_dict)
            downstream_seq_prop_primary = compute_context_prop(
                downstream_ref, downstream_seq, downstream_context, enforce, iupac_dict)
            # 1 means complete primary or control site, fraction means partial
            # primary or control site
            site_primary = upstream_seq_prop_primary * downstream_seq_prop_primary
            site_control = 1 - site_primary
            chars_seq = iupac_dict[base]
            chars_mut = iupac_dict[mutto]
            correct_mut = sum(
                [x in chars_mut for x in chars_seq]) / len(chars_seq)
            matchval_primary = correct_mut * site_primary
            matchval_control = correct_mut * site_control
        if keep_gaps:
            if "-" in upstream_ref + upstream_seq + downstream_ref + downstream_seq:
                site_primary = matchval_primary = site_control = matchval_control = 0
        if match == "strict" and (
            int(site_primary) != site_primary
            or int(matchval_primary) != matchval_primary
            or int(site_control) != site_control
            or int(matchval_control) != matchval_control
        ):  # ignore if multistate strict mode and not complete overlap
            site_primary = matchval_primary = site_control = matchval_control = 0
    return site_primary, matchval_primary, site_control, matchval_control


def summarize_matches(
    refseq,
    queryseq,
    start,
    finish,
    potentialre,
    enforce,
    mutto_orig,
    upstream_orig,
    downstream_orig,
    iupac_dict,
    match,
    keep_gaps,
    seqs,
    name,
    positionsfile,
):
    sites_primary = matches_primary = sites_control = matches_control = 0
    if finish is None:
        potentials = potentialre.finditer(refseq, start)
    else:
        potentials = potentialre.finditer(refseq, start, finish)
    for mymatch in potentials:
        site_primary, matchval_primary, site_control, matchval_control = (
            find_match_weight(
                refseq,
                queryseq,
                mymatch.start(),
                mymatch.end(),
                mutto_orig,
                upstream_orig,
                downstream_orig,
                enforce,
                iupac_dict,
                match,
                keep_gaps,
            )
        )
        sites_primary += site_primary
        matches_primary += matchval_primary
        sites_control += site_control
        matches_control += matchval_control
        if positionsfile is not None:
            if site_primary != 0:
                positionsfile.write(
                    str(seqs)
                    + ","
                    + name
                    + ","
                    + str(mymatch.start() + 1)
                    + ",primary,"
                    + str(site_primary)
                    + ","
                    + str(matchval_primary)
                    + "\n"
                )
            if site_control != 0:
                positionsfile.write(
                    str(seqs)
                    + ","
                    + name
                    + ","
                    + str(mymatch.start() + 1)
                    + ",control,"
                    + str(site_control)
                    + ","
                    + str(matchval_control)
                    + "\n"
                )
    return sites_primary, matches_primary, sites_control, matches_control


def calc_fisher(primarysites, primaries, controlsites, controls):
    return fisher_exact(
        [[primaries, primarysites - primaries],
            [controls, controlsites - controls]],
        alternative="greater",
    )


def check_positive(value):
    try:
        ival = int(value)
    except Exception:
        raise argparse.ArgumentTypeError("must be a positive integer")
    if ival != float(value) or ival < 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return ival


def calc_pval_ratio(primarysites, primaries, controlsites, controls):
    odds_ratio, pval = calc_fisher(
        primarysites, primaries, controlsites, controls)
    try:
        ratio = "%0.2f" % (primaries * controlsites /
                           (1.0 * primarysites * controls))
    except Exception:
        if primaries * controlsites > 0:
            ratio = "inf"
        else:
            ratio = "undef"
    return pval, ratio


def read_seq(fa, chars, error, line=None):
    if line is None:
        line = fa.readline()
        if line[0] != ">":
            raise ValueError("Input alignment must be in FASTA format.")
    name = line[1:].strip()
    seq = ""
    line = fa.readline()
    while line and line[0] != ">":
        seq += line
        line = fa.readline()
    seq = seq.replace("\n", "").upper()
    # check sequence
    check_chars(seq, chars, error)
    return name, seq, line


def parse_args(args, iupac_dict):
    parser = argparse.ArgumentParser(
        prog="Hypermut 3.0",
        description="Identify mutations in a user-defined context")
    parser.add_argument(
        "fasta",
        type=str,
        help="Alignment file in fasta format")
    parser.add_argument(
        "mutationfrom",
        type=str.upper,
        help="Base in the reference to consider as a site of interest for nucleotide substitution",
    )
    parser.add_argument(
        "mutationto",
        type=str.upper,
        help="Base in the query to consider a nucleotide substitution of interest",
    )
    parser.add_argument(
        "--upstreamcontext",
        "-u",
        type=str.upper,
        default="",
        help="Upstream nucleotide context of interest (default: no upstream context)",
    )
    parser.add_argument(
        "--downstreamcontext",
        "-d",
        type=str.upper,
        default="",
        help="Downstream nucleotide context of interest (default: no downstream context)",
    )
    parser.add_argument(
        "--prefix",
        "-p",
        type=str,
        default="",
        help="Prefix for output files (default: no prefix).",
    )
    parser.add_argument(
        "--enforce",
        "-e",
        type=str,
        choices=[
            "A",
            "D",
            "B"],
        default="D",
        help="What sequence to enforce the context on: ancestor/reference (A), descendant/query (D, default), or both (B)",
    )
    parser.add_argument(
        "--match",
        "-m",
        type=str,
        choices=[
            "strict",
            "partial"],
        default="strict",
        help="Whether to include only complete matches (strict, default), or also include partial matches (not completely overlapping bases between query and context, partial)",
    )
    parser.add_argument(
        "--keepgaps",
        "-k",
        action="store_true",
        help="Flag indicating to keep gaps in the alignment when identifying pattern matches (default without flag is to skip gaps)",
    )
    parser.add_argument(
        "--begin",
        "-b",
        type=check_positive,
        default=0,
        help="Position at which to start searching for mutations (default: 0). Note that the context may fall outside of these positions.",
    )
    parser.add_argument(
        "--finish",
        "-f",
        type=check_positive,
        help="Position at which to end searching for mutations (default: end of sequence). Note that the context may fall outside of these positions.",
    )
    args = parser.parse_args(args)

    # only allow partial matches when context is enforced on query sequence
    # only
    check_partial_enforce(args.match, args.enforce)
    # check input patterns
    check_input_patterns(
        args.mutationfrom,
        args.mutationto,
        args.upstreamcontext,
        args.downstreamcontext,
        iupac_dict,
    )
    return args


def loop_through_sequences(fa, args, iupac_dict, summaryfile, positionsfile):
    # reference sequence
    ref_chars = query_chars = list(iupac_dict.keys())
    ref_error = query_error = (
        "Sequences must contain only IUPAC characters or - (for gap)"
    )
    if args.match == "partial":
        ref_chars = list("ACGT-")
        ref_error = "In partial match mode, the reference sequence must contain only the following characters: ACGT-"
    name, refseq, line = read_seq(fa, ref_chars, ref_error)
    seqs = 0
    while line:
        seqs += 1
        name, sequence, line = read_seq(fa, query_chars, query_error, line)
        primarysites, primaries, controlsites, controls = summarize_matches(
            refseq,
            sequence,
            args.begin,
            args.finish,
            re.compile(args.mutationfrom),
            args.enforce,
            args.mutationto,
            args.upstreamcontext.split("|"),
            args.downstreamcontext.split("|"),
            iupac_dict,
            args.match,
            args.keepgaps,
            seqs,
            name,
            positionsfile,
        )

        pval, ratio = calc_pval_ratio(
            primarysites, primaries, controlsites, controls)

        if summaryfile is not None:
            summaryfile.write(
                name
                + ","
                + str(primaries)
                + ","
                + str(primarysites)
                + ","
                + str(controls)
                + ","
                + str(controlsites)
                + ","
                + ratio
                + ",%.6g" % float(pval)
                + "\n"
            )


iupac_dict = {
    "A": list("A"),
    "C": list("C"),
    "G": list("G"),
    "T": list("T"),
    "R": list("AG"),
    "Y": list("CT"),
    "S": list("GC"),
    "W": list("AT"),
    "K": list("GT"),
    "M": list("AC"),
    "B": list("CGT"),
    "D": list("AGT"),
    "H": list("ACT"),
    "V": list("ACG"),
    "N": list("ACGT"),
    "-": list("-"),
}

if __name__ == "__main__":
    args = parse_args(sys.argv[1:], iupac_dict)

    # write args file
    af = open(args.prefix + "args.csv", "w")
    af.write("arg_name,arg_value\n")
    for arg in vars(args):
        af.write(str(arg) + "," + str(getattr(args, arg)) + "\n")
    af.close()

    # prep for writing summary file
    sf = open(args.prefix + "summary.csv", "w")
    sf.write(
        "seq_name,primary_matches,potential_primaries,control_matches,potential_controls,rate_ratio,fisher_p\n"
    )
    # prep for writing positions file
    pf = open(args.prefix + "positions.csv", "w")
    pf.write("seq_num,seq_name,potential_mut_site,context,prop_context,mut_match\n")

    # open fasta file for reading
    fa = open(args.fasta, "r")
    loop_through_sequences(fa, args, iupac_dict, sf, pf)

    sf.close()
    pf.close()
