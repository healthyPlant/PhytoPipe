#!/usr/bin/env python
''' This script is based on Broad Institute Viral-ngs assembly.py order_and_orient
    It's used to align and oreint contigs to the reference genome.
'''

__author__ = "xiaojun.hu@usda.gov"

# built-ins
import argparse
import logging
import random
import numpy
import os, sys
import os.path
import shutil
import subprocess
import functools
import operator
import concurrent.futures
import csv
from itertools import zip_longest    # pylint: disable=E0611

#sys.path.append("..")
#export PYTHONPATH="$PYTHONPATH:/path_to_myapp/myapp/myapp/"
import tool_mummer
import util_misc
import util_file

# third-party
import Bio.AlignIO
import Bio.SeqIO
import Bio.Data.IUPACData

log = logging.getLogger(__name__)


class IncompleteAssemblyError(Exception):
    def __init__(self, actual_n, expected_n):
        super(IncompleteAssemblyError, self).__init__(
            'All computed scaffolds are incomplete. Best assembly has {} contigs, expected {}'.format(
                actual_n, expected_n)
        )

class BadInputError(RuntimeError):

    '''Indicates that an invalid input was given to a command'''

    def __init__(self, reason):
        super(BadInputError, self).__init__(reason)

def check_input(condition, error_msg):
    '''Check input to a command'''
    if not condition:
        raise BadInputError(error_msg)

def _order_and_orient_orig(inFasta, inReference, outFasta,
        outAlternateContigs=None,
        breaklen=None, # aligner='nucmer', circular=False, trimmed_contigs=None,
        maxgap=200, minmatch=10, mincluster=None,
        min_pct_id=0.6, min_contig_len=200, min_pct_contig_aligned=0.3):
    ''' This step cleans up the de novo assembly with a known reference genome.
        Uses MUMmer (nucmer or promer) to create a reference-based consensus
        sequence of aligned contigs (with runs of N's in between the de novo
        contigs).
    '''
    mummer = tool_mummer.MummerTool()
    #if trimmed_contigs:
    #    trimmed = trimmed_contigs
    #else:
    #    trimmed = util.file.mkstempfname('.trimmed.contigs.fasta')
    #mummer.trim_contigs(inReference, inFasta, trimmed,
    #        aligner=aligner, circular=circular, extend=False, breaklen=breaklen,
    #        min_pct_id=min_pct_id, min_contig_len=min_contig_len,
    #        min_pct_contig_aligned=min_pct_contig_aligned)
    #mummer.scaffold_contigs(inReference, trimmed, outFasta,
    #        aligner=aligner, circular=circular, extend=True, breaklen=breaklen,
    #        min_pct_id=min_pct_id, min_contig_len=min_contig_len,
    #        min_pct_contig_aligned=min_pct_contig_aligned)
    mummer.scaffold_contigs_custom(
        inReference,
        inFasta,
        outFasta,
        outAlternateContigs=outAlternateContigs,
        extend=True,
        breaklen=breaklen,
        min_pct_id=min_pct_id,
        min_contig_len=min_contig_len,
        maxgap=maxgap,
        minmatch=minmatch,
        mincluster=mincluster,
        min_pct_contig_aligned=min_pct_contig_aligned
    )
    #if not trimmed_contigs:
    #    os.unlink(trimmed)
    return 0

def _call_order_and_orient_orig(inReference, outFasta, outAlternateContigs, **kwargs):
    return _order_and_orient_orig(inReference=inReference, outFasta=outFasta, outAlternateContigs=outAlternateContigs, **kwargs)

def order_and_orient(inFasta, inReference, outFasta,
        outAlternateContigs=None, outReference=None,
        breaklen=None, # aligner='nucmer', circular=False, trimmed_contigs=None,
        maxgap=200, minmatch=10, mincluster=None,
        min_pct_id=0.6, min_contig_len=200, min_pct_contig_aligned=0.3, n_genome_segments=0, 
        outStats=None, threads=None):
    ''' This step cleans up the de novo assembly with a known reference genome.
        Uses MUMmer (nucmer or promer) to create a reference-based consensus
        sequence of aligned contigs (with runs of N's in between the de novo
        contigs).
    '''

    #chk = util.cmd.check_input
    chk = check_input #modify by xiaojun

    ref_segments_all = [tuple(Bio.SeqIO.parse(inRef, 'fasta')) for inRef in util_misc.make_seq(inReference)] #modify by xiaojun
    chk(ref_segments_all, 'No references given')
    chk(len(set(map(len, ref_segments_all))) == 1, 'All references must have the same number of segments')
    if n_genome_segments:
        if len(ref_segments_all) == 1:
            chk((len(ref_segments_all[0]) % n_genome_segments) == 0,
                'Number of reference sequences in must be a multiple of the number of genome segments')
        else:
            chk(len(ref_segments_all[0]) == n_genome_segments,
                'Number of reference segments must match the n_genome_segments parameter')
    else:
        n_genome_segments = len(ref_segments_all[0])
    ref_segments_all = functools.reduce(operator.concat, ref_segments_all, ())

    n_refs = len(ref_segments_all) // n_genome_segments
    log.info('n_genome_segments={} n_refs={}'.format(n_genome_segments, n_refs))
    ref_ids = []

    #change util.file.tempfnames to util_file.tempfnames by xiaojun
    with util_file.tempfnames(suffixes=[ '.{}.ref.fasta'.format(ref_num) for ref_num in range(n_refs)]) as refs_fasta, \
         util_file.tempfnames(suffixes=[ '.{}.scaffold.fasta'.format(ref_num) for ref_num in range(n_refs)]) as scaffolds_fasta, \
         util_file.tempfnames(suffixes=[ '.{}.altcontig.fasta'.format(ref_num) for ref_num in range(n_refs)]) as alt_contigs_fasta:

        for ref_num in range(n_refs):
            this_ref_segs = ref_segments_all[ref_num*n_genome_segments : (ref_num+1)*n_genome_segments]
            ref_ids.append(this_ref_segs[0].id)
            Bio.SeqIO.write(this_ref_segs, refs_fasta[ref_num], 'fasta')

        with concurrent.futures.ProcessPoolExecutor(max_workers=util_misc.sanitize_thread_count(threads)) as executor:  #change util.misc to util_misc by xiaojun
            retvals = executor.map(functools.partial(_call_order_and_orient_orig, inFasta=inFasta,
                breaklen=breaklen, maxgap=maxgap, minmatch=minmatch, mincluster=mincluster, min_pct_id=min_pct_id,
                min_contig_len=min_contig_len, min_pct_contig_aligned=min_pct_contig_aligned),
                refs_fasta, scaffolds_fasta, alt_contigs_fasta)
        for r in retvals:
            # if an exception is raised by _call_order_and_orient_contig, the
            # concurrent.futures.Executor.map function documentations states
            # that the same exception will be raised when retrieving that entry
            # of the retval iterator. This for loop is intended to reveal any
            # CalledProcessErrors from mummer itself.
            assert r==0

        scaffolds = [tuple(Bio.SeqIO.parse(scaffolds_fasta[ref_num], 'fasta')) for ref_num in range(n_refs)]
        base_counts = [sum([len(seg.seq.ungap('N')) for seg in scaffold]) \
            if len(scaffold)==n_genome_segments else 0 for scaffold in scaffolds]
        best_ref_num = numpy.argmax(base_counts)
        if len(scaffolds[best_ref_num]) != n_genome_segments:
            raise IncompleteAssemblyError(len(scaffolds[best_ref_num]), n_genome_segments)
        log.info('base_counts={} best_ref_num={}'.format(base_counts, best_ref_num))
        shutil.copyfile(scaffolds_fasta[best_ref_num], outFasta)
        if outAlternateContigs:
            shutil.copyfile(alt_contigs_fasta[best_ref_num], outAlternateContigs)
        if outReference:
            shutil.copyfile(refs_fasta[best_ref_num], outReference)
        if outStats:
            ref_ranks = (-numpy.array(base_counts)).argsort().argsort()
            with open(outStats, 'w') as stats_f:
                stats_w = csv.DictWriter(stats_f, fieldnames='ref_num ref_name base_count rank'.split(), delimiter='\t')
                stats_w.writeheader()
                for ref_num, (ref_id, base_count, rank) in enumerate(zip(ref_ids, base_counts, ref_ranks)):
                    stats_w.writerow({'ref_num': ref_num, 'ref_name': ref_id, 'base_count': base_count, 'rank': rank})


def parser_order_and_orient(parser=argparse.ArgumentParser()):
    parser.add_argument('-i','--input',dest='inFasta',  help='Input de novo assembly/contigs, FASTA format.') 
    parser.add_argument('-r','--ref',dest='inReference',  nargs='+',
        help=('Reference genome for ordering, orienting, and merging contigs, FASTA format.  Multiple filenames may be listed, each '
              'containing one reference genome. Alternatively, multiple references may be given by specifying a single filename, '
              'and giving the number of reference segments with the nGenomeSegments parameter.  If multiple references are given, '
              'they must all contain the same number of segments listed in the same order.') 
    )
    parser.add_argument('-o','--output',dest='outFasta', 
        help="""Output assembly, FASTA format, with the same number of
                chromosomes as inReference, and in the same order."""
    ) 
    parser.add_argument(
        '--outAlternateContigs',
        help="""Output sequences (FASTA format) from alternative contigs that mapped,
                but were not chosen for the final output.""",
        default=None
    )

    parser.add_argument('--nGenomeSegments', dest='n_genome_segments', type=int, default=0,
                        help="""Number of genome segments.  If 0 (the default), the `inReference` parameter is treated as one genome.
                        If positive, the `inReference` parameter is treated as a list of genomes of nGenomeSegments each.""")

    parser.add_argument('--outReference', help='Output the reference chosen for scaffolding to this file')
    parser.add_argument('--outStats', help='Output stats used in reference selection')
    #parser.add_argument('--aligner',
    #                    help='nucmer (nucleotide) or promer (six-frame translations) [default: %(default)s]',
    #                    choices=['nucmer', 'promer'],
    #                    default='nucmer')
    #parser.add_argument("--circular",
    #                    help="""Allow contigs to wrap around the ends of the chromosome.""",
    #                    default=False,
    #                    action="store_true",
    #                    dest="circular")
    parser.add_argument(
        "--breaklen",
        "-b",
        help="""Amount to extend alignment clusters by (if --extend).
                        nucmer default 200, promer default 60.""",
        type=int,
        default=None,
        dest="breaklen"
    )
    parser.add_argument(
        "--maxgap",
        "-g",
        help="""Maximum gap between two adjacent matches in a cluster.
                        Our default is %(default)s.
                        nucmer default 90, promer default 30. Manual suggests going to 1000.""",
        type=int,
        default=200,
        dest="maxgap"
    )
    parser.add_argument(
        "--minmatch",
        "-l",
        help="""Minimum length of an maximal exact match.
                        Our default is %(default)s.
                        nucmer default 20, promer default 6.""",
        type=int,
        default=10,
        dest="minmatch"
    )
    parser.add_argument(
        "--mincluster",
        "-c",
        help="""Minimum cluster length.
                        nucmer default 65, promer default 20.""",
        type=int,
        default=None,
        dest="mincluster"
    )
    parser.add_argument(
        "--min_pct_id",
        "-p",
        type=float,
        default=0.6,
        help="show-tiling: minimum percent identity for contig alignment (0.0 - 1.0, default: %(default)s)"
    )
    parser.add_argument(
        "--min_contig_len",
        type=int,
        default=200,
        help="show-tiling: reject contigs smaller than this (default: %(default)s)"
    )
    parser.add_argument(
        "--min_pct_contig_aligned",
        "-v",
        type=float,
        default=0.3,
        help="show-tiling: minimum percent of contig length in alignment (0.0 - 1.0, default: %(default)s)"
    )
    #parser.add_argument("--trimmed_contigs",
    #                    default=None,
    #                    help="optional output file for trimmed contigs")

    #util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None))) #comment by xiaojun
    #util.cmd.attach_main(parser, order_and_orient, split_args=True) #comment by xiaojun

    #return parser
    return parser.parse_args()


def main():
    ### Input arguments
    parser = parser_order_and_orient()
    inFasta = parser.inFasta
    inReference = parser.inReference
    outFasta = parser.outFasta
    try:
        order_and_orient(inFasta, inReference, outFasta)
    except:
        print("Refined " + inFasta + " contigs failed")
    """
        outAlternateContigs=None, outReference=None,
        breaklen=None, # aligner='nucmer', circular=False, trimmed_contigs=None,
        maxgap=200, minmatch=10, mincluster=None,
        min_pct_id=0.6, min_contig_len=200, min_pct_contig_aligned=0.3, n_genome_segments=0, 
        outStats=None, threads=None)
    """

if __name__ == '__main__':
	main()
