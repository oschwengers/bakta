#!/usr/bin/env python3

#
# author: M Radey (email: marad_at_uw.edu)
#

import os, sys, glob, argparse
import multiprocessing as mp
from subprocess import run
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from tqdm import tqdm
from shutil import which
from colorama import init as cinit
from colorama import Fore, Back, Style

version     = "1.0.0"

cinit(autoreset=True)

def do_args():
    desc = "Make user region GFFs to use with the Bakta --regions feature"
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__), description=desc)
    parser.add_argument("seqs", help="specifies the path to a Fasta file with the CDS " +\
        "sequences to be searched for in the genomes.")
    parser.add_argument("genomes", help="specifies the path to a directory with the Fasta " +\
        "sequences of the genomes to be processed.")
    parser.add_argument("outdir", help="specifies the output directory")
    parser.add_argument("-t", "--threads", type=int, default=mp.cpu_count(),\
        help="specifies the number of threads to use. The default is %(default)s.")
    parser.add_argument("-m", "--minident", type=float, default=0.9,\
        help="specifies the minimum identity value to use for matching, where sequence " +\
        "matches below the threshold will be excluded. The default is %(default)s.")
    return parser.parse_args()

def search_seqs_task(vals):
    db, genome, outdir, minident, diamond = vals
    outfile = str(outdir) + "/" + Path(genome).stem + ".tab"
    # str(minident), '--outfmt', '"6 qseqid qstart qend qstrand pident"', '--out', outfile]
    args1 = [diamond, 'blastx', '--quiet', '--threads', '1', '--db', db, '-q', genome, '--id',\
        str(minident), '--outfmt', '6', 'qseqid', 'sseqid', 'qstart', 'qend', 'qstrand', 'pident',\
        'qseq', '--out', outfile]
    run(args1, shell=False)

def search_seqs(seqs, genomes, outdir, minident, threads):
    diamond = which("diamond")
    if diamond == None:
        print(Fore.RED + "ERROR: Could not find diamond program.")
        sys.exit()
    if not outdir.exists(): os.mkdir(str(outdir))
    mydb = str(outdir) + "/db"
    print(Fore.CYAN + "Making sequence database...")
    args1 = [diamond, 'makedb', '--quiet', '--in', str(seqs), '-d', mydb]
    run(args1, shell=False)
    infiles = glob.glob(str(genomes) + "/*.fa")
    infiles += glob.glob(str(genomes) + "/*.fna")
    infiles += glob.glob(str(genomes) + "/*.fsa")
    infiles += glob.glob(str(genomes) + "/*.fasta")
    print(Fore.CYAN + "Searching genomes...")
    with tqdm(total=len(infiles), colour='green') as pbar:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = []
            for infile in infiles:
                futures.append(executor.submit(search_seqs_task, (mydb, infile, str(outdir),\
                    minident, diamond)))
            for future in as_completed(futures):
                if future.cancelled():
                    print(Fore.RED + "ERROR: Future was cancelled.")
                elif not future.exception() == None:
                    print(Fore.RED + "ERROR: Future returned exception:")
                    raise future.exception()
                pbar.update(1)
    os.remove(mydb + ".dmnd")

def process_results(outdir):
    print(Fore.CYAN + "Processing results...")
    infiles = glob.glob(str(outdir) + "/*.tab")
    for infile in tqdm(infiles, colour='green'):
        outfile = str(Path(infile).parent) + "/" + Path(infile).stem + ".gff3"
        with open(outfile, 'wt') as o:
            o.write('##gff-version 3\n')
            names = []
            seqs = []
            with open(infile, 'rt') as n:
                for line in n:
                    qseqid, sseqid, qstart, qend, qstrand, pident,\
                        qseq = line.rstrip('\n').split('\t')
                    o.write('\t'.join([qseqid, 'user', 'CDS', qstart, qend, '.', qstrand,\
                        '0', 'gene=' + sseqid]) + '\n') 
                    names.append(sseqid)
                    seqs.append(qseq)
            o.write('##FASTA\n')
            for name, seq in zip(names, seqs):
                o.write('>' + name + '\n' + seq + '\n')
        os.remove(infile)

def main():
    # setup
    args = do_args()
    print(Fore.GREEN + Style.BRIGHT + "\u2022 Using " + str(args.threads) + " processor cores.")
    args.seqs = Path(args.seqs).absolute()
    args.genomes = Path(args.genomes).absolute()
    args.outdir = Path(args.outdir).absolute()
    search_seqs(args.seqs, args.genomes, args.outdir, args.minident, args.threads)
    process_results(args.outdir)
    print(Fore.CYAN + "Done.")
    return 0

if __name__ == "__main__":
   sys.exit(main())

