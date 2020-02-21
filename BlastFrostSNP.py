import os, sys, numpy as np
from copy import copy
import click, csv, subprocess
from collections import OrderedDict, defaultdict

def readFasta(fname) :
    seq = {}
    with open(fname) as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                seq[name] = []
            else :
                seq[name].extend(line.strip().split())
    for n in seq :
        seq[n] = ''.join(seq[n]).upper()
    return seq

@click.command()
@click.option('-q', '--query', help='FASTA files that contain gene sequences')
@click.option('-s', '--sites', help='Sites to detect in the genes')
@click.option('-o', '--output', help='prefix for the intemediate files and outputs')
@click.option('-b', '--blastfrost', help='executable file for BlastFrost', default='/home/zhemin/software/BlastFrost/build/BlastFrost')
@click.option('-g', '--gfa', help='gfa file for the BiFrost output. color file will be parsed automatically')
@click.option('-d', '--dist', default=2, help='Allowed number of mismatches in the target region. Default: 2')
@click.option('-D', '--extdist', default=2, help='Allowed number of mismatches in surrounding region. Default: 2')
@click.option('-k', '--kmer', default=21, type=int, help='Kmer that was used in BiFrost.')
def main(query, sites, output, blastfrost, gfa, dist, extdist, kmer) :
    sequence = readFasta(query)
    sitelist = OrderedDict()
    with open(sites, 'rt') as fin :
        for part in csv.reader(fin, delimiter='\t') :
            n, c, s0, e0 = part[0], part[1], int(part[2]), int(part[3])
            assert c in sequence, 'Contig {0} is not present'.format(c)
            assert s0 < e0, '{0}: start needs to be smaller than end'.format(n)
            assert e0 <= len(sequence[c]), '{0}: end site is larger than the contig size'.format(n)
            sitelist[n] = [c, s0, e0]
    nodebug = 1
    if nodebug == 1 :
        with open(output + '.kmers', 'wt') as fout :
            for sitename, (contig, start, end) in sitelist.items() :
                seq = sequence[contig]
                varied_seq = list(seq[start - 1:end])
                varied_seqs = [copy(varied_seq)]
                for vi, vs in enumerate(varied_seq) :
                    ori_base = vs
                    vpool = {'A', 'C', 'G', 'T'} - {ori_base}
                    for alt_s in vpool:
                        varied_seq[vi] = alt_s
                        varied_seqs.append(copy(varied_seq))
                        if dist >= 2 :
                            for vi2, vs2 in enumerate(varied_seq[vi+1:]) :
                                ori_v2 = vs2
                                vpool2 = {'A', 'C', 'G', 'T'} - {ori_v2}
                                for alt_v2 in vpool2 :
                                    varied_seq[vi+1+vi2] = alt_v2
                                    varied_seqs.append(copy(varied_seq))
                                varied_seq[vi + 1 + vi2] = ori_v2
                    varied_seq[vi] = ori_base
                s0 = max(end - kmer, 0)
                e0 = min(start, len(seq) - kmer)
                for s in np.arange(s0, e0) :
                    word = list(seq[s:s+kmer])
                    varied_site = start - s - 1
                    word[varied_site:(end - s)] = []
                    for varied_seq in varied_seqs :
                        # d=0
                        nseq = word[:varied_site] + varied_seq + word[varied_site:]
                        nname = '{0}__{1}_0_'.format(sitename, ''.join(varied_seq), s)
                        fout.write('>{0}\n{1}\n'.format(nname, ''.join(nseq)))

                        if extdist >= 1 :
                            for wi, ws in enumerate(word) :
                                ori_w = ws
                                wpool = {'A', 'C', 'G', 'T'} - {ori_w}
                                for alt_w in wpool :
                                    word[wi] = alt_w
                                    # d = 1
                                    nseq = word[:varied_site] + varied_seq + word[varied_site:]
                                    nname = '{0}__{1}_1_{2}-{3}{4}'.format(sitename, ''.join(varied_seq), s, wi, alt_w)
                                    fout.write('>{0}\n{1}\n'.format(nname, ''.join(nseq)))

                                    if extdist == 2 and wi < varied_site :
                                        for wi2, ws2 in enumerate(word[varied_site:]) :
                                            ori_w2 = ws2
                                            wpool2 = {'A', 'C', 'G', 'T'} - {ori_w2}
                                            for alt_w2 in wpool2:
                                                word[varied_site+wi2] = alt_w2
                                                nseq = word[:varied_site] + varied_seq + word[varied_site:]
                                                nname = '{0}__{1}_2_{2}-{3}{4}-{5}{6}'.format(sitename, ''.join(varied_seq), s, wi, alt_w, wi2, alt_w2)
                                                fout.write('>{0}\n{1}\n'.format(nname, ''.join(nseq)))
                                            word[varied_site+wi2] = ori_w2
                                word[wi] = ori_w
    print('{0} -d 0 -k {1} -t 20 -g {2} -f {3} -o {4} -s 5 -q {5}'.format(blastfrost, kmer, gfa, gfa[:-3] + 'bfg_colors', output, output + '.kmers'))
    subprocess.Popen('{0} -d 0 -k {1} -t 20 -g {2} -f {3} -o {4} -s 5 -q {5}'.format(blastfrost, kmer, gfa, gfa[:-3]+'bfg_colors', output, output+'.kmers').split()).communicate()

    variants = defaultdict(dict)
    with open('{0}_{1}.kmers.search'.format(output, os.path.basename(output)), 'rt') as fin :
        for part in csv.reader(fin, delimiter='\t'):
            sitename, sitetype = part[0].split('__')
            site_var, d = sitetype.split('_', 2)[:2]
            key = (sitename, part[1])
            variants[key][site_var] = variants[key].get(site_var, 0) + dist + 1 - int(d)

    with open('{0}_snv.search'.format(output), 'wt') as fout :
        for (sitename, genome), vars in sorted(variants.items()) :
            fout.write('{0}\t{1}\t{2}\n'.format(sitename, genome, '\t'.join([ '{0}:{1}'.format(*v) for v in sorted(vars.items(), key=lambda v:-v[1]) ])))
if __name__ == '__main__' :
    main()