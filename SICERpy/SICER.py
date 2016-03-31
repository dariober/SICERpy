#!/usr/bin/env python

import argparse
import subprocess
import tempfile
import os
import shutil
import sys
import atexit
import errno

parser = argparse.ArgumentParser(description= """
DESCRIPTION

Run the SICER pipeline using a control and an input file.
    
SEE ALSO:

https://github.com/dariober/SICERpy
    
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))


parser.add_argument('--treatment', '-t',
                   required= False,
                   default= '-',
                   help='''Treatment (pull-down) file in bam format
                   ''')


parser.add_argument('--control', '-c',
                   required= True,
                   help='''Control (input) file in bam format
                   ''')

parser.add_argument('--species', '-s',
                   required= True,
                   help='''Species to use. See or edit lib/GenomeData.py for available species. 
                   ''')

parser.add_argument('--effGenomeSize', '-gs',
                   required= False,
                   default= 0.74,     
                   type= float,              
                   help='''Effective Genome as fraction of the genome size. It depends on read length. Default 0.74.
                   ''')

parser.add_argument('--requiredFlag', '-f',
                   required= False,
                   default= 0,
                   type= int,
                   help='''Keep reads with these bits set in flag. Same as `samtools view -f`. Default 0
                   ''')

parser.add_argument('--filterFlag', '-F',
                   required= False,
                   default= 4,
                   type= int,
                   help='''Discard reads with these bits set in flag. Same as `samtools view -F`. Default 4. 
You probably want to discard also second-in-pair reads, secondary and supplementary alignments, reads failing QC.
                   ''')

parser.add_argument('--mapq', '-q',
                   required= False,
                   default= 5,
                   type= int,
                   help='''Discard reads with mapping quality lower than this. Default 5.
                   ''')

parser.add_argument('--redThresh', '-rt',
                   required= False,
                   default= 0,
                   type= int,
                   help='''Redundancy threshold to keep reads mapping to the same position on the same strand. Default 0 (do not filter for redundancy). 
                   ''')

parser.add_argument('--windowSize', '-w',
                   required= False,
                   default= 200,
	               type= int,
                   help='''Size of the windows to scan the genome. WINDOW_SIZE is the smallest possible island. Default 200.
                   ''')

parser.add_argument('--gapSize', '-g',
                   required= False,
                   default= 3,
    	           type= int,
                   help='''Multiple of window size used to determine the gap size. Must be an integer. Default: 3.
                   ''')

parser.add_argument('--fragSize', '-fs',
                   required= False,
                   default= 150,
	           type= int,
                   help='''Size of the sequenced fragment. The center of the the fragment will be taken as half the fragment size. Default 150.
                   ''')

parser.add_argument('--keeptmp',
                   action= 'store_true',
                   help='''For debugging: Do not delete temp directory at the end of run.
                   ''')

parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')


args= parser.parse_args()

# ------------------------------------------------------------------------------

## Get dir where working scripts are.
srcDir= os.path.join(os.path.dirname(os.path.realpath(__file__)), 'src')
python= sys.executable ## Path to python itself
pythonpath= os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lib')


## Set tmp dir
tmpdir= tempfile.mkdtemp(prefix= 'tmp_sicer_', dir= os.getcwd())
if not args.keeptmp:
    atexit.register(shutil.rmtree, tmpdir)

## Create output dir see also 
## http://stackoverflow.com/questions/273192/how-to-check-if-a-directory-exists-and-create-it-if-necessary
#try:
#    os.makedirs(args.outputDir)
#except OSError as exception:
#    if exception.errno != errno.EEXIST:
#        raise
## Work with full paths
# outputDir= os.path.abspath(args.outputDir)

treatment= os.path.abspath(args.treatment)
control= os.path.abspath(args.control)

## Basename
#if not args.basename:
#    basename= os.path.basename(treatment).split('.')[0]
#else:
#    basename= args.basename

## Remove reduntant reads
## ======================
#if args.redThresh > 0:
sys.stderr.write("\n*** Preprocess raw files to remove reduntant reads\n")

filteredSampleBam= os.path.join(tmpdir, os.path.basename(treatment).split('.')[0] + '.removed.bam')
filteredControlBam= os.path.join(tmpdir, os.path.basename(control).split('.')[0] + '.removed.bam')

procs= []
for inBam, outBam in zip([treatment, control], [filteredSampleBam, filteredControlBam]):
    # Create a separate dir for each process so tmp files don't bother each other
    tmpRedDir= os.path.join(tmpdir, 'tmp_' + os.path.basename(inBam) + '_dir')
    os.makedirs(tmpRedDir)
    cmd= """cd %(tmpRedDir)s
PYTHONPATH=%(pythonpath)s
%(python)s %(script)s -s %(species)s -t %(redThresh)s -b %(inBam)s -o %(outBam)s -f %(requiredFlag)s -F %(filterFlag)s -q %(mapq)s""" \
        %{'tmpRedDir': tmpRedDir,
          'pythonpath': pythonpath, 
          'python': python, 
          'script': os.path.join(srcDir, 'remove_redundant_reads_bam.py'), 
          'species': args.species, 
          'redThresh': args.redThresh,
          'inBam': inBam, 
          'outBam': outBam,
          'requiredFlag': args.requiredFlag,
          'filterFlag': args.filterFlag,
          'mapq': args.mapq};
    sys.stderr.write(cmd + '\n\n')
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    procs.append(p)

for p in procs:
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        sys.stderr.write(stderr + '\n')
        sys.exit(p.returncode)
#else:
#    filteredSampleBam= treatment
#    filteredControlBam= control

## Partion the genome in windows
## =============================
sys.stderr.write('\n*** Partion the genome in windows\n')
summaryGraph= os.path.join(tmpdir, 'summary.bedgraph')
cmd= """cd %(tmpdir)s
PYTHONPATH=%(pythonpath)s
%(python)s %(script)s -s %(species)s -b %(filteredSampleBam)s -w %(windowSize)s -i %(fragSize)s -o %(summaryGraph)s""" \
            %{'tmpdir': tmpdir,
              'pythonpath': pythonpath, 
              'python': python, 
              'script': os.path.join(srcDir, 'run-make-graph-file-by-chrom_bam.py'), 
              'species': args.species, 
              'filteredSampleBam': filteredSampleBam,
              'windowSize': args.windowSize,
              'fragSize': args.fragSize,
              'summaryGraph': summaryGraph};
sys.stderr.write(cmd + '\n')

p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
stdout, stderr= p.communicate()
if p.returncode != 0:
    sys.stderr.write(stderr + '\n')
    sys.exit(p.returncode)

## Find candidate islands exhibiting clustering
## ============================================
sys.stderr.write('\n*** Find candidate islands exhibiting clustering\n')
island= os.path.join(tmpdir, 'scoreisland.bed')
cmd= """PYTHONPATH=%(pythonpath)s
%(python)s %(script)s -s %(species)s -b %(summaryGraph)s -w %(windowSize)s -g %(gapSize)s -t %(effGenomeSize)s -e %(evalue)s  -f %(island)s""" \
            %{'pythonpath': pythonpath, 
              'python': python, 
              'script': os.path.join(srcDir, 'find_islands_in_pr.py'), 
              'species': args.species,
              'summaryGraph': summaryGraph, 
              'windowSize': args.windowSize,
              'gapSize': args.gapSize * args.windowSize,
              'effGenomeSize': args.effGenomeSize,
              'evalue': 1000,
              'island': island};
sys.stderr.write(cmd + '\n')
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
stdout, stderr= p.communicate()
if p.returncode != 0:
    sys.stderr.write(stderr + '\n')
    sys.exit(p.returncode)

## Calculate significance of candidate islands using the control library
## =====================================================================
sys.stderr.write('\n*** Calculate significance of candidate islands using the control library\n')
islandSig= os.path.join(tmpdir, 'island-summary.bed')
cmd= """PYTHONPATH=%(pythonpath)s
%(python)s %(script)s -s %(species)s  -a %(filteredSampleBam)s -b %(filteredControlBam)s -d %(island)s -f %(fragSize)s -t %(effGenomeSize)s -o %(islandSig)s""" \
            %{'pythonpath': pythonpath, 
              'python': python, 
              'script': os.path.join(srcDir, 'associate_tags_with_chip_and_control_w_fc_q_bam.py'), 
              'species': args.species,
              'filteredSampleBam': filteredSampleBam, 
              'filteredControlBam': filteredControlBam,
              'island': island,
              'fragSize': args.fragSize,
              'effGenomeSize': args.effGenomeSize,
              'islandSig': islandSig};

sys.stderr.write(cmd + '\n')
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
stdout, stderr= p.communicate()
if p.returncode != 0:
    sys.stderr.write(stderr + '\n')
    sys.exit(p.returncode)

## Finally print to stdout
fin= open(islandSig)
for line in fin:
    sys.stdout.write(line)
fin.close()

sys.exit()



