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
                   required= True,
                   help='''Treatment (pull-down) file in bed format
                   ''')


parser.add_argument('--control', '-c',
                   required= True,
                   help='''Control (input) file in bed format
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

parser.add_argument('--redThresh', '-rt',
                   required= False,
                   default= 1,
                   type= int,
                   help='''Redundancy threshold to keep reads mapping to the same position on the same strand. 
Set to 0 to skip filtering and use input bed as is. Default 1. 
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

parser.add_argument('--fragSize', '-f',
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

#parser.add_argument('--fdr', '-fdr',
#                   required= False,
#                   default= 0.01,
#    	           type= float,
#                   help='''False discovery rate to report results.
#                   ''')
#parser.add_argument('--outputDir', '-o',
#                   required= False,
#                   default= os.getcwd(),
#                   help='''Output directory. Created if does not exist. 
#                   ''')
#parser.add_argument('--basename', '-n',
#                   required= False,
#                   help='''Basename for output files. Default is taken from treatment filename.
#                   ''')

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
if args.redThresh > 0:
    sys.stderr.write("\n*** Preprocess raw files to remove reduntant reads\n")

    filteredSampleBed= os.path.join(tmpdir, os.path.basename(treatment).split('.')[0] + '.removed.bed')
    filteredControlBed= os.path.join(tmpdir, os.path.basename(control).split('.')[0] + '.removed.bed')
    
    procs= []
    for inBed, outBed in zip([treatment, control], [filteredSampleBed, filteredControlBed]):
        # Create a separate dir for each process so tmp files don't bother each other
        tmpRedDir= os.path.join(tmpdir, 'tmp_' + os.path.basename(inBed) + '_dir')
        os.makedirs(tmpRedDir)
        cmd= """cd %(tmpRedDir)s
PYTHONPATH=%(pythonpath)s
%(python)s %(script)s -s %(species)s -t %(redThresh)s -b %(inBed)s -o %(outBed)s""" \
            %{'tmpRedDir': tmpRedDir,
              'pythonpath': pythonpath, 
              'python': python, 
              'script': os.path.join(srcDir, 'remove_redundant_reads.py'), 
              'species': args.species, 
              'redThresh': args.redThresh,
              'inBed': inBed, 
              'outBed': outBed};
        sys.stderr.write(cmd + '\n\n')
        p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
        procs.append(p)
    
    for p in procs:
        stdout, stderr= p.communicate()
        if p.returncode != 0:
            sys.stderr.write(stderr + '\n')
            sys.exit(p.returncode) 
else:
    filteredSampleBed= treatment
    filteredControlBed= control

## Partion the genome in windows
## =============================
sys.stderr.write('\n*** Partion the genome in windows\n')
summaryGraph= os.path.join(tmpdir, 'summary.bedgraph')
cmd= """cd %(tmpdir)s
PYTHONPATH=%(pythonpath)s
%(python)s %(script)s -s %(species)s -b %(filteredSampleBed)s -w %(windowSize)s -i %(fragSize)s -o %(summaryGraph)s""" \
            %{'tmpdir': tmpdir,
              'pythonpath': pythonpath, 
              'python': python, 
              'script': os.path.join(srcDir, 'run-make-graph-file-by-chrom.py'), 
              'species': args.species, 
              'filteredSampleBed': filteredSampleBed,
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
%(python)s %(script)s -s %(species)s  -a %(filteredSampleBed)s -b %(filteredControlBed)s -d %(island)s -f %(fragSize)s -t %(effGenomeSize)s -o %(islandSig)s""" \
            %{'pythonpath': pythonpath, 
              'python': python, 
              'script': os.path.join(srcDir, 'associate_tags_with_chip_and_control_w_fc_q.py'), 
              'species': args.species,
              'filteredSampleBed': filteredSampleBed, 
              'filteredControlBed': filteredControlBed,
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

## POSSIBLY NOT NECESSARY

## Normalize summary graph by total island filtered reads per million
## ==================================================================
#sys.stderr.write('\n*** Normalize summary graph by total island filtered reads per million\n')
#normalizedSummary= os.path.join(outputDir, basename + '.normalized.graph')
#cmd= """%(python)s %(script)s -i %(summaryGraph)s -a 3 -t 1000000 -o %(normalizedSummary)s""" \
#            %{'python': python, 
#              'script': os.path.join(srcDir, 'normalize.py'),
#              'summaryGraph': summaryGraph,
#              'normalizedSummary': normalizedSummary};
#sys.stderr.write(cmd + '\n')
#
#p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
#stdout, stderr= p.communicate()
#if p.returncode != 0:
#    sys.stderr.write(stderr + '\n')
#    sys.exit(p.returncode)
#
## Convert the normalized summary graph into wig vstep format
#sys.stderr.write('\n*** Convert the normalized summary graph into wig vstep format\n')

## Identify significant islands using FDR criterion
## ================================================
#sys.stderr.write('\n*** Identify significant islands using FDR criterion\n')
#islandFdr= os.path.join(outputDir, basename + '.island-fdr')
#cmd= """%(python)s %(script)s -i %(islandSig)s -p %(fdr)s -c 7 -o %(islandFdr)s""" \
#            %{'python': python, 
#              'script': os.path.join(srcDir, 'filter_islands_by_significance.py'), 
#              'islandSig': islandSig,
#              'fdr': args.fdr,
#              'islandFdr': islandFdr};
#
#sys.stderr.write(cmd + '\n')
#p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
#stdout, stderr= p.communicate()
#if p.returncode != 0:
#    sys.stderr.write(stderr + '\n')
#    sys.exit(p.returncode)


