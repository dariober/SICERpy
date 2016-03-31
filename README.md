## SICERpy: A friendly version of SICER peak caller

[SICER](http://home.gwu.edu/~wpeng/Software.htm) is a popular peak caller particularly suitable for detection of broad ChIP-Seq marks. 
However, I found the original master script, `SICER.sh`, clunky to use. 

`SICER.py` is a re-write of `SICER.sh` aimed at improving the following aspects:

* **Friendly API** `SICER.py` takes command line arguments using the usual and convenient format `SICER.py --option arg`. In contrast, the original script
takes eleven (!) positional arguments. Also, each step of the pipeline is checked for clean exit so that if something goes wrong you now right away.

* **Input is BAM** As opposed to `SICER.sh` which requires bed format. Options to filter reads by flag and mapping quality are also provided.
To convert bed to bam see [bedtools bedtobam](http://bedtools.readthedocs.org/en/latest/content/tools/bedtobam.html)

* **Parallel execution of multiple SICER runs** The original script caused temporary files from multiple instances of SICER to overwrite each other. 
`SICER.py` instead uses unique temporary directories so it can be run in parallel in the same directory from the same input files.

* **Less redundant output** By default, the only output produced is the table of candidate islands with their statistical significance.

* **Faster execution** Some steps of the pipeline have been parallelised.

## Requirements and Installation

`SICER.py` requires python 2.6+ with `scipy` package as per the original version. 
In addition it requires the [pysam](http://pysam.readthedocs.org/en/latest/) package.

To install, first download the `SICERpy` directory and change permission of `SICER.py` to be executable: `chmod 755 SICER.py`.
Then use sicer in one of the following ways:

* *Best*: Create a symlink from `SICER.py` to a directory somewhere on your PATH and execute `SICER.py [options]`, e.g. 

```
ln -s /path/to/SICERpy/SICER.py ~/bin/
``` 

* Call using full path `/path/to/SICERpy/SICER.py [options]`
* Add the whole directory `SICERpy` to your PATH variable.

## Usage examples

These input files can be found in the `ex/` directory.

```
SICER.py -t ex/test.bam -c ex/control.bam -s hg19 -rt 0 > peaks.bed 2> sicer.log
```

Here, the table of candidate islands is sent to `peaks.bed`. The log is captured by stderr. 
Note that by setting the redundancy threshold to 0 we use directly the input files without further filtering. 
This is useful if the bam file have been already de-duplicated with external tools like picard/MarkDuplicates.

Apply some filtering to discard reads with low mapping quality and with certain bits set (see [explain samflag](https://broadinstitute.github.io/picard/explain-flags.html)):

```
SICER.py -t ex/test.bam -c ex/control.bam -s hg19 -F 3972 -q 5 > peaks.bed
```

**Important** if working with paired-end reads discard the second-in-pair to avoid double counting!

The output is in bed format with columns:

* chrom
* start
* end
* ChIP island read count
* CONTROL island read count
* p value
* Fold change
* FDR threshold

The output file `peaks.bed` can be easily parsed to get the enriched regions. For example, the get regions with FDR < 0.01:

```
awk '$8 < 0.01' peaks.bed > peaks.01.bed
```

## Help on options

*The one shown here might be out of date*

```
SICER.py -h
usage: SICER.py [-h] [--treatment TREATMENT] --control CONTROL --species
                SPECIES [--effGenomeSize EFFGENOMESIZE]
                [--requiredFlag REQUIREDFLAG] [--filterFlag FILTERFLAG]
                [--mapq MAPQ] [--redThresh REDTHRESH]
                [--windowSize WINDOWSIZE] [--gapSize GAPSIZE]
                [--fragSize FRAGSIZE] [--keeptmp] [--version]

DESCRIPTION

Run the SICER pipeline using a control and an input file.
    
SEE ALSO:

https://github.com/dariober/SICERpy
    

optional arguments:
  -h, --help            show this help message and exit
  --treatment TREATMENT, -t TREATMENT
                        Treatment (pull-down) file in bam format
                                           
  --control CONTROL, -c CONTROL
                        Control (input) file in bam format
                                           
  --species SPECIES, -s SPECIES
                        Species to use. See or edit lib/GenomeData.py for available species. 
                                           
  --effGenomeSize EFFGENOMESIZE, -gs EFFGENOMESIZE
                        Effective Genome as fraction of the genome size. It depends on read length. Default 0.74.
                                           
  --requiredFlag REQUIREDFLAG, -f REQUIREDFLAG
                        Keep reads with these bits set in flag. Same as `samtools view -f`. Default 0
                                           
  --filterFlag FILTERFLAG, -F FILTERFLAG
                        Discard reads with these bits set in flag. Same as `samtools view -F`. Default 4. 
                        You probably want to discard also second-in-pair reads, secondary and supplementary alignments, reads failing QC.
                                           
  --mapq MAPQ, -q MAPQ  Discard reads with mapping quality lower than this. Default 5.
                                           
  --redThresh REDTHRESH, -rt REDTHRESH
                        Redundancy threshold to keep reads mapping to the same position on the same strand. Default 0 (do not filter for redundancy). 
                                           
  --windowSize WINDOWSIZE, -w WINDOWSIZE
                        Size of the windows to scan the genome. WINDOW_SIZE is the smallest possible island. Default 200.
                                           
  --gapSize GAPSIZE, -g GAPSIZE
                        Multiple of window size used to determine the gap size. Must be an integer. Default: 3.
                                           
  --fragSize FRAGSIZE, -fs FRAGSIZE
                        Size of the sequenced fragment. The center of the the fragment will be taken as half the fragment size. Default 150.
                                           
  --keeptmp             For debugging: Do not delete temp directory at the end of run.
                                           
  --version             show program's version number and exit
```

## See also

* Original files: http://home.gwu.edu/~wpeng/Software.htm

Original papers describing SICER: 

* [Spatial Clustering for Identification of ChIP-Enriched Regions (SICER) to Map Regions of Histone Methylation Patterns in Embryonic Stem Cells](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4152844/)
* [A clustering approach for identification of enriched domains from histone modification ChIP-Seq data](http://bioinformatics.oxfordjournals.org/content/25/15/1952.full)

* Google group for SICER: https://groups.google.com/forum/#!forum/sicer-users

## TODO

* Clean up the mess inside SICERpy. Lots of original code and scripts are not necessary anymore.
* Handle cases where window size is larger than chrom size. Currently
  the solution is to create new genomes without the small chromosomes.
* In fact, you don't need the genome data at all since it can be
  extracted from the bam header!
* Re-write the other `SICER*.sh` scripts

