# OAS-NLS-sims-dist
Code used in Sams et al. 2016 to generate distributions of Neandertal-like-sites
    under a neutral coalescent model with Neandertal introgression. It can be used
    on linux/unix platforms.

Sams et al. 2016 is currently on biorxiv at http://biorxiv.org/content/early/2016/05/04/051466

This package includes two main programs in the /scripts directory.

1. nea_gravel_freq_sims.py - This program generates the distributions of derived
                                    allele counts at Neandertal-Like-Sites described
                                    in Sams et al. 2016.

2. nea_gravel_haps_sims.py - This program generates the distributions of the H(D/A)
                                    statistic at Neandertal-Like-Sites described
                                    in Sams et al. 2016.

## Options

The options for each program can be accessed with the `-h` flag:

```
./nea_gravel_freq_sims.py -h
```

Most options should be straightforward to interpret, but please contact Aaron Sams
    with any issues or questions.

## Example use

To use, first cd to the /scripts directory:

```
cd scripts
```
Analyses must be run from here unless dependencies are in your PATH, in which case see below.


To generate a set of 1000 Neandertal-Like-Site derived allele counts for Europeans and
    East Asians, simply run the frequency program with default parameters, specifying
    only the index values of the results.

```
./nea_gravel_freq_sims.py 1 1000
```

Note that nea_gravel_freq_sims.py will by default write two output files to stderr and
    stdout.

In contrast, nea_gravel_haps_sims.py requires the user to specify an output file prefix
    with the `--outtag` flag.

## Supporting files

Several files specific to the analyses in Sams et al. 2016 are included in the
    /supporting directory.

These include:

1. oas_recrates.txt - The file used to specify a recombination map in haplotype simulations.
2. 1000g_OAS_NLS_freqs.txt - A file with frequencies for all Neandertal-Like-Sites in the
    OAS gene region analyzed in our paper.
3. CEU_OAS_H_results.txt - A file with H-scan results for all Neandertal-Like-Sites in the
    OAS gene region analyzed in our paper. H(D/A) is not specified in this file but can
    be calculated with the values in the Hder and Hanc columns.

These files are provided merely for those who may wish to replicate our results or explore
    simulation parameters that we did not.

## Software Dependencies

This package includes an empty /bin directory which you should fill with
    five executables. Alternatively, these executables can be in your machine's PATH,
    in which case, you must use the `--local` flag.

From the [**ms**](http://home.uchicago.edu/rhudson1/source/mksamples.html) software package (Hudson, 2002):

1. ms
2. sample_stats

From the [**macs**](https://github.com/gchen98/macs) software package (Chen et al. 2009):

3. macs
4. msformatter

Finally (5) **H-scan** (Messer), which can be downloaded from:
    https://messerlab.org/resources/

Additionally, the haplotype simulator relies on **numpy** (http://www.numpy.org),
    so your python distribution must have it installed as well.

### References:

1. Sams, AJ et al. (2016). Adaptively introgressed Neandertal haplotype at the OAS locus functionally impacts innate immune responses in humans. Biorxiv http://dx.doi.org/10.1101/051466
2. Hudson, R (2002). Generating samples under a Wright–Fisher neutral model of genetic variation. Bioinformatics.
3. Chen, GK et al. (2009). Fast and flexible simulation of DNA sequence data. Genome Research, 19(1), 136–142. http://doi.org/10.1101/gr.083634.108
4. Messer, PW. H-scan software is available from the Messer lab https://messerlab.org/resources/
