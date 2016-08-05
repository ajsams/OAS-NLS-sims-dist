#!/usr/local/bin/python
#----------------------------------------------------------------------------
# DEPENDENCIES
#----------------------------------------------------------------------------

import argparse
import sys
import os
import time
import gzip

from math import e,log
from random import uniform, randint
from shlex import split as shsplit
from subprocess import Popen, PIPE
from collections import defaultdict, Counter

import numpy as np

from numpy.random import binomial

#----------------------------------------------------------------------------
# DEFINITIONS
#----------------------------------------------------------------------------

phelp='''

    Simulation of Neandertal-like-site Haplotype length scores under gravel model
        with Neandertal admixture.

        filename: xFNAMEx
        author: Aaron J. Sams
        date: 2016-07-07

    This program runs macs repeatedly for a single segregating site with
        the given input parameters until a locus with Neandertal-like-sites
        is generated. Once generated, sample statistics are calculated.

    Model Choices:
        one-pulse     -> Single pulse of Neandertal ancestry to common
                         ancestor of Europeans/Asians.
        two-pulse-all -> Includes first pulse plus an additional pulse into
                         both Europeans and Asians.
        two-pulse-eur -> Second pulse only in Europeans.
        two-pulse-asn -> Second pulse only in Asians.

'''

if hasattr('__main__','__file__'):
    fname = os.path.abspath(sys.modules['__main__'].__file__)
else:
    fname = 'None'
phelp=phelp.replace('xFNAMEx', fname)

subs = ('xCHRLENx','xTHETAx','xRFLAGx','xRVALx','xTOTSAMPx',
            'xAFRSAMPx','xEURSAMPx','xASNSAMPx','xNEASAMPx','xAFRN0x',
            'xEURN0x','xASNN0x','xNEANx','xAFRG1x','xEURG1x',
            'xASNG1x','xM1Ax','xM1Bx','xM1Cx','xAFRN1x',
            'xEURG2x','xASNG2x','xEURN2x','xASNN2x','xRI2Ex',
            'xRI2Ax','xM2x','xRI1x','xAFRN2x','xTG1x',
            'xTG2x','xTIE2Ex','xTIE2Ax','xTIS2Ex','xTIS2Ax',
            'xTS1x','xTIE1x','xTIS1x','xTS2x','xTG3x',
            'xTS3x')

modp1 = '''xTOTSAMPx xCHRLENx -t xTHETAx xRFLAGx xRVALx \
-I 4 xAFRSAMPx xEURSAMPx xASNSAMPx xNEASAMPx 0 \
-n 1 xAFRN0x \
-n 2 xEURN0x \
-n 3 xASNN0x \
-n 4 xNEANx \
-eg 0 1 xAFRG1x \
-eg 1e-12 2 xEURG1x \
-eg 2e-12 3 xASNG1x \
-em 3e-12 1 2 xM1Ax \
-em 4e-12 2 1 xM1Ax \
-em 5e-12 1 3 xM1Bx \
-em 6e-12 3 1 xM1Bx \
-em 7e-12 2 3 xM1Cx \
-em 8e-12 3 2 xM1Cx \
-en xTG1x 1 xAFRN1x \
-eg xTG1x 2 xEURG2x \
-eg xTG1x 3 xASNG2x \
-en xTG2x 2 xEURN2x \
-en xTG2x 3 xASNN2x \
-ej xTS1x 3 2 \
-em xTS1x 1 2 xM2x \
-em xTS1x 2 1 xM2x \
-em xTIE1x 2 4 xRI1x \
-em xTIS1x 2 4 0 \
-ej xTS2x 2 1 \
-en xTG3x 1 xAFRN2x \
-ej xTS3x 4 1 \
-h 1e3 \
'''

modp2all = '''xTOTSAMPx xCHRLENx -t xTHETAx xRFLAGx xRVALx \
-I 4 xAFRSAMPx xEURSAMPx xASNSAMPx xNEASAMPx 0 \
-n 1 xAFRN0x \
-n 2 xEURN0x \
-n 3 xASNN0x \
-n 4 xNEANx \
-eg 0 1 xAFRG1x \
-eg 1e-12 2 xEURG1x \
-eg 2e-12 3 xASNG1x \
-em 3e-12 1 2 xM1Ax \
-em 4e-12 2 1 xM1Ax \
-em 5e-12 1 3 xM1Bx \
-em 6e-12 3 1 xM1Bx \
-em 7e-12 2 3 xM1Cx \
-em 8e-12 3 2 xM1Cx \
-en xTG1x 1 xAFRN1x \
-eg xTG1x 2 xEURG2x \
-eg xTG1x 3 xASNG2x \
-en xTG2x 2 xEURN2x \
-en xTG2x 3 xASNN2x \
-em xTIE2Ex 2 4 xRI2Ex \
-em xTIE2Ax 3 4 xRI2Ax \
-em xTIS2Ex 2 4 0 \
-em xTIS2Ax 3 4 0 \
-ej xTS1x 3 2 \
-em xTS1x 1 2 xM2x \
-em xTS1x 2 1 xM2x \
-em xTIE1x 2 4 xRI1x \
-em xTIS1x 2 4 0 \
-ej xTS2x 2 1 \
-en xTG3x 1 xAFRN2x \
-ej xTS3x 4 1 \
-h 1e3 \
'''

modp2eur = '''xTOTSAMPx xCHRLENx -t xTHETAx xRFLAGx xRVALx \
-I 4 xAFRSAMPx xEURSAMPx xASNSAMPx xNEASAMPx 0 \
-n 1 xAFRN0x \
-n 2 xEURN0x \
-n 3 xASNN0x \
-n 4 xNEANx \
-eg 0 1 xAFRG1x \
-eg 1e-12 2 xEURG1x \
-eg 2e-12 3 xASNG1x \
-em 3e-12 1 2 xM1Ax \
-em 4e-12 2 1 xM1Ax \
-em 5e-12 1 3 xM1Bx \
-em 6e-12 3 1 xM1Bx \
-em 7e-12 2 3 xM1Cx \
-em 8e-12 3 2 xM1Cx \
-en xTG1x 1 xAFRN1x \
-eg xTG1x 2 xEURG2x \
-eg xTG1x 3 xASNG2x \
-en xTG2x 2 xEURN2x \
-en xTG2x 3 xASNN2x \
-em xTIE2Ex 2 4 xRI2Ex \
-em xTIS2Ex 2 4 0 \
-ej xTS1x 3 2 \
-em xTS1x 1 2 xM2x \
-em xTS1x 2 1 xM2x \
-em xTIE1x 2 4 xRI1x \
-em xTIS1x 2 4 0 \
-ej xTS2x 2 1 \
-en xTG3x 1 xAFRN2x \
-ej xTS3x 4 1 \
-h 1e3 \
'''

modp2asn = '''xTOTSAMPx xCHRLENx -t xTHETAx xRFLAGx xRVALx \
-I 4 xAFRSAMPx xEURSAMPx xASNSAMPx xNEASAMPx 0 \
-n 1 xAFRN0x \
-n 2 xEURN0x \
-n 3 xASNN0x \
-n 4 xNEANx \
-eg 0 1 xAFRG1x \
-eg 1e-12 2 xEURG1x \
-eg 2e-12 3 xASNG1x \
-em 3e-12 1 2 xM1Ax \
-em 4e-12 2 1 xM1Ax \
-em 5e-12 1 3 xM1Bx \
-em 6e-12 3 1 xM1Bx \
-em 7e-12 2 3 xM1Cx \
-em 8e-12 3 2 xM1Cx \
-en xTG1x 1 xAFRN1x \
-eg xTG1x 2 xEURG2x \
-eg xTG1x 3 xASNG2x \
-en xTG2x 2 xEURN2x \
-en xTG2x 3 xASNN2x \
-em xTIE2Ax 3 4 xRI2Ax \
-em xTIS2Ax 3 4 0 \
-ej xTS1x 3 2 \
-em xTS1x 1 2 xM2x \
-em xTS1x 2 1 xM2x \
-em xTIE1x 2 4 xRI1x \
-em xTIS1x 2 4 0 \
-ej xTS2x 2 1 \
-en xTG3x 1 xAFRN2x \
-ej xTS3x 4 1 \
-h 1e3 \
'''

#----------------------------------------------------------------------------
# FUNCTIONS
#----------------------------------------------------------------------------

###
# Generate macs command line
###

def scale_time(tyears, yearspergen, scalingfactor):
    '''Scales time in years 'tyears', in units of generations
    per 'scalingfactor' (typically 4*N0), given a 'yearspergen'.
    '''
    return (tyears/float(yearspergen))/float(scalingfactor)

def fit_decay_rate(n0, nt, t):
    '''Provides the rate (r) of exponential decay that fits
    the input parameters.
        n0 = initial size
        nt = size at time t
        t = time
    '''
    return log(n0/float(nt))/float(t)

def split_times(ts1, ts2, ts3):
    '''If split times are not given (False), randomly sample them.
        ts1 = European/Asian split
        ts2 = African/non-African split
        ts3 = Neandertal/modern human split
    '''
    if not ts1:
        ts1 = uniform(36000,55000)
    if not ts2:
        ts2 = 70000
    if not ts3:
        ts3 = 700000
    return float(ts1), float(ts2), float(ts3)

def introgression_params(tis1, tie1, ri1, tis2e, tie2e, ri2e,
                         tis2a, tie2a, ri2a, ts1, ts2, model):
    '''If introgression times are not given, they
    are randomly sampled.  tis is sampled from
    Uniform(ts1, ts2). Unless values are specified,
    there is always at least 500 years of introgression.
    If ri (rate of introgression) not specified, sample.
    '''
    # First pulse always happens
    if not isinstance(tis1, float):
        #+1 and minus 1 ensure that time of events in macs are not identical.
        #uniform(ts1 + 1, min([65000,ts2]) - 1)
        tis1 = uniform(ts1,65000)
    if not isinstance(tie1, float):
        tie1 = tis1 - 500
    if not isinstance(ri1, float):
        ri1 = 7.5e-4

    # Second 500 year pulse to Europe
    if model in ['two-pulse-all','two-pulse-eur']:
        if not isinstance(tis2e, float):
            tis2e = ts1 - 500
        if not isinstance(tie2e, float):
            tie2e = tis2e - 500
        if not isinstance(ri2e, float):
            ri2e = uniform(0, 1.5e-5)

    # Second pulse to Asia
    if model in ['two-pulse-all','two-pulse-asn']:
        if not isinstance(tis2a, float):
            tis2a = ts1 - 500
        if not isinstance(tie2a, float):
            tie2a = tis2a - 500
        if not isinstance(ri2a, float):
            ri2a = uniform(0, 2.5e-4)
    return tis1, tie1, ri1, tis2e, tie2e, ri2e, tis2a, tie2a, ri2a

def eurasian_N(eurasnN, eurasnNratio, eurasnNsum):
    '''If size of ancestral eurasian pop 'eurasnN' not specified,
    samples.  Also samples initial European N 'eurN2' if not given.
    '''
    if not eurasnN:
        eurasnN=uniform(5000,15000)
    if not eurasnNratio:
        eurasnNratio=uniform(1.0,2.5)
    if not eurasnNsum:
        eurasnNsum=uniform(8000,25000)
    eurN2 = eurasnNratio * eurasnNsum
    asnN2 = eurasnNsum/float(eurasnNratio + 1.0)
    return eurasnN, eurN2, asnN2

def modstring(model):
    '''Return the parameter substitution string for the
    specified model.
    '''
    if model == 'two-pulse-all':
        return modp2all
    elif model == 'two-pulse-eur':
        return modp2eur
    elif model == 'two-pulse-asn':
        return modp2asn
    else:
        return modp1

def substitute(model, subs, params):
    '''Creates the string input to use with subprocess
    given input parameters.
    '''
    times = set([])
    timeparams = ('xTG1x','xTG2x','xTIE2Ex','xTIE2Ax','xTIS2Ex','xTIS2Ax',
                    'xTS1x','xTIE1x','xTIS1x','xTS2x','xTG3x','xTS3x')
    input = modstring(model)
    sublist = [[subs[i],params[i]] for i in range(len(subs))]
    for i,[sub,param] in enumerate(sublist):
        if sub in timeparams:
            while True:
                lastinput = input
                if param in times:
                    param += 1e-8
                    times.add(param)
                else:
                    times.add(param)
                input = input.replace(sub,str(param),1)
                if lastinput == input:
                    break
        else:
            input = input.replace(sub,str(param))
    return input

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description=phelp)
    p.add_argument('range', type = int, nargs = 2,
                    help='First and last index (inclusive) to use in \
                    out and log files.')
    p.add_argument('--hotspotfile', '-r', type=str, default=False,
                    help='Path to recombination rate file that specifies \
                    cM/bp in intervals on the chromosome scaled by total \
                    length of chromosome.')
    p.add_argument('--chromlen','-c', type=int, default=600000,
                    help='Length of simulated chromosome in base pairs.')
    p.add_argument('--chrstart','-s',type=float, default=0.4166666,
                    help='Start of focal region (as fraction of total chrom len.')
    p.add_argument('--chrend','-d',type=float, default=0.6333334,
                    help='End of focal region (as fraction of total chrom len.')
    p.add_argument('--thin','-t', type = float, default=0.0535,
                    help='Fraction of SNPs to be randomly thinned out of dataset.')
    p.add_argument('--minfreq', type = float, default=0.28,
                    help='Minimum frequency NLS to gather.')
    p.add_argument('--outtag','-o', type=str,
                    default='OAS_haps', help='Output file tag. \
                    Defaults to OAS_haps.')
    p.add_argument('--mutrate','-m', type=float, default=2.5e-08,
                    help='Neutral mutation rate for simulations.')
    p.add_argument('--recrate', type=float, default=1.3e-08,
                    help='Recombination rate for simulations.')
    p.add_argument('--afrsamp', type=int, default=216,
                    help='Haploid sample size for African sample.')
    p.add_argument('--eursamp', type=int, default=198,
                    help='Haploid sample size for European sample.')
    p.add_argument('--asnsamp', type=int, default=0,
                    help='Haploid sample size for Asian sample.')
    p.add_argument('--neasamp', type=int, default=2,
                    help='Haploid sample size for Neandertal sample.')
    p.add_argument('--tsplit1', type=float, default=False,
                    help='Time (in years) from present to split between \
                    European/Asian pops.')
    p.add_argument('--tsplit2', type=float, default=False,
                    help='Time (in years) from present to split between \
                    Eurasian/African pops.')
    p.add_argument('--tsplit3', type=float, default=False,
                    help='Time (in years) from present to split between \
                    African/Neandertal pops.')
    p.add_argument('--tintstart1', type=float, default=False,
                    help='Time (in years) from present to start of \
                    introgression.')
    p.add_argument('--tintend1', type=float, default=False,
                    help='Time (in years) from present to end of \
                    introgression.')
    p.add_argument('--intrate1', type=float, default=False,
                    help='Fraction of Non-African population per \
                    generation made up of Neandertals during introgression.')
    p.add_argument('--tintstart2eur', type=float, default=False,
                    help='Time (in years) from present to start of \
                    introgression 2 into Europe.')
    p.add_argument('--tintend2eur', type=float, default=False,
                    help='Time (in years) from present to end of \
                    introgression 2 into Europe.')
    p.add_argument('--intrate2eur', type=float, default=False,
                    help='Fraction of European population per generation \
                    made up of Neandertals during introgression.')
    p.add_argument('--tintstart2asn', type=float, default=False,
                    help='Time (in years) from present to start of \
                    introgression 2 into Asia.')
    p.add_argument('--tintend2asn', type=float, default=False,
                    help='Time (in years) from present to end of \
                    introgression 2 into Asia.')
    p.add_argument('--intrate2asn', type=float, default=False,
                    help='Fraction of Asian population per generation \
                    made up of Neandertals during introgression.')
    p.add_argument('--eurasnN', type=float, default=False,
                    help='Initial population size of the combined \
                    eurasian population.')
    p.add_argument('--eurasnNratio', type=float, default=False,
                    help='Ratio of eurNe / asnNe after split.')
    p.add_argument('--eurasnNsum', type=int, default=False,
                    help='Sum of eurNe and asnNe after split.')
    p.add_argument('--model', type=str, default='one-pulse',
                    choices=['one-pulse','two-pulse-all', 'two-pulse-eur',
                    'two-pulse-asn'],
                    help='Define model to use for simulation.')
    p.add_argument('--local', action='store_true', default=False,
                    help='If flagged, uses executable instead of \
                    ~/OAS-NLS-sims-dist/bin/executable')
    p.add_argument('--focpop', type=str, choices = ['EUR','ASN'], default = 'EUR',
                    help='Focal population for which to extract results.')
    p.add_argument('--focregstart', type=int, default=113100000,
                    help='True position of start of simulated chromosome.')
    args = p.parse_args()
    return args

def generate_macs_input(args):
    '''Uses input arguments to generate the input for macs.'''
    ###
    # Scaling parameters
    ###
    Nscale = float(7310)                          # Initial Ne
    unit = 4*Nscale                               # Scaling unit
    yearspergen = 25.0
    theta0 = float(unit*args.mutrate)             # theta-zero
    rho0 = float(unit*args.recrate)
    totsamp = sum([args.afrsamp, args.eursamp, args.asnsamp, args.neasamp])

    ###
    # Unscaled sizes, times, and growth rates
    ###
    # Pop sizes
    afrN0 = 424000                                # pop 1 (time 0)
    afrN1 = 14474                                 # pop 1 (after pop decline at tg1)
    afrN2 = 7310                                  # pop 1 (initial)
    eurN0 = 512000                                # pop 2 (time 0)
    asnN0 = 1370990                               # pop 3 (time 0)
    neaN = 1500                                   # pop 4 (time 0)
    eurasnN, eurN2, asnN2 = eurasian_N(args.eurasnN,   # OOA, EUR, and ASN N after split.
                                       args.eurasnNratio, args.eurasnNsum)
    eurN1 = max([eurN2, 9475])                    # eur and asn Ns after
    asnN1 = max([asnN2, 8879])                    # first growth (at 5115)

    # Event times (g = growth, s = split, i = introgression)
    tg1 = 5115
    tg2 = 23000
    tg3 = 148000              # time African pop doubles
    ts1, ts2, ts3 = split_times(args.tsplit1, args.tsplit2, args.tsplit3)
    tis1, tie1, ri1, tis2e, tie2e, ri2e, tis2a, tie2a, ri2a = introgression_params(
                         args.tintstart1, args.tintend1, args.intrate1,
                         args.tintstart2eur, args.tintend2eur, args.intrate2eur,
                         args.tintstart2asn, args.tintend2asn, args.intrate2asn,
                         ts1, ts2, args.model)
    # Migration rates
    m1a = 2.5e-5 * (23000/ts1)                    # migration bt Africa & Europe
    m1b = 7.8e-6 * (23000/ts1)                    # migration bt Africa & Asia
    m1c = 3.11e-5 * (23000/ts1)                   # migration bt Europe & Asia
    m2 = 1.5e-4                                   # migration bt Africa & Eurasia

    ###
    # Ne and Unit (4Ne) Scaled Parameters
    ###

    afrN0s = afrN0/Nscale
    afrN1s = afrN1/Nscale
    afrN2s = afrN2/Nscale
    eurN0s = eurN0/Nscale
    eurN1s = eurN1/Nscale
    eurN2s = eurN2/Nscale
    asnN0s = asnN0/Nscale
    asnN1s = asnN1/Nscale
    asnN2s = asnN2/Nscale
    neaNs = neaN/Nscale

    tg1s = scale_time(tg1, yearspergen, unit)
    tg2s = scale_time(tg2, yearspergen, unit)
    tg3s = scale_time(tg3, yearspergen, unit)
    ts1s = scale_time(ts1, yearspergen, unit)
    ts2s = scale_time(ts2, yearspergen, unit)
    ts3s = scale_time(ts3, yearspergen, unit)
    tis1s = scale_time(tis1, yearspergen, unit)
    tie1s = scale_time(tie1, yearspergen, unit)
    tis2es = scale_time(tis2e, yearspergen, unit)
    tie2es = scale_time(tie2e, yearspergen, unit)
    tis2as = scale_time(tis2a, yearspergen, unit)
    tie2as = scale_time(tie2a, yearspergen, unit)

    m1as = m1a*unit
    m1bs = m1b*unit
    m1cs = m1c*unit
    m2s = m2*unit

    # Growth rates
    afrG1 = fit_decay_rate(afrN0s, afrN1s, tg1s)
    eurG1 = fit_decay_rate(eurN0s, eurN1s, tg1s)
    asnG1 = fit_decay_rate(asnN0s, asnN1s, tg1s)
    eurG2 = 0.0
    asnG2 = 0.0
    if eurN1s > eurN2s: #if sampled eurN2 is smaller than later eurN1s
        eurG2 = fit_decay_rate(eurN1s, eurN2s, tg2s)
    if asnN1s > asnN2s:
        asnG2 = fit_decay_rate(asnN1s, asnN2s, tg2s)
    inputa = '/'.join(os.getcwd().split('/')[:-1])+'/bin/macs '
    if args.local:
        inputa = 'macs '

    # Configure recombination
    if args.hotspotfile:
        rval = str(args.hotspotfile)
        rflag = '-R'
    else:
        rval = str(rho0)
        rflag = '-r'

    # Generate macs input string
    params = (args.chromlen, theta0, rflag, rval, totsamp,
              args.afrsamp, args.eursamp, args.asnsamp, args.neasamp, afrN0s,
              eurN0s, asnN0s, neaNs, afrG1, eurG1,
              asnG1, m1as, m1bs, m1cs, afrN1s,
              eurG2, asnG2, eurN2s, asnN2s, ri2e,
              ri2a, m2, ri1, afrN2s, tg1s,
              tg2s, tie2es, tie2as, tis2es, tis2as,
              ts1s, tie1s, tis1s, ts2s, tg3s,
              ts3s)

    inputb = substitute(args.model, subs, params)
    input = inputa+inputb
    return input

###
# Probabilistically thin data
###

def build_prob_dict(totaln):
    '''Builds dictionary for probability thinning.
        totaln: Total number of haplotypes in a sample.
    '''
    pdict={}
    for i in range(totaln+1):
        if i in [0, totaln]:
            pdict[i] = 0.0 #filters out non-segregating sites in sample
        elif i in [1, totaln-1]:
            pdict[i] = 0.25
        elif i in [2, totaln-2]:
            pdict[i] = 0.5
        elif i in [3, totaln-3]:
            pdict[i] = 0.75
        elif i in [4, totaln-4]:
            pdict[i] = 0.8
        elif i in [5, totaln-5]:
            pdict[i] = 0.9
        elif i in [6, totaln-6]:
            pdict[i] = 0.95
        elif i in [7, totaln-7]:
            pdict[i] = 0.96
        elif i in [8, totaln-8]:
            pdict[i] = 0.97
        elif i in [9, totaln-9]:
            pdict[i] = 0.98
        else:
            pdict[i] = 0.99
    return pdict

def prob_thin_mask(counts,pdict,rnd_pct=0.05):
    '''Generate a boolean mask given:
        pdict: Probability dictionary from build_prob_dict.
        counts: A list of derived allele counts.
        rnd_pct: Random probability of thinning out.
    Returned outmask[i] is True if site i should be kept.
    Can be added (with and) to nlsmask to generate a final mask.
    Note: This function also filters out non-segregating sites.
    '''
    cond = lambda x: binomial(1,x)
    cond2 = lambda: binomial(1,1.0-rnd_pct)
    test = lambda x: True if cond(x) and cond2() else False
    outmask = []
    for i,count in enumerate(counts):
        outmask.append(test(pdict[count]))
    return np.array(outmask,dtype=bool)

###
# Run ms command and capture output
###


def run_macs2ms(input, local=False):
    '''Run ms 'input' in the shell and gather input.
        input: macs command line string.
        local: Specifies if local machine has msformatter
            in path (True) or in directory (False).
    '''
    msformat = '/'.join(os.getcwd().split('/')[:-1])+'/bin/msformatter '

    if local:
        msformat = 'msformatter'
    args = shsplit(input)
    args2 = shsplit(msformat)
    p1 = Popen(args, stdout=PIPE, stderr=PIPE)
    p2 = Popen(args2, stdin = p1.stdout, stdout=PIPE, stderr=PIPE)
    p1.stdout.close()
    p1.stderr.close() #This stderr is the error output from macs, which gives the trees.
    msout, mserr = p2.communicate()
    return msout, mserr

def parse_ms_output(stdout,afrsz,eursz,asnsz,neasz,focpdict,focal='EUR',
                    stt=0.0,stp=1.0,rnd_pct=0.05,minfreq=0.2):
    '''Parse the ms stdout.  Needs:
        stdout: string from subprocess.
        afrsz: African sample size.
        eursz: European sample size.
        asnsz: Asian sample size.
        focpdict: Probability thinning dict for focal population.
        focal: Focal population for NLS analysis (EUR or ASN)
        stt: Focal region start (on uniform chr)
        stp: Focal region end (on uniform chr)
        rnd_pct: Random thinning percentage.
        minfreq: Smallest frequency NLS to process.
    Returns header and 5 sets of tuples for (eur, asn):
        headerd: Dictionary of items from header.
        POPhaps: Haplotypes from each population.
        POPcts: Counts of derived alleles for each site in population.
        POPtest: True if pop passed filters for downstream analysis, else False.
        POPjointmask: Boolean mask, True if site passes all filters for focal site analysis.
        POPprobmask: Boolean mask, True if site passed probability thinning filter.
    '''
    # Make header dictionary
    sout = stdout.split('\n')
    if len(sout) < 6:
        return (False, False, False, False, False, False)
    outit = iter(sout)
    headerd = {'input': outit.next(), 'seed': outit.next(), 'blank': outit.next(),
                'hash': outit.next(), 'segs': int(outit.next().split()[1]),
                'positions': np.array([float(p) for p in outit.next().split()[1:]])}
    # Focal positions mask True for those that are in focal region
    headerd['subpos-mask'] = np.ma.masked_inside(headerd['positions'],stt,stp).mask

    # Initialize hap lists
    afrhaps = []
    eurhaps = []
    asnhaps = []
    neahaps = []

    scount = headerd['segs']
    # Add values to arrays
    for i in range(afrsz):
        line = np.fromiter(outit.next(), dtype=int, count = scount)
        afrhaps.append(line)
    for i in range(eursz):
        line = np.fromiter(outit.next(), dtype=int, count = scount)
        eurhaps.append(line)
    for i in range(asnsz):
        line = np.fromiter(outit.next(), dtype=int, count = scount)
        asnhaps.append(line)
    for i in range(neasz):
        line = np.fromiter(outit.next(), dtype=int, count = scount)
        neahaps.append(line)

    afrhaps = np.array(afrhaps)
    neahaps = np.array(neahaps)
    assert focal in ['EUR','ASN'], 'Focal population is mis-specified!'
    if focal == 'EUR':
        fochaps = np.array(eurhaps)
    else:
        fochaps = np.array(asnhaps)

    # Derived allele counts
    afrcts = np.sum(afrhaps, axis=0)
    neacts = np.sum(neahaps, axis=0)
    foccts = np.sum(fochaps, axis=0)

    # Frequency masks
    if focal == 'EUR':
        focsz = eursz
    else:
        focsz = asnsz
    focfmask = np.ma.masked_greater_equal(foccts, minfreq * focsz).mask

    # Test for NLS at each segsite for each pop
    afrnlscond = np.ma.masked_equal(afrcts, 0.0).mask
    neanlscond = np.ma.masked_greater(neacts, 0.0).mask
    focnlscond = np.ma.masked_greater(foccts, 0.0).mask

    focnlsmask = reduce(np.logical_and, [afrnlscond, focnlscond, neanlscond])

    # Probability thinning and segsite mask
    focprobmask = prob_thin_mask(foccts, focpdict)

    # Combine all masks to find positions that pass all
    focjointmask = reduce(np.logical_and, [focfmask, focnlsmask, focprobmask, headerd['subpos-mask']])
    foctest = any(focjointmask)

    if not foctest:
        return (False, False, False, False, False, False)

    return (foctest, headerd, fochaps, foccts, focjointmask, focprobmask)

def subpop_ms(headerd, pophaps, popcts, popjointmask, popprobmask):
    '''Convert output for single pop from parse_ms_output
    to a new ms file string to use in sample_stats.
        headerd: header dict from parse_ms_output.
        pophaps: POPhaps list from parse_ms_output.
        popcts: POPcts dict from parse_ms_output.
        popjointmask: Boolean mask for sites that passed all filters.
        popprobmask: Boolean mask for sites that passed segsite/thinning filter.
    Returns:
        msout: ms output file formatted string to use in
            ms_to_samplestats.
        msderouts: A list of ms outputs for derived haps at each site
            that passed joint mask.
        msancouts: A list of ms outputs for ancestral haps at each site
            that passed joint mask.
        posouts: A list of NLS positions in H-scan format (converted to
            1e6 chr len).
        dercts: List of derived allele counts for each NLS position.
        i: Segregating site count in subpopulation sample.
    '''
    n = len(pophaps)
    i = np.sum(popprobmask)
    msout = []
    # Filter all data by prob-thinning/segsite mask
    mprobmask = np.array([popprobmask]*n)
    subpophaps = pophaps[mprobmask].reshape(n,i)
    jprobmask = popjointmask[popprobmask]
    subpopcts = popcts[popprobmask]
    positions = headerd['positions'][popprobmask]

    # Re-format new header and generate full subpop ms output
    newsegs = 'segsites: '+str(i)
    newpositions = 'positions: '+' '.join([str(p) for p in positions])
    fullhaps = [''.join(np.array(hap,dtype=str)) for hap in subpophaps]
    msout = '\n'.join([headerd['input'],headerd['seed'],headerd['blank'],
                headerd['hash'],newsegs,newpositions]+fullhaps)

    # Generate ms output(s) for each site that passed joint mask

    msderouts = []
    msancouts = []
    posouts = [] #positions of NLS in H-scan format
    dercts = []
    for index,test in enumerate(jprobmask):
        if test:
            # Store H-scan version of position
            posouts.append(int(np.floor(positions[index]*1e6))) #important to use np.round

            # Generate der and anc hap sets
            sortedhaps = sorted(subpophaps, key=lambda x: x[index], reverse=True)
            derhaps = np.array(sortedhaps[:subpopcts[index]])
            anchaps = np.array(sortedhaps[subpopcts[index]:])

            # Output ms output without filtering to segsites in sub-haps (derhaps, anchaps)
            nd = len(derhaps)
            na = len(anchaps)
            dercts.append(nd)

            # Generate der and anc ms outputs
            dersegs = 'segsites: '+str(i)
            ancsegs = 'segsites: '+str(i)
            derpositions = 'positions: '+' '.join([str(p) for p in positions])
            ancpositions = 'positions: '+' '.join([str(p) for p in positions])

            fullderhaps = [''.join(np.array(hap,dtype=str)) for hap in derhaps]
            fullanchaps = [''.join(np.array(hap,dtype=str)) for hap in anchaps]

            derout = '\n'.join([headerd['input'],headerd['seed'],headerd['blank'],
                headerd['hash'],dersegs,derpositions]+fullderhaps)
            ancout = '\n'.join([headerd['input'],headerd['seed'],headerd['blank'],
                headerd['hash'],ancsegs,ancpositions]+fullanchaps)
            msderouts.append(derout)
            msancouts.append(ancout)
    return msout, msderouts, msancouts, posouts, dercts, i

###
# Get stats from ms outputs
###

def ms_to_samplestats(msout, local=False):
    '''Converts eurms and asnms objects to files and pipes
    to sample_stats program to calculate sample statistics.
        msout: eurms or asnms from parse_ms_output
        local: Specifies if local machine has sample_stats
            in path (True) or in directory (False).
    Returns:
        pi, seg sites, Fay and Wu's H, Tajima's D, Heterozygosity.
    '''
    sformat = '/'.join(os.getcwd().split('/')[:-1])+'/bin/sample_stats '

    if local:
        sformat = 'sample_stats'
    p1 = Popen(sformat, stdin = PIPE, stdout = PIPE, stderr = PIPE)
    out, err = p1.communicate(input=msout)
    out = out.strip().split('\t')
    pi, ss, td, fwh, h = [out[i] for i in range(1,10,2)]
    return pi, ss, td, fwh, h

def safelog(x, base=e, default=-float("inf")):
    try:
        return log(x, base)
    except (OverflowError, ValueError):
        return default

def ms_to_hscan(msout, focalpos, local=False):
    '''Use ms output to run H-scan.
        msout: ms output (e.g. from subpop ms)
        focalpos: Focal position of SNP on 1mb chrom.
             NOTE: Should multiply ms position by 1 million
                    and round to nearest integer.
    '''
    # Write msout to temporary file
    temp = os.getcwd()+'/temp.'+str(randint(1e8, 999999999))+'.txt'
    f = open(temp,'w')
    f.write(msout)
    f.close()
    # Run H-scan and collect output
    hformat = '/'.join(os.getcwd().split('/')[:-1])+'/bin/H-scan '

    if local:
        hformat = 'H-scan'
    args1 = shsplit(hformat+' -m -i '+temp+' -l '+str(focalpos)+' -r '+str(focalpos))
    args2 = shsplit('tail -n1')
    p1 = Popen(args1, stdout=PIPE, stderr=PIPE)
    p2 = Popen(args2, stdin = p1.stdout, stdout=PIPE, stderr=PIPE)
    p1.stdout.close()
    p1.stderr.close() #This stderr from H-scan
    hout, herr = p2.communicate()
    os.remove(temp)
    return hout.split('\t'), herr #(pos, H), error

def all_hpairs(derouts, ancouts, posouts, dercts, chrlen, focregstart, local=False):
    '''Get H-scan values for ancestral and derived haps
    at all NLS positions in 'posouts'."
        derouts: List of derived hap ms outputs per site.
        ancouts: List of ancestral hap ms outputs per site.
        posouts: List of positions of NLS.
        dercts: List of derived allele counts of NLS.
        chrstt: Start position on uniform chrom of focal region.
        chrend: End position on uniform chrom of focal region.
        focregstart: True start position in OAS region on chr12.
        chrlen: Length of simulated chromosome in bp.
    Returns:
        hpairs: List of (Hder, Hanc, ln(Hder/Hanc)) values.
    '''
    hpairs = []
    convlen = lambda x: str(int(np.round(float(x) * (chrlen/1000000.0)))) # Rescale H output.
    convpos = lambda pos: focregstart + int(np.round((pos/1e6)*chrlen)) # Rescale H position to chromosome length.
    for i,position in enumerate(posouts):
        newpos = convpos(position)
        (pos, hder),err = ms_to_hscan(derouts[i], position, local)
        (pos, hanc),err = ms_to_hscan(ancouts[i], position, local)
        hpairs.append([str(newpos), str(dercts[i]), convlen(hder), convlen(hanc), str(safelog(float(hder)/float(hanc)))])
    return hpairs

###
# Generate sims until NLS in focal region, process stats, write stats.
###


def grab_nls(args):
    '''Run macs until a Neandertal-like site (NLS)
    is generated.
    Input:
        args: arguments from get_args()
        seed: macs random seed (optional)

    Returns:
        POPparams: [seed, #sim tries, segregating-site-count,
                    der allele count, pi, segsites, Tajima's D,
                    Fay and Wu's H, thetaH-pi, macs input] for POP output.
        POPhpairs: (Hder, Hanc) for each NLS in valid region for each pop.
    '''
    if args.focpop == 'EUR':
        pid = 0
        pdict = build_prob_dict(args.eursamp)
    else:
        pid = 1
        pdict = build_prob_dict(args.asnsamp)
    tries = 0
    done = False

    while not done:
        macsin = generate_macs_input(args)
        msoutf = run_macs2ms(macsin, args.local)[0]
        params = parse_ms_output(msoutf, args.afrsamp,
                                    args.eursamp,
                                    args.asnsamp,
                                    args.neasamp,
                                    pdict,args.focpop,
                                    stt=args.chrstart,stp=args.chrend,
                                    rnd_pct=args.thin,minfreq=args.minfreq)
        # process Focal
        if params[0]:
            msout, derouts, ancouts, hpositions, dercts, segsites = subpop_ms(params[1],
                                                                    params[2],
                                                                    params[3],
                                                                    params[4],
                                                                    params[5])
            stats = list(ms_to_samplestats(msout, args.local))
            hpairs = all_hpairs(derouts, ancouts, hpositions, dercts,
                                    args.chromlen, args.focregstart,
                                    args.local)
            params = [params[1]['seed'],tries,
                        segsites]+stats+[params[1]['input']]
            done = True

        else:
            tries += 1
    return params, hpairs, msoutf

#----------------------------------------------------------------------------
# RUN
#----------------------------------------------------------------------------

if __name__=="__main__":
    args = get_args()
    outf = open(args.outtag+'.hstats.txt','w')
    logf = open(args.outtag+'.stats-params.txt','w')
    msof = gzip.open(args.outtag+'.msout.txt.gz','w')

    outhead = '\t'.join(['#sim-index','nls-index','Position','DAC','Hd','Ha','Hda'])
    loghead = '\t'.join(['#index','ms-seed','tries','segsites',
                        'pi','ss','D','thetaH','thetaH-pi','ms-input'])
    print >> outf, outhead
    print >> logf, loghead

    for i in range(args.range[0], args.range[1]+1):
        params,hpairs,msout = grab_nls(args)
        print >> logf, '\t'.join([str(i)]+[str(p) for p in params])
        print >> msof, msout
        for j,hvals in enumerate(hpairs):
            print >> outf, '\t'.join([str(i),str(j)]+hvals)

    outf.close()
    logf.close()
    msof.close()

#----------------------------------------------------------------------------
# NOTES
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
