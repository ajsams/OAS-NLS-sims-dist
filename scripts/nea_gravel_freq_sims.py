#!/usr/local/bin/python
#----------------------------------------------------------------------------
# DEPENDENCIES
#----------------------------------------------------------------------------

import argparse
import sys
import os

from math import log
from random import uniform
from shlex import split as shsplit
from subprocess import Popen, PIPE

#----------------------------------------------------------------------------
# DEFINITIONS
#----------------------------------------------------------------------------

phelp='''


    Simulation of Neandertal-like-site frequencies under gravel model
        with Neandertal admixture.

        filename: xFNAMEx
        author: Aaron J. Sams
        date: 2016-07-07

    This program runs ms repeatedly for a single segregating site with
        the given input parameters and generates a derived allele count
        distribution for Neandertal-like sites for three samples
        (European, East Asian, Neandertal {1 ind}).

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

subs = ('xTOTSAMPx','xAFRSAMPx','xEURSAMPx','xASNSAMPx','xNEASAMPx',
            'xAFRN0x','xEURN0x','xASNN0x','xNEANx','xAFRG1x',
            'xEURG1x','xASNG1x','xM1Ax','xM1Bx','xM1Cx',
            'xAFRN1x','xEURG2x','xASNG2x','xEURN2x','xASNN2x',
            'xRI2Ex','xRI2Ax','xM2x','xRI1x','xAFRN2x',
            'xTG1x','xTG2x','xTIE2Ex','xTIE2Ax','xTIS2Ex',
            'xTIS2Ax','xTS1x','xTIE1x','xTIS1x','xTS2x',
            'xTG3x','xTS3x')

modp1 = '''xTOTSAMPx 1 -s 1 \
-I 4 xAFRSAMPx xEURSAMPx xASNSAMPx xNEASAMPx 0 \
-n 1 xAFRN0x \
-n 2 xEURN0x \
-n 3 xASNN0x \
-n 4 xNEANx \
-eg 0 1 xAFRG1x \
-eg 0 2 xEURG1x \
-eg 0 3 xASNG1x \
-em 0 1 2 xM1Ax \
-em 0 2 1 xM1Ax \
-em 0 1 3 xM1Bx \
-em 0 3 1 xM1Bx \
-em 0 2 3 xM1Cx \
-em 0 3 2 xM1Cx \
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
'''

modp2all = '''xTOTSAMPx 1 -s 1 \
-I 4 xAFRSAMPx xEURSAMPx xASNSAMPx xNEASAMPx 0 \
-n 1 xAFRN0x \
-n 2 xEURN0x \
-n 3 xASNN0x \
-n 4 xNEANx \
-eg 0 1 xAFRG1x \
-eg 0 2 xEURG1x \
-eg 0 3 xASNG1x \
-em 0 1 2 xM1Ax \
-em 0 2 1 xM1Ax \
-em 0 1 3 xM1Bx \
-em 0 3 1 xM1Bx \
-em 0 2 3 xM1Cx \
-em 0 3 2 xM1Cx \
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
'''

modp2eur = '''xTOTSAMPx 1 -s 1 \
-I 4 xAFRSAMPx xEURSAMPx xASNSAMPx xNEASAMPx 0 \
-n 1 xAFRN0x \
-n 2 xEURN0x \
-n 3 xASNN0x \
-n 4 xNEANx \
-eg 0 1 xAFRG1x \
-eg 0 2 xEURG1x \
-eg 0 3 xASNG1x \
-em 0 1 2 xM1Ax \
-em 0 2 1 xM1Ax \
-em 0 1 3 xM1Bx \
-em 0 3 1 xM1Bx \
-em 0 2 3 xM1Cx \
-em 0 3 2 xM1Cx \
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
'''

modp2asn = '''xTOTSAMPx 1 -s 1 \
-I 4 xAFRSAMPx xEURSAMPx xASNSAMPx xNEASAMPx 0 \
-n 1 xAFRN0x \
-n 2 xEURN0x \
-n 3 xASNN0x \
-n 4 xNEANx \
-eg 0 1 xAFRG1x \
-eg 0 2 xEURG1x \
-eg 0 3 xASNG1x \
-em 0 1 2 xM1Ax \
-em 0 2 1 xM1Ax \
-em 0 1 3 xM1Bx \
-em 0 3 1 xM1Bx \
-em 0 2 3 xM1Cx \
-em 0 3 2 xM1Cx \
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
'''

#----------------------------------------------------------------------------
# FUNCTIONS
#----------------------------------------------------------------------------

###
# Generate ms command line
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
    input = modstring(model)
    sublist = ((subs[i],params[i]) for i in range(len(subs)))
    for i,(sub,param) in enumerate(sublist):
        input = input.replace(sub,str(param))
    return input

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description=phelp)
    p.add_argument('range', type = int, nargs = 2,
                    help='First and last index (inclusive) to use in \
                    out and log files.')
    p.add_argument('--outfile','-o', type=argparse.FileType('w'),
                    default=sys.stdout, help='Output file name. \
                    Defaults to stdout.')
    p.add_argument('--logfile','-l', type=argparse.FileType('w'),
                    default=sys.stderr, help='Log file name.  \
                    Defaults to stderr.')
    p.add_argument('--afrsamp', type=int, default=216,
                    help='Haploid sample size for African sample.')
    p.add_argument('--eursamp', type=int, default=250,
                    help='Haploid sample size for European sample.')
    p.add_argument('--asnsamp', type=int, default=250,
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
    args = p.parse_args()
    return args

def generate_ms_input(args):
    '''Uses input arguments to generate the input for ms.'''
    ###
    # Scaling parameters
    ###
    Nscale = float(7310)                          # Initial Ne
    unit = 4*Nscale                               # Scaling unit
    yearspergen = 25.0
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
    inputa = '/'.join(os.getcwd().split('/')[:-1])+'/bin/ms '
    if args.local:
        inputa = 'ms '

    params = (totsamp, args.afrsamp, args.eursamp, args.asnsamp, args.neasamp,
              afrN0s, eurN0s, asnN0s, neaNs, afrG1,
              eurG1, asnG1, m1as, m1bs, m1cs,
              afrN1s, eurG2, asnG2, eurN2s, asnN2s,
              ri2e, ri2a, m2, ri1, afrN2s,
              tg1s, tg2s, tie2es, tie2as, tis2es,
              tis2as, ts1s, tie1s, tis1s, ts2s,
              tg3s, ts3s)

    inputb = substitute(args.model, subs, params)
    input = inputa+inputb
    return input


###
# Run ms command and capture output
###

def run_ms(input):
    '''Run ms 'input' in the shell and gather input.'''
    args = shsplit(input)
    process = Popen(args, stdout = PIPE, stderr = PIPE)
    return process.communicate()

def parse_ms_output(stdout,afrsz,eursz,asnsz,neasz):
    '''Parse the ms stdout.  Needs:
        'stdout': string from subprocess.
        'afrsz': African sample size.
        'eursz': European sample size.
        'asnsz': Asian sample size.
    Returns:
        seed: random seed from ms run.
        afr, eur, asn, nea: counts of derived alleles per pop.
    '''
    outit = iter(stdout.split('\n'))
    _ = outit.next()
    seed = outit.next()
    afr = 0
    eur = 0
    asn = 0
    nea = 0
    for i in range(4):
        _ = outit.next()
    for i in range(afrsz):
        afr += int(outit.next())
    for i in range(eursz):
        eur += int(outit.next())
    for i in range(asnsz):
        asn += int(outit.next())
    for i in range(neasz):
        nea += int(outit.next())
    return seed, afr, eur, asn, nea

###
# Generate multiple ms runs.
###

def is_nls(afr, eur, asn, nea):
    '''Returns True if data reflect a Neandertal-like
    site.
    afr, eur, asn, nea are derived allele counts.
    '''
    if not nea:
        return False
    else:
        if afr:
            return False
        else:
            if (eur or asn):
                return True
            else:
                return False

def grab_nls(args):
    '''Run ms until a Neandertal-like site (NLS)
    is generated. Returns:
        input: ms command line input
        seed: ms random seed
        (eur, asn, nea): population derived allele counts
        tries: Number of non-NLS simulations needed to
            generate NLS.
    '''
    tries = 0
    while True:
        input = generate_ms_input(args)
        stdout,stderr = run_ms(input)
        seed, afr, eur, asn, nea = parse_ms_output(stdout,
                                                   args.afrsamp,
                                                   args.eursamp,
                                                   args.asnsamp,
                                                   args.neasamp)
        test = is_nls(afr,eur,asn,nea)
        if test:
            return input, seed, [eur, asn, nea], tries
        else:
            tries += 1

#----------------------------------------------------------------------------
# RUN
#----------------------------------------------------------------------------

if __name__=="__main__":
    args = get_args()
    print >> args.outfile, '\t'.join(['#index','EUR','ASN','NEA'])
    print >> args.logfile, '\t'.join(['#index','ms-seed','tries','ms-input'])
    for i in range(args.range[0],args.range[1]+1):
        input, seed, popcounts, tries = grab_nls(args)
        print >> args.outfile, '\t'.join([str(j) for j in [i]+popcounts])
        print >> args.logfile, '\t'.join([str(i), seed, str(tries), input])
    os.remove('seedms')
    args.outfile.close()
    args.logfile.close()

#----------------------------------------------------------------------------
# NOTES
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
