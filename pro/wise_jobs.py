#!/usr/bin/env python

import sys, pyfits, os, subprocess

NWORKERS = int(sys.argv[1])
ALLSKY = True

if ALLSKY:
  indfile = os.environ['WISE_DATA'] + '/index-allsky-L1b.fits'
else:
  indfile = os.environ['WISE_DATA'] + '/index-metadata-L1b.fits'
  
mjdhdu = pyfits.open(indfile)

mjd = (mjdhdu[1].data)['MJD']
dt = (max(mjd) - min(mjd))/NWORKERS

ngroup = NWORKERS/4

basename = 'wisejobs'

outpath = sys.argv[2]

tpad = 0.000001

for i in range(0, ngroup):
    thisname = basename + str(i)
    outfile = open(thisname, 'w')
    outlines = []
    for j in range(0, 4):
        ind = 4*i+j
        mjdmin = str(mjd[0] + ind*dt - tpad) 
        mjdmax = str(mjd[0] + (ind+1)*dt + tpad) #pad to ensure no files left behind
        line = 'echo cpu, tpool_nthreads=2 \& wise_clean_batch, [' + mjdmin + ', ' + mjdmax + '], outpath = \\"' + outpath  + '\\"' + ALLSKY*', /allsky' + ' |idl &> b' + str(ind) + '.log &\n'
        outlines.append(line)
    outlines.append('./wise_continue.pl\n')
    outfile.writelines(outlines)
    outfile.close()
    p = subprocess.Popen('chmod 755 ' + thisname, shell=True)
    p.communicate()
