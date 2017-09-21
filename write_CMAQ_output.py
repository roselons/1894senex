#!/usr/bin/env python

import os, sys, argparse
import numpy as np

sys.path.append('/home/rpernak/revised_python_libraries')
import utils

parser = argparse.ArgumentParser(\
  description='Read in files generated with ' + \
  'CMAQ_grid_point_collect.py and write output to a single text file.')
parser.add_argument('--indir', type=str, default='CMAQ_npz', \
  help='Directory with .npz files from CMAQ_grid_point_collect.py.')
parser.add_argument('--outfile', type=str, \
  default='CMAQ_SENEX_jun2013_Surface.txt', \
  help='Path to text file to which the output is written.')
args = parser.parse_args()

inDir = args.indir; utils.file_check(inDir)
nFiles = utils.ls('ls %s/*.npz' % inDir)

outFile = args.outfile
outFP = open(outFile, 'w')
for nFile in nFiles:
  print os.path.basename(nFile)

  nDat = np.load(nFile)
  concArr = nDat['mConc']
  latArr = nDat['mLats']
  lonArr = nDat['mLons']
  pArr = nDat['modelP']
  rvmrArr = nDat['rvmr']

  for conc, rvmr, lat, lon, p in \
    zip(concArr, rvmrArr, latArr, lonArr, pArr):
    outLine = '%10.3f%10.3f%10.3f%10.3f%10.3f\n' % \
      (conc[0], rvmr, lat, lon, p[0])
    outFP.write(outLine)
  # end loop over arrays

# end loop over nFiles

outFP.close()