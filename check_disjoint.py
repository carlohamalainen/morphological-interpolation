"""

Morphological image interpolation based on the paper
"A Morphology-Based Approach for Interslice Interpolation of Anatomical
Slices From Volumetric Images", Albu et al., Biomedical Engineering
(Volume:55 ,  Issue: 8 ), 2008.

http://dx.doi.org/10.1109/TBME.2008.921158

In the code below, ABL refers to this paper.

Copyright (c) 2015, Carlo Hamalainen
The Centre for Advanced Imaging
The University of Queensland
All rights reserved.

See LICENSE for full license details.

"""


import math
import numpy as np
import os
import sys

import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe

from nipype.utils.filemanip import split_filename
from nipype.interfaces.utility import Merge, Rename, IdentityInterface, Select, Function

from nipype.interfaces.base import (
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    StdOutCommandLineInputSpec,
    StdOutCommandLine,
    File,
    Directory,
    InputMultiPath,
    OutputMultiPath,
    traits,
    isdefined,
    BaseInterface,
    BaseInterfaceInputSpec,
)

from pyminc.volumes.factory import volumeFromFile, volumeFromDescription
import skimage.morphology as morph
import scipy.ndimage.measurements

from utils import *

from morph import select_structure, overlap

import os.path

import multiprocessing as mp

def check_overlap((f_i, f_j,)):
    data_i = volumeFromFile(f_i).data
    data_j = volumeFromFile(f_j).data

    data_i[np.where(data_i > 0)] = 1
    data_j[np.where(data_j > 0)] = 1

    data_i = data_i.astype('uint8')
    data_j = data_j.astype('uint8')

    for k in range(data_i.shape[1]):
        if overlap(data_i[:, k, :], data_j[:, k, :]):
            print (f_i, f_j)
            return (f_i, f_j)

    print f_i, f_j, None

    return None


# Files that we have:

files = []
for i in range(95):
    # OUTPUT OF INDIVIDUAL INTERPOLATIONS; POSSIBLY NOT DISJOINT:
    # filename_i = '/mnt/home/carlo/raw_interpolated_components/component_%d.mnc' % i

    # FOR TESTING, CHECKING THE OUTPUT OF THE ALLEGEDLY DISJOINT COMPONENTS:
    filename_i = '/scratch/intersection_output/component_%d.mnc' % i

    if os.path.exists(filename_i):
        files.append(filename_i)


# Pairs of files to check:

file_pairs = []

for (i, f_i) in enumerate(files):
    for (j, f_j) in enumerate(files):
        if j <= i: continue
        print f_i, f_j
        file_pairs.append((f_i, f_j))

pool = mp.Pool(processes=8)
results = pool.map(check_overlap, file_pairs)

for r in results:
    if r is None: continue
    print 'overlap', r[0], r[1]


"""
for i in range(95):
    filename_i = '/?????/morph_00_init/raw_interpolated_components/component_%d.mnc' % i

    if not os.path.exists(filename_i): continue

    data_i = volumeFromFile(filename_i).data.astype('uint8')

    for j in range(i + 1, 95):
        filename_j = '/?????/morph_00_init/raw_interpolated_components/component_%d.mnc' % j

        if not os.path.exists(filename_j): continue

        data_j = volumeFromFile(filename_j).data.astype('uint8')

        if overlap(data_i, data_j):
            o = 'OVERLAP! :('
        else:
            o = 'disjoint, yay :)'

        print i, j, o
"""
