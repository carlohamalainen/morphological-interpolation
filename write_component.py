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

from morph import select_structure

input_file = sys.argv[1]
prefix     = sys.argv[2]

data = volumeFromFile(input_file).data

# cs = sorted(uniq(data.flatten()))
# assert cs[0] == 0
# cs = cs[1:]

# d = select_structure(data, cs[1])
d = data

for i in range(d.shape[1]):
    print i
    write_slice(d[:, i, :], '%s_%04d.png' % (prefix, i,))
