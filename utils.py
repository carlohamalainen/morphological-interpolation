import math
import numpy as np
import os

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

import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt

from pyminc.volumes.factory import volumeFromFile, volumeFromDescription
import skimage.morphology as morph
import scipy.ndimage.measurements

def uniq(X): return list(set(X))

def write_slice(x, f):
    plt.imshow(x, interpolation='none')
    plt.savefig(f)
    plt.close()

def show_slice(x):
    plt.imshow(x, interpolation='none')
    plt.show()

def load_pickled_array(fname):
    x = np.load(fname)
    assert len(x.keys()) == 1
    return x[x.keys()[0]]

def show_npz(fname):
    show_slice(load_pickled_array(fname))

def load_pklz(f):
    import pickle
    import gzip

    return pickle.load(gzip.open(f))
