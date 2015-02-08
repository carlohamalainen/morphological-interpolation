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

def nr_components(x):
    """
    Use scipy.ndimage.measurements.label to find the number of
    components in the binary (?) image.

    http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.measurements.label.html
    """

    label, num_features = scipy.ndimage.measurements.label(x)
    return num_features

def components(x):
    """
    Return a list of the components of x.
    """

    label, num_features = scipy.ndimage.measurements.label(x)

    return [(label == i).astype('uint8') for i in range(1, num_features + 1)]

def overlap(x, y):
    """
    Do two slices overlap?

    FIXME This is slow?
    """

    return int(np.bitwise_and(x, y).max())

def is_subset(x, y):
    """
    Is x a subset of y, element-wise?
    """

    return int(np.bitwise_and(x, np.bitwise_xor(x, y)).max()) == 0

def conditional_dilation(B, A, structure_component):
    """
    Conditional dilation of set B with respect to reference set A.

    This is ABL equation (2).
    """

    assert B.dtype == 'uint8'

    B_dilation = morph.binary_dilation(B, selem=structure_component).astype('uint8')

    return np.bitwise_and(B_dilation, A)

def transition_sequence(start_set, ref_set, structure_component):
    """
    Construct a sequence of conditional dilations starting from
    'start_set' and using 'ref_set' as the reference set at each step.
    This process is guaranteed to terminate due to reference Serra [25]
    in ABL.
    """

    dilations = [conditional_dilation(start_set, ref_set, structure_component)]

    while not np.array_equal(dilations[-1], ref_set):
        new_dilation = conditional_dilation(dilations[-1], ref_set, structure_component)
        dilations.append(new_dilation)

    return dilations

def distance(A, B):
    """
    Distance between two images, defined as the cardinality of
    the symmetric difference (Equation 7 in ABL).
    """

    flat = np.bitwise_xor(A, B).flatten()
    result = float(flat.sum())
    del flat

    return result

def interp(A, B, structure_component):
    _, result = _interp(A, B, structure_component)
    return result

def _interp(A, B, structure_component):
    """
    Interpolate between two slices A and B.

    The seq array corresponds to Equation 5 in ABL, and the final
    answer seq[index_med] uses Equation 6.
    """

    assert overlap(A, B)

    # A \intersection B -> A
    # print r'A \intersection B -> A'
    A_and_B_to_A = transition_sequence(np.bitwise_and(A, B), A, structure_component)

    # A \intersection B -> B
    # print r'A \intersection B -> B'
    A_and_B_to_B = transition_sequence(np.bitwise_and(A, B), B, structure_component)

    l_A = len(A_and_B_to_A)
    l_B = len(A_and_B_to_B)

    # print 'l_A:', l_A
    # print 'l_B:', l_B

    seq = []

    for i in range(1, max(l_A, l_B) + 1):
        # print 'interp: i =', i

        if i <= min(l_A, l_B):
            # A_and_B_to_A[l_A - i] \union A_and_B_to_B[i - 1]
            result = np.bitwise_or(A_and_B_to_A[l_A - i], A_and_B_to_B[i - 1])
        elif l_A <= i <= l_B:
            # The A sequence ran out but the B sequence is longer, so use that (as
            # we want the end of seq to be most like B).
            result = A_and_B_to_B[i - 1]
        elif l_B <= i <= l_A:
            # The B sequence ran out so fudge it with the opposite end of the A sequence.
            # This is not the same definition as the paper which uses the constant
            # A_and_B_to_A[l_A - 1] instead.
            # result = A_and_B_to_A[l_A - 1]
            result = A_and_B_to_A[l_A - i + 1]
        else:
            assert False, 'Unknown case in interp().'

        assert result.shape == A.shape
        assert result.shape == B.shape
        seq.append(result)

    # Equation 6:
    distances = np.array([np.abs(distance(s, A) - distance(s, B)) for s in seq])
    index_med = np.argmin(distances)
    result = seq[index_med]

    # FIXME Avoid getting stuck?
    # assert distance(result, A) != 0
    # assert distance(result, B) != 0

    return seq, result

def complement_in(x, y):
    """
    Compute the complement of x with respect to y.

    >>> complement_in(np.array([[1, 0], [0, 1]]), np.array([[1, 1], [0, 0]]))
    array([[0, 1],
           [0, 0]])
    """

    return np.bitwise_and(np.logical_not(x), y)

def compute_Ref(S, i, A):
    assert 0 <= i < len(S)

    def drop_ith(x, i):
        """
        Drop the i-th element of list x.
        """

        result = x[:i] + x[i+1:]
        assert len(result) == len(x) - 1
        return result

    return reduce(np.bitwise_or, [complement_in(s, A) for s in drop_ith(S, i)])

def split(A, Bs, structure_component):
    for B in Bs: assert A.shape == B.shape

    S   = [np.bitwise_and(B, A) for B in Bs]
    Ref = [compute_Ref(S, i, A) for i in range(len(S))]

    for s in S: assert is_subset(s, A)

    p = 0

    last_few = []
    max_last_few_len = 10

    not_finished = True

    while not_finished:
        new_S   = [conditional_dilation(S[i], Ref[i], structure_component) for i in range(len(S))]
        new_Ref = [compute_Ref(new_S, i, A)                                for i in range(len(S))]

        for s in new_S: assert is_subset(s, A)

        _d = distance(A, reduce(np.bitwise_or, S))

        # Detect a loop. Not mentioned in ABL.
        if len(last_few) < max_last_few_len:
            # Not enough data yet...
            last_few.append(_d)
        else:
            last_few.append(_d)
            last_few.pop(0)

            if len(uniq(last_few)) <= 2:
                not_finished = False

        # Detect perfect convergence:
        if _d == 0: not_finished = False

        S   = new_S
        Ref = new_Ref
        p  += 1

    # Make sure the regions are disjoint.
    for i in range(len(S)):
        for j in range(len(S)):
            if i < j:
                intersection = np.bitwise_and(S[i], S[j])
                S[i] -= intersection
                S[j] -= intersection
                assert S[i].min() == 0
                assert S[j].min() == 0
                assert S[i].max() == 1
                assert S[j].max() == 1

    for i in range(len(S)):
        for j in range(len(S)):
            if i < j:
                assert np.bitwise_and(S[i], S[j]).max() == 0

    return S


def multi_component_dilation(A, B, structure_component):
    """
    Equation 9 in APL...
    Dilate each of the components in the multi-blob thing.
    """

    print 'multi_component_dilation...'
    assert nr_components(A) == 1

    for b in B:
        assert b.shape == A.shape
        assert overlap(A, b)

    dilations = [morph.binary_dilation(b, selem=structure_component) for b in B]
    union_dilations = reduce(np.bitwise_or, dilations)

    # Now check if A \subseteq union_dilations...

    if is_subset(A, union_dilations):
        print 'multi_component_dilation: using interp(A, union_dilations)'
        result = interp(A, union_dilations)
    else:
        print 'multi_component_dilation: splitting A and doing %d one-to-one interps' % len(B)

        # FIXME should nr_components(A_split) > 1?
        A_split = split(A, B, structure_component)

        assert len(A_split) == len(B)

        result = reduce(np.bitwise_or, [interp(A_split[i], B[i]) for i in range(len(B))])

        # The thing that we produce has to be different to B (otherwise we'll get stuck).
        assert distance(result, reduce(np.bitwise_or, B)) != 0

    return result

def interp_one_to_one(A, B, structure_component):
    print 'interp_one_to_one...'
    assert A.shape == B.shape

    result = interp(A, B, structure_component)
    assert result.shape == A.shape

    return result

def interp_one_to_many(A, B, structure_component):
    print 'interp_one_to_many...'
    assert A.shape == B.shape

    assert nr_components(A) == 1
    assert nr_components(B)  > 1

    Bs = components(B)
    for b in Bs: assert b.shape == A.shape

    result = multi_component_dilation(A, Bs, structure_component)
    assert result.shape == A.shape

    return result

def interp_many_to_many(A, B, structure_component):
    print 'interp_many_to_many'

    assert A.shape == B.shape

    comp_A = components(A)
    comp_B = components(B)

    targets = {}

    for i in range(len(comp_A)):
        targets[i] = [b for b in comp_B if overlap(comp_A[i], b)]

    results = []

    for (i, Bs) in targets.iteritems():
        results.append(blap(comp_A[i], reduce(np.bitwise_or, Bs), structure_component))

    return reduce(np.bitwise_or, results)

def interp_zero_to_one(A, B, structure_component):
    print 'interp_zero_to_one'
    assert A.shape == B.shape

    assert nr_components(A) == 0
    assert nr_components(B) == 1

    components_B = components(B)

    centroids = [find_centroid(b) for b in components_B]

    results = []
    for i in range(len(centroids)):
        new = np.copy(A).astype('uint8')
        new[tuple(centroids[i])] = 1 # single pixel extreme point

        # FIXME Equation 12 is used here. We can't use interp() because that assumes that the two
        # regions overlap, but here the centroid might be outside the blob.
        results.append(extreme_dilation(new, components_B[i], structure_component))

    result = reduce(np.bitwise_or, results)

    # Edge case - interpolation doesn't make the image
    # any smaller, so truncate now to zero. Not mentioned in ABL.
    if distance(result, B) == 0: result = np.zeros(result.shape, dtype='uint8')

    return result

def interp_zero_to_zero(A, B, structure_component):
    print 'interp_zero_to_zero'
    assert A.shape == B.shape
    return np.copy(A)

def closest_point(pt, A):
    A_indexes = np.argwhere(A > 0)

    def euclidean_distance(x, y):
        return math.sqrt((x[0] - y[0])**2 + (x[1] - y[1])**2)

    distances = [euclidean_distance(pt, a) for a in A_indexes]

    return A_indexes[np.argmin(distances)]

def extreme_dilation(E, X, structure_component):
    print 'extreme_dilation'
    assert E.shape == X.shape

    assert nr_components(E) == 1
    assert len(np.argwhere(E > 0)) == 1 # single pixel
    assert nr_components(X) == 1

    # As defined in ABL, the conditional dilations won't work
    # if E and X do not overlap. So translate the extreme point.
    if not overlap(E, X):
        original_pt = tuple(np.argwhere(E > 0)[0])
        closest = closest_point(original_pt, X)

        E_translated = np.zeros(E.shape, dtype='uint8')
        E_translated[tuple(closest)] = 1

        assert overlap(E_translated, X)

        E = E_translated

    assert E.dtype == 'uint8'
    assert X.dtype == 'uint8'

    # This is my own heuristic for the size of the structuring element.
    r = int(0.95*math.sqrt(X.flatten().sum()/math.pi))
    structure_component = morph.disk(radius=min(9, r))

    d = conditional_dilation(E, X, structure_component)

    seq = []
    seq.append(d)

    while distance(d, X) > 0:
        d = conditional_dilation(d, X, structure_component)
        seq.append(d)

    distances = np.array([np.abs(distance(s, E) - distance(s, X)) for s in seq])
    index_med = np.argmin(distances)

    result = seq[index_med]
    del seq
    assert result.shape == E.shape

    return result

def interp_zero_to_many(A, B, structure_component):
    print 'interp_zero_to_many'
    assert A.shape == B.shape

    assert nr_components(A) == 0
    assert nr_components(B)  > 1

    components_B = components(B)

    centroids = [find_centroid(b) for b in components_B]

    results = []
    for i in range(len(centroids)):
        new = np.copy(A).astype('uint8')
        new[tuple(centroids[i])] = 1 # single pixel extreme point

        results.append(extreme_dilation(new, components_B[i], structure_component))

    return reduce(np.bitwise_or, results)

def blap(A, B, structure_component):
    assert A.shape == B.shape

    nr_A = nr_components(A)
    nr_B = nr_components(B)

    # One to one:
    if nr_A == 1 and nr_B == 1:
        result = interp_one_to_one(A, B, structure_component)

    # One to many:
    elif nr_A == 1 and nr_B > 1:
        result = interp_one_to_many(A, B, structure_component)

    # Many to one:
    elif nr_B == 1 and nr_A > 1:
        result = interp_one_to_one(B, A, structure_component)

    # Many to many:
    elif nr_A > 1 and nr_B > 1:
        result = interp_many_to_many(A, B, structure_component)

    # Zero to one:
    elif nr_A == 0 and nr_B == 1:
        result = interp_zero_to_one(A, B, structure_component)

    # Zero to many:
    elif nr_A == 0 and nr_B > 1:
        result = interp_zero_to_many(A, B, structure_component)

    # One to zero:
    elif nr_A == 1 and nr_B == 0:
        result = interp_zero_to_one(B, A, structure_component)

    # Many to zero:
    elif nr_A >  1 and nr_B == 0:
        result = interp_zero_to_many(B, A, structure_component)

    # Zero to zero:
    elif nr_A == 0 and nr_B == 0:
        result = interp_zero_to_zero(A, B, structure_component)

    # Unhandled case:
    else:
        assert False

    return result

def select_structure(x, float_val, tol=0.5):
    """
    Return a binary image consisting of the pixels p in x such
    that float_val - tol < p < float_val + p.
    """

    return np.where((x < float_val + tol) & (x > float_val - tol), 1, 0)

def find_centroid(x):
    """
    Find the centroid of a 2D shape by averaging the Euclidean points.
    http://en.wikipedia.org/wiki/Centroid#Of_a_finite_set_of_points
    """

    assert nr_components(x) == 1

    locs = np.argwhere(x > 0)

    c = locs.mean(axis=0)
    assert c.shape == (2,)

    c = (int(c[0]), int(c[1]))

    assert 0 <= c[0] <= x.shape[0]
    assert 0 <= c[1] <= x.shape[1]

    return c

"""

class InterpolationInputSpec(BaseInterfaceInputSpec):
    volume = File(
                exists=True,
                desc='volume to be interpolated',
                mandatory=True)

    interpolation_dimension = traits.Int(
                                desc='dimension to be interpolated',
                                mandatory=True)

    k = traits.Int(
            desc='number of levels to interpolate. Produces 2^k - 1 new slices',
            mandatory=True)

class InterpolationOutputSpec(TraitedSpec):
    interpolated_volume = File(exists=True, desc='interpolated volume')

class Interpolation(BaseInterface):
    input_spec  = InterpolationInputSpec
    output_spec = InterpolationOutputSpec

    def _run_interface(self, runtime):
        fname = self.inputs.volume

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        fname = self.inputs.volume
        _, base, _ = split_filename(fname)
        outputs["thresholded_volume"] = os.path.abspath(base + '_thresholded.nii')
        return outputs


"""

class NumpySliceToPNGInputSpec(BaseInterfaceInputSpec):
    in_file = File(
                exists=True,
                desc='Numpy slice in npz format',
                mandatory=True)

class NumpySliceToPNGOutputSpec(TraitedSpec):
    out_file = File(desc='rendering of Numpy slice as a PNG', exists=True)

class NumpySliceToPNG(BaseInterface):
    input_spec  = NumpySliceToPNGInputSpec
    output_spec = NumpySliceToPNGOutputSpec

    def _run_interface(self, runtime):
        in_file = self.inputs.in_file

        _, base, _ = split_filename(in_file)

        self.png = base + '.png'

        write_slice(load_pickled_array(in_file), self.png)

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(self.png)
        return outputs

class InterpolateBetweenSlicesInputSpec(BaseInterfaceInputSpec):
    slices = InputMultiPath(
                    traits.File,
                    desc='input slices (only supply two!!!)',
                    exists=True,
                    mandatory=True)

    level = traits.Int(
            desc='level of this interpolation slice',
            mandatory=True)

class InterpolateBetweenSlicesOutputSpec(TraitedSpec):
    interpolated_slice = File(desc='interpolation between slice0 and slice1', exists=True)

class InterpolateBetweenSlices(BaseInterface):
    input_spec  = InterpolateBetweenSlicesInputSpec
    output_spec = InterpolateBetweenSlicesOutputSpec

    def _run_interface(self, runtime):
        assert len(self.inputs.slices) == 2
        slice0_fname = self.inputs.slices[0]
        slice1_fname = self.inputs.slices[1]
        level        = self.inputs.level # FIXME Unused?

        slice0 = load_pickled_array(slice0_fname)
        slice1 = load_pickled_array(slice1_fname)

        _, base0, _ = split_filename(slice0_fname)
        _, base1, _ = split_filename(slice1_fname)

        self.output_file = base0 + '_to_' + base1 + '.npz'

        # FIXME This should be a parameter?
        structure_component = morph.disk(radius=9)

        result = blap(slice0, slice1, structure_component)

        np.savez_compressed(self.output_file, result)
        print '==>', self.output_file

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['interpolated_slice'] = os.path.abspath(self.output_file)
        return outputs

def select_function(x, i): return x[i]

def create_stage(stage_nr, workflow, inputs, inputs_nr_slices, slice_name):
    """
    Don't use this directly, see build_workflow() instead.

    Create an interpolation stage. Mutates the 'workflow' argument.
    """

    # Selectors into 'inputs'
    select_inputs = {}
    for i in range(inputs_nr_slices):
        fi = Function(input_names=['x', 'i'], output_names=['out_file'], function=select_function)

        select_inputs[i] = pe.Node(interface=fi, name='select_inputs_%s_%i' % (slice_name, i,))
        select_inputs[i].inputs.i = i
        workflow.connect(inputs, 'out_files', select_inputs[i], 'x')

    # Interpolations.
    interp_nodes = []
    for i in range(inputs_nr_slices - 1):
        interp_node = pe.Node(interface=InterpolateBetweenSlices(), name='interp_%s_%d_%08d' % (slice_name, stage_nr, i,))

        select_node = pe.Node(interface=Select(index=[i, i + 1]), name='select_%s_%d_%d' % (slice_name, i, i + 1,))
        workflow.connect(inputs, 'out_files', select_node, 'inlist')
        workflow.connect(select_node, 'out', interp_node, 'slices')
        interp_node.inputs.level = stage_nr

        interp_nodes.append(interp_node)

    # Rename slices.
    renamers = []
    k = 0
    rename = pe.Node(interface=Rename(), name='rename_%s_%d_%08d' % (slice_name, stage_nr, k,))
    rename.inputs.format_string = 'slice_%08d.npz' % k
    workflow.connect(select_inputs[0], 'out_file', rename, 'in_file')
    renamers.append(rename)
    k += 1

    for i in range(len(interp_nodes)):
        rename = pe.Node(interface=Rename(), name='rename_%s_%d_%08d' % (slice_name, stage_nr, k,))
        rename.inputs.format_string = 'slice_%08d.npz' % k
        workflow.connect(interp_nodes[i], 'interpolated_slice', rename, 'in_file')
        renamers.append(rename)
        k += 1

        rename = pe.Node(interface=Rename(), name='rename_%s_%d_%08d' % (slice_name, stage_nr, k,))
        rename.inputs.format_string = 'slice_%08d.npz' % k
        workflow.connect(select_inputs[i + 1], 'out_file', rename, 'in_file')
        renamers.append(rename)
        k += 1

    # Could skip this unless we want to see intermediate steps.
    datasink = pe.Node(nio.DataSink(), name='sinker_%s_%d' % (slice_name, stage_nr,))
    for (i, rename) in enumerate(renamers):
        workflow.connect(rename, 'out_file', datasink, 'slices.@%d' % i)

    # If we want to do another stage, use the out_file's of renamers.
    return renamers


class SlicesToMincInputSpec(BaseInterfaceInputSpec):
    input_files = InputMultiPath(
                        traits.File,
                        desc='input slices',
                        exists=True,
                        mandatory=True)

    slice_dimension = traits.Int(desc='dimension that changes along slices', mandatory=True)

    dimension_names = traits.List(traits.Str,   desc='dimension names',                mandatory=True)
    sizes           = traits.List(traits.Int,   desc='dimension sizes',                mandatory=True)
    starts          = traits.List(traits.Float, desc='dimension starts',               mandatory=True)
    separations     = traits.List(traits.Float, desc='dimension separations (steps)',  mandatory=True)

class SlicesToMincOutputSpec(TraitedSpec):
    out_file = File(desc='MINC file', exists=True)

class SlicesToMinc(BaseInterface):
    input_spec  = SlicesToMincInputSpec
    output_spec = SlicesToMincOutputSpec

    def _run_interface(self, runtime):
        input_files = self.inputs.input_files

        _, base0, _ = split_filename(input_files[0])
        _, basen, _ = split_filename(input_files[-1])

        self.output_file = base0 + '_' + basen + '.mnc'

        # Clobbers the output file!
        try:    os.remove(self.output_file)
        except: pass

        v = volumeFromDescription(
                    self.output_file,
                    self.inputs.dimension_names,
                    self.inputs.sizes,
                    self.inputs.starts,
                    self.inputs.separations,
                    volumeType='ushort')

        for (i, fname) in enumerate(input_files):
            s = load_pickled_array(fname)

            dim = self.inputs.slice_dimension

            if dim == 0:
                v.data[i, :, :] = s
            elif dim == 1:
                v.data[:, i, :] = s
            elif dim == 2:
                v.data[:, :, i] = s
            else:
                assert False, 'Invalid slice dimension %d' % dim

        if dim == 0:
            assert v.data.shape[0] == len(input_files)
        elif dim == 1:
            assert v.data.shape[1] == len(input_files)
        elif dim == 2:
            assert v.data.shape[2] == len(input_files)

        v.writeFile()
        v.closeVolume()

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(self.output_file)
        return outputs

class ExtractFeatureInputSpec(BaseInterfaceInputSpec):
    in_file = File(
                exists=True,
                desc='MINC file',
                mandatory=True)

    float_val = traits.Float(desc='float value of feature', mandatory=True)
    tol       = traits.Float(desc='tolerance; we find pixels p such that float_val - tol < p < float_val + p.', default_value=0.5)

    dim_to_interpolate = traits.Int(desc='dimension to interpolate', mandatory=True)

class ExtractFeatureOutputSpec(TraitedSpec):
    out_files = InputMultiPath(
                        traits.File,
                        desc='slices in npz format',
                        exists=True,
                        mandatory=True)

class ExtractFeature(BaseInterface):
    input_spec  = ExtractFeatureInputSpec
    output_spec = ExtractFeatureOutputSpec

    def _run_interface(self, runtime):
        in_file             = self.inputs.in_file
        float_val           = self.inputs.float_val
        tol                 = self.inputs.tol
        dim_to_interpolate  = self.inputs.dim_to_interpolate

        self.out_files = []
        volume = volumeFromFile(in_file)
        data = volume.data

        data = select_structure(data, float_val, tol)

        _, base, _ = split_filename(in_file)

        assert dim_to_interpolate in [0, 1, 2]
        assert len(data.shape) == 3

        if dim_to_interpolate == 0:
            slices = [data[i, :, :] for i in range(data.shape[0])]
        elif dim_to_interpolate == 1:
            slices = [data[:, i, :] for i in range(data.shape[1])]
        elif dim_to_interpolate == 2:
            slices = [data[:, :, i] for i in range(data.shape[2])]

        for (i, s) in enumerate(slices):
            fname = 'slice_%08d.npz' % i
            np.savez_compressed(fname, s)
            self.out_files.append(os.path.abspath(fname))

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_files'] = self.out_files
        return outputs


class MergeMincsInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiPath(
                        traits.File,
                        desc='MINC files to merge',
                        exists=True,
                        mandatory=True)

    labels = traits.List(traits.Float, desc='labels for each MINC file', mandatory=True)

    dimension_names = traits.List(traits.Str,   desc='dimension names',                mandatory=True)
    sizes           = traits.List(traits.Int,   desc='dimension sizes',                mandatory=True)
    starts          = traits.List(traits.Float, desc='dimension starts',               mandatory=True)
    separations     = traits.List(traits.Float, desc='dimension separations (steps)',  mandatory=True)

class MergeMincsOutputSpec(TraitedSpec):
    out_file = File(desc='Merged MINC file', exists=True)

class MergeMincs(BaseInterface):
    input_spec  = MergeMincsInputSpec
    output_spec = MergeMincsOutputSpec

    def _run_interface(self, runtime):
        in_files = self.inputs.in_files
        labels   = self.inputs.labels

        self.out_file = os.path.abspath('merged.mnc')

        v = volumeFromDescription(
                    self.out_file,
                    self.inputs.dimension_names,
                    self.inputs.sizes,
                    self.inputs.starts,
                    self.inputs.separations,
                    volumeType='ushort')

        components = {}

        for (i, lab) in enumerate(labels):
            components[i] = volumeFromFile(in_files[i]).data.astype('uint8')

            v.data += lab*components[i]

        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                assert not overlap(components[i], components[j])

        v.writeFile()
        v.closeVolume()

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = self.out_file
        return outputs

def build_workflow(workflow, initial_inputs, slice_name, slice_size, dim_to_interpolate, dimension_names, nr_interpolation_steps, starts, new_separations, new_sizes):

    current_nr_slices = slice_size

    stage_outfiles = initial_inputs
    for k in range(1, nr_interpolation_steps + 1):
        stage_outfiles = create_stage(k, workflow, stage_outfiles, current_nr_slices, slice_name)
        current_nr_slices = calc_new_sizes(current_nr_slices, 1)

    final_outfiles = stage_outfiles

    png_sink = pe.Node(nio.DataSink(), name='pngs_%s' % (slice_name,))

    final_pngs = []

    for (i, node) in enumerate(final_outfiles):
        png = pe.Node(
                    interface=NumpySliceToPNG(),
                    name='final_npz_to_png_%s_%d' % (slice_name, i,))
        final_pngs.append(png)

        workflow.connect(node, 'out_file', png, 'in_file')
        workflow.connect(png, 'out_file', png_sink, 'pngs.@%d' % i)

    # Merge final npz slices into a single thing.
    merge_final_npz = pe.Node(interface=Merge(len(final_outfiles)), name='merge_final_npz_%s' % (slice_name,))

    for (i, f) in enumerate(final_outfiles):
        workflow.connect(final_outfiles[i], 'out_file', merge_final_npz, 'in%d' % (i + 1))

    slices_to_minc = pe.Node(
                            interface=SlicesToMinc(
                                            slice_dimension=dim_to_interpolate,
                                            dimension_names=dimension_names,
                                            sizes=new_sizes,
                                            starts=starts,
                                            separations=new_separations),
                            name='slices_to_minc_%s' % (slice_name,))

    workflow.connect(merge_final_npz, 'out', slices_to_minc, 'input_files')

    minc_sink = pe.Node(interface=nio.DataSink(), name='minc_%s' % (slice_name,))
    workflow.connect(slices_to_minc, 'out_file', minc_sink, 'minc_file')

    return minc_sink

def calc_new_sizes(old_size, nr_levels):
    """
    Each interpolation stage produces a new slice between each original
    slice. This function calculates the total number of slices after
    nr_levels stages.

    >>> calc_new_sizes(4, 2)
    13

    >>> calc_new_sizes(10, 3)
    73
    """

    new_sizes = old_size

    for k in range(nr_levels):
        new_sizes += new_sizes - 1

    return new_sizes

def go():
    # Top-level parameters:
    input_file = 'data/small.mnc'
    dim_to_interpolate = 1
    nr_interpolation_steps = 1
    workflow_name = 's62'

    ################

    input_file = os.path.abspath(input_file)

    volume = volumeFromFile(input_file)
    data = volume.data

    assert len(data.shape) == 3
    assert dim_to_interpolate in [0, 1, 2]
    assert nr_interpolation_steps >= 1

    old_separations = volume.separations

    new_separations = old_separations[:]
    new_separations[dim_to_interpolate] /= 2**nr_interpolation_steps

    new_sizes = volume.sizes[:3]
    for s in new_sizes: assert s > 0
    new_sizes[dim_to_interpolate] = calc_new_sizes(new_sizes[dim_to_interpolate], nr_interpolation_steps)

    workflow = pe.Workflow(name=workflow_name)

    cs = sorted(uniq(data.flatten()))
    assert cs[0] == 0
    cs = cs[1:]

    cs = cs[21:][:2] # FIXME Just for testing...

    extractors = []

    minc_sink = {}

    for c in cs:
        slice_name = ('feature_%0.4f' % c).replace('.', '_')
        print 'building workflow for:', slice_name

        extract_node = pe.Node(
                            interface=ExtractFeature(
                                            in_file=input_file,
                                            float_val=c,
                                            tol=0.5,
                                            dim_to_interpolate=dim_to_interpolate),
                            name='extract_feature_' + slice_name)

        extractors.append(extract_node)

        slice_size = data.shape[dim_to_interpolate]

        minc_sink[slice_name] = build_workflow(
                                        workflow, extract_node, slice_name,
                                        slice_size, dim_to_interpolate, volume.dimnames,
                                        nr_interpolation_steps, volume.starts, new_separations, new_sizes)

    merge_minc_names = pe.Node(interface=Merge(len(minc_sink)), name='merge_minc_names')

    for (i, minc) in enumerate(minc_sink.itervalues()):
        workflow.connect(minc, 'out_file', merge_minc_names, 'in%d' % (i + 1))

    actual_final_merge = pe.Node(interface=MergeMincs(), name='actual_final_merge_mincs')

    workflow.connect(merge_minc_names, 'out', actual_final_merge, 'in_files')

    actual_final_merge.inputs.labels = map(float, cs)

    actual_final_merge.inputs.dimension_names = volume.dimnames
    actual_final_merge.inputs.sizes           = new_sizes
    actual_final_merge.inputs.starts          = volume.starts
    actual_final_merge.inputs.separations     = new_separations

    merged_minc_sink = pe.Node(interface=nio.DataSink(), name='merged_minc_sink')
    workflow.connect(actual_final_merge, 'out_file', merged_minc_sink, 'merged_minc_file')

    # workflow.run()
    workflow.run(plugin='MultiProc', plugin_args={'n_procs' : 4})




