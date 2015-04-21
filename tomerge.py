if False:
    from multiprocessing import Pool
    pool = Pool(processes=4)

    component_nrs = [x for x in range(0, 29 + 1) if x != 8]
    print pool.map(do_component, component_nrs)

if False:
    for component_a in range(0, 29 + 1):
        if component_a == 8: continue

        vol_a  = volumeFromFile('/scratch/final_smoothed_disjoint_components/component_%d.mnc' % component_a)

        for component_b in range(component_a + 1, 29 + 1):
            if component_b == 8: continue
            vol_b  = volumeFromFile('/scratch/final_smoothed_disjoint_components/component_%d.mnc' % component_b)

            print component_a, component_b

            for j in range(vol_a.data.shape[1]):
                assert not overlap(vol_a.data[:, j, :].astype('uint8'), vol_b.data[:, j, :].astype('uint8'))


if False:
    data = volumeFromFile('data/small.mnc').data
    cs = sorted(uniq(data.flatten()))
    assert cs[0] == 0
    cs = cs[1:]
    del data

    new_volume = volumeLikeFile('/scratch/final_smoothed_disjoint_components/component_0.mnc', '/scratch/final_smoothed_disjoint_components/final.mnc', dtype='float')
    new_volume.data[:] = 0

    for component_nr in range(0, 29 + 1):
        if component_nr == 8: continue

        print component_nr

        component = volumeFromFile('/scratch/final_smoothed_disjoint_components/component_%d.mnc' % component_nr)

        new_volume.data[:] += cs[component_nr]*component.data.astype('uint8')

    # FIXME Hardcoded...
    os.system('rm -fr /export/nif02/uqchamal/scratch_morph')
    os.system('mkdir /export/nif02/uqchamal/scratch_morph')
    workflow.base_dir = '/export/nif02/uqchamal/scratch_morph'

    new_volume.writeFile()
    new_volume.closeVolume()




if False:
    failed = [ 0, 2, 13, 45, 80, 94 ]

    input_file = '/mnt/home/carlo/00-init-label-vol.mnc'
    input_file = os.path.abspath(input_file)

    volume = volumeFromFile(input_file)
    data = volume.data

    volume.closeVolume()

    cs = sorted(uniq(data.flatten()))
    assert cs[0] == 0
    cs = cs[1:]

    for x in failed:
        print 'failed component: index %d, value %g' % (x, cs[x],)


# FINAL AS FLOAT
if True:
    data = volumeFromFile('/mnt/home/carlo/00-init-label-vol.mnc').data
    cs = sorted(uniq(data.flatten()))
    assert cs[0] == 0
    cs = cs[1:]
    del data

    # We skip the last one because it was dodgy.
    cs = cs[:-1]
    assert len(cs) == 94

    os.system('rm -fv /export/nif02/uqchamal/init_00_output/final_float.mnc')

    v = volumeFromFile('/scratch/intersection_output/component_0.mnc')

    new_volume = volumeFromDescription('/export/nif02/uqchamal/init_00_output/final_float.mnc',
                                       v.dimnames,
                                       v.sizes[:3],
                                       v.starts,
                                       v.separations,
                                       volumeType='float',
                                       dtype='float')


    new_volume.data[:] = 0

    del v

    for component_nr in range(len(cs)):
        # if component_nr == 8: continue

        print '%d\t%.4f' % (component_nr, cs[component_nr],)

        component = volumeFromFile('/scratch/intersection_output/component_%d.mnc' % component_nr)

        x = np.zeros(component.data.shape, dtype='float')
        x[component.data != 0] = 1

        new_volume.data[:] += cs[component_nr]*x

    # FIXME Hardcoded...
    # os.system('rm -fr /export/nif02/uqchamal/scratch_morph')
    # os.system('mkdir /export/nif02/uqchamal/scratch_morph')
    # workflow.base_dir = '/export/nif02/uqchamal/scratch_morph'

    new_volume.writeFile()
    new_volume.closeVolume()




# if __name__ == '__main__': go(int(sys.argv[1]))


# data_8 = np.load('data_8.npz')['arr_0']
# interp_3_to_4 = blap(data_8[:, 3, :], data_8[:, 4, :], morph.disk(radius=1))

# overlapping_pairs = [ (1, 2), (1, 18), (2, 10), (2, 18), (3, 16), (3, 18), (3, 23), (4, 15), (4, 17), (4, 29), (5, 6), (5, 7), (5, 22), (7, 22), (9, 10), (9, 12), (9, 14), (9, 21), (9, 22), (10, 12), (10, 14), (10, 16), (10, 18), (10, 19), (10, 20), (10, 21), (10, 22), (10, 23), (10, 25), (10, 26), (10, 28), (11, 12), (12, 15), (12, 19), (12, 21), (12, 29), (13, 18), (14, 16), (14, 19), (14, 22), (14, 23), (15, 21), (15, 29), (16, 18), (16, 20), (16, 21), (16, 23), (16, 26), (18, 20), (18, 21), (19, 22), (19, 26), (20, 21), (21, 29), (22, 26), (23, 24), (23, 25), (23, 26), (23, 28), (24, 25), (24, 26), (25, 27), (25, 28), (27, 28) ] + [(10, 15)]

if False:
    import multiprocessing as mp

    files = [ '/home/carlo/Desktop/component_0.mnc', '/home/carlo/Desktop/component_10.mnc', '/home/carlo/Desktop/component_11.mnc', '/home/carlo/Desktop/component_12.mnc', '/home/carlo/Desktop/component_13.mnc', '/home/carlo/Desktop/component_14.mnc', '/home/carlo/Desktop/component_15.mnc', '/home/carlo/Desktop/component_16.mnc', '/home/carlo/Desktop/component_17.mnc', '/home/carlo/Desktop/component_18.mnc', '/home/carlo/Desktop/component_19.mnc', '/home/carlo/Desktop/component_1.mnc', '/home/carlo/Desktop/component_20.mnc', '/home/carlo/Desktop/component_21.mnc', '/home/carlo/Desktop/component_22.mnc', '/home/carlo/Desktop/component_23.mnc', '/home/carlo/Desktop/component_24.mnc', '/home/carlo/Desktop/component_25.mnc', '/home/carlo/Desktop/component_26.mnc', '/home/carlo/Desktop/component_27.mnc', '/home/carlo/Desktop/component_28.mnc', '/home/carlo/Desktop/component_29.mnc', '/home/carlo/Desktop/component_2.mnc', '/home/carlo/Desktop/component_3.mnc', '/home/carlo/Desktop/component_4.mnc', '/home/carlo/Desktop/component_5.mnc', '/home/carlo/Desktop/component_6.mnc', '/home/carlo/Desktop/component_7.mnc', '/home/carlo/Desktop/component_9.mnc']

    overlaps = {}

    def check_overlap((f_i, f_j,)):
        data_i = volumeFromFile(f_i).data
        data_j = volumeFromFile(f_j).data

        data_i[np.where(data_i > 0)] = 1
        data_j[np.where(data_j > 0)] = 1

        data_i = data_i.astype('uint8')
        data_j = data_j.astype('uint8')

        for k in range(data_i.shape[1]):
            if overlap(data_i[:, k, :], data_j[:, k, :]):
                #assert i < j
                return (i, j)

        return None

    file_pairs = []

    for (i, f_i) in enumerate(files):
        for (j, f_j) in enumerate(files):
            if j <= i: continue
            print f_i, f_j
            file_pairs.append((f_i, f_j))

    pool = mp.Pool(processes=4)
    results = pool.map(check_overlap, file_pairs)

    print(results)

if False:
    overlapping_pairs = [(0, 11), (0, 17), (0, 38), (0, 48), (0, 64), (1, 2), (1, 6), (1, 8), (1, 11), (1, 18), (1, 22), (1, 47), (1, 74), (1, 78), (2, 6), (2, 8), (2, 18), (2, 19), (2, 22), (2, 39), (2, 47), (2, 65), (2, 69), (2, 73), (2, 74), (2, 81), (2, 84), (2, 85), (2, 86), (3, 8), (3, 22), (3, 55), (4, 11), (4, 45), (4, 70), (5, 7), (5, 11), (5, 14), (5, 23), (5, 24), (5, 26), (5, 42), (5, 55), (5, 56), (6, 18), (6, 22), (6, 69), (6, 74), (6, 84), (6, 85), (7, 12), (7, 20), (7, 26), (7, 38), (7, 42), (7, 43), (7, 44), (8, 18), (8, 22), (8, 23), (8, 24), (8, 26), (8, 28), (8, 30), (8, 47), (8, 55), (8, 74), (9, 18), (9, 19), (9, 21), (9, 30), (9, 47), (9, 61), (9, 62), (9, 63), (9, 65), (9, 67), (9, 76), (9, 89), (10, 45), (10, 48), (11, 17), (11, 23), (11, 30), (11, 35), (11, 42), (11, 45), (11, 62), (11, 63), (11, 64), (11, 70), (11, 74), (12, 13), (12, 20), (12, 22), (12, 26), (12, 27), (12, 28), (12, 29), (12, 32), (12, 34), (12, 43), (12, 44), (12, 50), (12, 60), (12, 69), (12, 80), (13, 20), (13, 26), (13, 27), (13, 28), (13, 52), (13, 55), (13, 69), (13, 80), (13, 87), (13, 88), (14, 26), (14, 55), (14, 56), (15, 31), (15, 74), (15, 78), (16, 23), (16, 24), (16, 25), (16, 30), (16, 35), (16, 45), (16, 49), (16, 62), (16, 64), (16, 68), (16, 71), (17, 23), (17, 24), (17, 25), (17, 28), (17, 35), (18, 19), (18, 47), (18, 61), (19, 61), (19, 65), (19, 67), (19, 81), (20, 34), (20, 43), (20, 50), (20, 51), (20, 52), (20, 60), (20, 80), (21, 30), (21, 62), (21, 64), (22, 23), (22, 27), (22, 28), (22, 55), (22, 69), (22, 74), (22, 85), (23, 24), (23, 25), (23, 28), (23, 30), (23, 35), (23, 38), (23, 58), (23, 74), (24, 25), (24, 28), (24, 30), (24, 49), (25, 35), (25, 49), (26, 27), (26, 28), (26, 32), (26, 42), (26, 44), (26, 55), (26, 56), (26, 80), (27, 28), (27, 43), (27, 55), (27, 69), (28, 55), (28, 74), (29, 34), (29, 57), (29, 60), (30, 35), (30, 47), (30, 49), (30, 61), (30, 62), (30, 63), (30, 64), (30, 74), (31, 83), (32, 33), (32, 34), (32, 56), (32, 72), (32, 80), (33, 34), (33, 53), (33, 56), (33, 72), (34, 53), (34, 57), (34, 60), (34, 77), (34, 80), (34, 82), (34, 88), (35, 38), (35, 62), (35, 64), (35, 74), (36, 53), (36, 66), (36, 71), (36, 77), (36, 80), (36, 87), (36, 88), (36, 89), (36, 91), (36, 92), (37, 45), (37, 59), (38, 42), (38, 44), (39, 40), (39, 46), (39, 54), (39, 65), (39, 81), (39, 86), (39, 91), (41, 47), (41, 55), (41, 85), (42, 44), (43, 51), (43, 60), (43, 75), (43, 82), (45, 48), (45, 49), (45, 58), (45, 59), (45, 70), (46, 65), (46, 86), (46, 90), (46, 91), (47, 49), (47, 55), (47, 61), (47, 63), (47, 74), (47, 85), (48, 58), (48, 70), (49, 55), (49, 77), (49, 89), (50, 51), (50, 52), (50, 69), (51, 52), (51, 60), (51, 82), (52, 60), (52, 77), (52, 80), (52, 82), (52, 92), (53, 55), (53, 56), (53, 57), (53, 66), (53, 77), (53, 88), (53, 91), (55, 56), (55, 66), (55, 69), (55, 77), (55, 80), (55, 85), (55, 87), (55, 88), (55, 89), (56, 66), (56, 72), (56, 88), (57, 60), (57, 77), (57, 82), (57, 91), (58, 64), (60, 77), (60, 82), (60, 91), (61, 63), (61, 65), (61, 74), (62, 63), (62, 64), (63, 65), (63, 74), (65, 67), (65, 71), (65, 81), (65, 90), (65, 91), (66, 71), (66, 77), (67, 81), (68, 71), (69, 73), (69, 78), (69, 82), (69, 84), (69, 85), (69, 86), (69, 87), (69, 88), (69, 91), (71, 76), (71, 77), (71, 81), (71, 85), (71, 86), (71, 87), (71, 88), (71, 89), (71, 91), (71, 92), (73, 74), (73, 91), (74, 78), (77, 88), (77, 89), (77, 91), (77, 92), (79, 91), (79, 93), (80, 87), (80, 88), (80, 92), (81, 85), (81, 86), (82, 91), (82, 92), (82, 93), (84, 85), (84, 87), (85, 86), (85, 87), (85, 88), (85, 91), (86, 88), (86, 91), (87, 88), (87, 91), (87, 92), (88, 89), (88, 91), (91, 92), (91, 93), (92, 93),]

    scores = {}

    for (feature_i, feature_j) in overlapping_pairs:
        # print (feature_i, feature_j)

        assert feature_i != feature_j

        if (feature_i, feature_j) not in scores: scores[(feature_i, feature_j)] = 0
        if (feature_j, feature_i) not in scores: scores[(feature_j, feature_i)] = 0

        data_i = volumeFromFile('/mnt/home/carlo/components_tmp_for_removing_overlaps/component_%d.mnc' % feature_i).data
        data_j = volumeFromFile('/mnt/home/carlo/components_tmp_for_removing_overlaps/component_%d.mnc' % feature_j).data

        # FIXME Why the !!!?! are the values in component_?.mnc not 0/1?????
        data_i[np.where(data_i > 0)] = 1
        data_j[np.where(data_j > 0)] = 1

        data_i = data_i.astype('uint8')
        data_j = data_j.astype('uint8')

        for k in range(data_i.shape[1]):
            if overlap(data_i[:, k, :], data_j[:, k, :]):

                slice_intersection = np.bitwise_and(data_i[:, k, :], data_j[:, k, :])

                # if np.max(slice_intersection.flatten()) == 0: continue

                # print feature_i, feature_j, k

                # slice_i_part = (+1)*(np.bitwise_and(data_i[:, k, :], slice_intersection).astype('float'))
                # slice_j_part = (+9)*(np.bitwise_and(data_j[:, k, :], slice_intersection).astype('float'))

                # save_this = slice_i_part + slice_j_part

                # write_slice(save_this, 'save_this_%d_%d_%04d.png' % (feature_i, feature_j, k))

                # save_this = 10*data_i.astype('float')[:, k, :] - 2*data_j.astype('float')[:, k, :]

                # write_slice(save_this, 'save_this_%d_%d_%04d.png' % (feature_i, feature_j, k))

                # assert False

                # Winner is largest feature?
                feature_i_size = float(data_i[:, k, :].flatten().sum())
                feature_j_size = float(data_j[:, k, :].flatten().sum())

                # Metrics.
                intersection_size = float(slice_intersection.flatten().sum())

                #if intersection_size/feature_i_size > 0.8:
                #    # feature_i is almost completely contained in the intersection, so it's
                #    # probably one of those dodgy extreme dilation situations.
                #    scores[(feature_i, feature_j)] += 1000 # FIXME Make this tuneable parameter.

                print 'feature %d vs feature %d, slice %d' % (feature_i, feature_j, k)
                print '\tintersection_size/feature_i_size: %.2f' % (intersection_size/feature_i_size)
                print '\tintersection_size/feature_j_size: %.2f' % (intersection_size/feature_j_size)
                print

                if intersection_size/feature_i_size > 0.5:
                    scores[(feature_i, feature_j)] += 100

                if intersection_size/feature_j_size > 0.5:
                    scores[(feature_j, feature_i)] += 100


                scores[(feature_i, feature_j)] += intersection_size/feature_i_size
                scores[(feature_j, feature_i)] += intersection_size/feature_j_size

                # if score[(i, j)] is big, we prefer to nuke i when j is overlapping.
                # Need to check if (i, j) and (j, i) in scores, take maximum.

                print '\tintersection_size/feature_j_size: %.2f' % (intersection_size/feature_j_size)


if False:

    data_i = volumeFromFile('/home/carlo/Desktop/component_27.mnc').data.astype('uint8')
    data_j = volumeFromFile('/home/carlo/Desktop/component_28.mnc').data.astype('uint8')

    intersection = np.bitwise_and(data_i, data_j)

    for k in range(data_i.shape[1]):
        if overlap(data_i[:, k, :], data_j[:, k, :]):
            print k
            write_slice(0*data_i[:, k, :] + data_j[:, k, :], 'metrics_%04d.png' % (k,))

if False:

    data_i = volumeFromFile('/scratch/morph/morpho_1/merged_minc_sink/merged_minc_file/merged.mnc').data.astype('uint8')
    data_j = volumeFromFile('/scratch/morph/morpho_2/merged_minc_sink/merged_minc_file/merged.mnc').data.astype('uint8')

    intersection = np.bitwise_and(data_i, data_j)

    for k in range(data_i.shape[1]):
        if overlap(data_i[:, k, :], data_j[:, k, :]):
            print k
            write_slice(np.bitwise_and(data_i[:, k, :], data_j[:, k, :]), 'overlap_%04d.png' % (k,))

def highest_score(scores):
    max_score = None
    max_idx   = None

    for ((i, j), s) in scores.iteritems():
        if max_score is None:
            max_score = s
            max_idx   = (i, j)

        else:
            if s > max_score:
                max_score = s
                max_idx   = (i, j)

    return max_idx, max_score

if False:
    os.system('rm -fr /scratch/intersection_input/*')
    # Make 0-1 versions of each file.
    for i in set(reduce(operator.add, [[x[0], x[1]] for x in overlapping_pairs], [])):
        print '0-1:', i

        v_source = volumeFromFile('/home/carlo/Desktop/component_%d.mnc' % i)

        v = volumeLikeFile('/home/carlo/Desktop/component_%d.mnc' % i, '/scratch/intersection_input/input_%d.mnc' % i, dtype='ubyte')
        v.data[:] = v_source.data[:]

        v.data[np.where(v.data > 0)] = 1
        v.writeFile()
        v.closeVolume()

        v_source.closeVolume()

if False:

    scores = pickle.load(open('scores_00_init.pkl', 'r'))
    print scores

    #os.system('rm -fr /scratch/intersection_output')
    #os.system('mkdir /scratch/intersection_output')
    #os.system('cp -v /mnt/home/carlo/raw_interpolated_components/* /scratch/intersection_output/')

    print highest_score(scores)

    while len(scores) > 0:
        (i, j), s = highest_score(scores)

        assert scores[(j, i)] < s

        v_i = volumeFromFile('/scratch/intersection_output/component_%d.mnc' % i, readonly=False)
        v_j = volumeFromFile('/scratch/intersection_output/component_%d.mnc' % j)

        # v_i loses to v_j.
        print 'removing intersection(%d, %d) from component %d' % (i, j, i,)

        """
        for k in range(v_i.data.shape[1]):
            if overlap(v_i.data[:, k, :].astype('uint8'), v_j.data[:, k, :].astype('uint8')):
                print k

                v_i.data[:, k, :] -= np.bitwise_and(v_i.data[:, k, :].astype('uint8'), v_j.data[:, k, :].astype('uint8'))
        """

        v_i.data[(v_i.data > 0) * (v_j.data > 0)] = 0

        v_i.volumeType = 'ubyte' # FIXME What???????
        v_i.writeFile()
        v_i.closeVolume()
        v_j.closeVolume()

        del scores[(i, j)]
        del scores[(j, i)]


if False:
    files = [ 'input_10.mnc', 'input_11.mnc', 'input_12.mnc', 'input_13.mnc', 'input_14.mnc', 'input_15.mnc', 'input_16.mnc', 'input_17.mnc', 'input_18.mnc', 'input_19.mnc', 'input_1.mnc', 'input_20.mnc', 'input_21.mnc', 'input_22.mnc', 'input_23.mnc', 'input_24.mnc', 'input_25.mnc', 'input_26.mnc', 'input_27.mnc', 'input_28.mnc', 'input_29.mnc', 'input_2.mnc', 'input_3.mnc', 'input_4.mnc', 'input_5.mnc', 'input_6.mnc', 'input_7.mnc', 'input_9.mnc' ]


    for i in range(len(files)):
        for j in range(i + 1, len(files)):
            assert i != j
            assert i < j

            if (i, j) in [(0, 5), (0, 15)]: continue

            v_i = volumeFromFile('/scratch/intersection_output/' + files[i])
            v_j = volumeFromFile('/scratch/intersection_output/' + files[j])

            for k in range(v_i.data.shape[1]):
                print i, j, k
                assert not overlap(v_i.data[:, k, :].astype('uint8'), v_j.data[:, k, :].astype('uint8'))

            v_i.closeVolume()
            v_j.closeVolume()




if False:
    v = volumeFromFile('/home/carlo/input_2.mnc')

    data = v.data[:, 13, :]
    data = data.astype('uint8')

    # Invert the image.
    data = invert(data)

    # Dilate it.
    data = morph.binary_dilation(data, selem=morph.disk(radius=2)).astype('uint8')
    data = morph.binary_dilation(data, selem=morph.disk(radius=2)).astype('uint8')
    data = morph.binary_dilation(data, selem=morph.disk(radius=2)).astype('uint8')

    # Invert back.
    data = invert(data)

    show_slice(data)




