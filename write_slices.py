from morph import *

def write_slices_png(component_nr):
    minc_file = '/mnt/home/carlo/raw_interpolated_components/component_%d.mnc' % (component_nr,)

    v = volumeFromFile(minc_file)

    d = v.data

    for j in range(d.shape[1]):
        if np.max(d[:, j, :]) > 0:
            print component_nr, j
            write_slice(d[:, j, :], '/mnt/home/carlo/component_slices_of_all_final_components/component_%04d_%04d.png' % (component_nr, j,))

if __name__ == '__main__': write_slices_png(int(sys.argv[1]))
