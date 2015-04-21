from morph import *

def write_component_as_png(component_nr):
    input_file = '/mnt/home/carlo/00-init-label-vol.mnc'
    input_file = os.path.abspath(input_file)

    volume = volumeFromFile(input_file)
    data = volume.data

    cs = sorted(uniq(data.flatten()))
    assert cs[0] == 0
    cs = cs[1:]

    v = volume
    d = v.data

    d = select_structure(d, cs[component_nr], tol=0.5)

    for j in range(d.shape[1]):
        if np.max(d[:, j, :]) > 0:
            print component_nr, j
            write_slice(d[:, j, :], '/mnt/home/carlo/pngs/component_%04d_%04d.png' % (component_nr, j,))

if __name__ == '__main__': write_component_as_png(int(sys.argv[1]))
