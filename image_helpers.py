import rawpy
from os.path import join
import matplotlib.pyplot as plt


def read_image(filename):
    image = rawpy.imread(filename).postprocess(output_bps=16,
                                               no_auto_bright=True,
                                               adjust_maximum_thr=0.)
    return image[:, :, 2].astype('float64')


def get_image_data(filename, bgname):
    ring_data = read_image(filename)
    bg_data = read_image(bgname)

    return ring_data - bg_data


if __name__ == "__main__":
    folder = "Images"
    shot_number = 15676
    fname = join(folder, "{0:07d}_000.nef".format(shot_number))
    bg_fname = join(folder, "{0:07d}_001.nef".format(shot_number))
    print fname
    print bg_fname
    fp_data = get_image_data(fname, bg_fname)

    fig, ax = plt.subplots()
    ax.imshow(fp_data, cmap='gray')
    ax.set_aspect(1.0)
    plt.show()
