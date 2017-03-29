import rawpy
from os.path import join, isfile
import matplotlib.pyplot as plt


def read_image(filename, color=None):
    """ This function reads in .nef image files and returns a 2D numpy array for a given color.
        Required packages: os.path.join; rawpy
    
        Args:
            filename(str): Filename to be read in. It needs to be a .nef raw image file. If 
                            '.nef' is left off of filename string, it will be 
                            automatically appended.
            color(str or int; optional): Returns RGB color provided. Default is None,
                            which returns all colors. Allowed strings
                            are 'r', 'g', 'b', 'red', 'green', and 'blue' (not case 
                            sensitive). Allowed ints are 0=Red, 1=Green, and 2=Blue.
        Returns:
            2D numpy.ndarray with indices corresponding to pixels from nef image
            unless color is None, which returns 3D numpy.ndarray with last dimension being RGB.
        Raises:
            Exception: if 'filename' has a format other than .nef
            Exception: if 'filename' does not exist
            Exception: if color is string type and not in allowed strings (see above Args description)
            Exception: if color is int type and not 0, 1, or 2
            
                
        (Dependency Note: To install rawpy package, make sure libraw-dev is installed for linux 
                            systems (using apt), then run "pip install rawpy". You might also 
                            have to upgrade libgcc as well, "conda install libgcc".)
    """
    if filename[-4::] != '.nef':
        if '.' in filename:
            raise Exception('Invalid file format. {0} is not a .nef file!'.format(filename))
        else:
            filename = filename + '.nef'

    if not isfile(filename):
        raise Exception('{0} does not exist!'.format(filename))

    image = rawpy.imread(filename).postprocess(demosaic_algorithm=rawpy.DemosaicAlgorithm.LINEAR,
                                               output_color=rawpy.ColorSpace.raw,
                                               output_bps=16,
                                               no_auto_bright=True,
                                               adjust_maximum_thr=0.,
                                               gamma=(1, 1))

    if color is None:
        return image[:, :, :].astype('float64')
    elif type(color) is str:
        color = color.lower()
        allowedcolors = ['r', 'g', 'b',
                         'red', 'green', 'blue']
        if color not in allowedcolors:
            raise Exception('{0} is not a valid color option. Pick from {1} (not case sensitive)'.format(color, allowedcolors))
        else:
            if color in ['r', 'red']:
                cix = 0
            if color in ['g', 'green']:
                cix = 1
            if color in ['b', 'blue']:
                cix = 2
    else:
        if color not in [0, 1, 2]:
            raise Exception('{0} is not a valid RGB[0,1,2] index, try again.'.format(int(color)))
        else:
            cix = color
    return image[:, :, cix].astype('float64')


def get_image_data(filename, bgname, color=None):
    """ This function reads in .nef image files and returns a 2D numpy array for a given color
        with the 'bgname' image subtracted from the 'filename' image.
            Dependency: read_image

            Args:
                filename(str): Image filename to be read in. It needs to be a .nef raw image file. If 
                                '.nef' is left off of filename string, it will be 
                                automatically appended.
                bgname(str): Background image filename. Same requirements as filename.
                color(str or int; optional): Returns RGB color provided. Default is None,
                                which returns all colors. Allowed strings
                                are 'r', 'g', 'b', 'red', 'green', and 'blue' (not case
                                sensitive). Allowed ints are 0=Red, 1=Green, and 2=Blue.
            Returns:
                2D numpy.ndarray of image with background subtracted
                unless color is None, which returns 3D numpy.ndarray with last dimension being RGB.
    """
    ring_data = read_image(filename, color=color)
    bg_data = read_image(bgname, color=color)

    return ring_data - bg_data



if __name__ == "__main__":
    folder = "Images"
    shot_number = 15676
    fname = join(folder, "{0:07d}_000.nef".format(shot_number))
    bg_fname = join(folder, "{0:07d}_001.nef".format(shot_number))
    print fname
    print bg_fname
    #fp_data = get_image_data(fname, bg_fname, color='b')
    fp_data = read_image(fname, color='b')

    fig, ax = plt.subplots()
    i = ax.imshow(fp_data, cmap='gray')
    fig.colorbar(i)
    ax.set_aspect(1.0)
    plt.show()
