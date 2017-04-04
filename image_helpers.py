import rawpy    #for reading the .nef files
from os.path import join, isfile    #for checking to make sure the file exists in read_image()
import matplotlib.pyplot as plt     #plotting
import numpy as np  #always handy
from mpl_toolkits.axes_grid1 import make_axes_locatable #also plotting

def read_image(filename, color=None):
    """ This function reads in .nef image files and returns a 2D numpy array for a given color.
        Required packages: os.path.join; rawpy
    
        Args:
            filename(str): Filename to be read in. It needs to be a .nef raw image file. If 
                            '.nef' is left off of filename string, it will be 
                            automatically appended.
            color(str or int; optional): Returns RGB color provided. Default is None,
                            which returns all colors. Allowed strings
                            are 'rR', 'g', 'b', 'red', 'green', and 'blue' (not case
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
    #### Just checking the filetype and appending .nef if left off
    if filename[-4::] != '.nef':
        if '.' in filename:
            raise Exception('Invalid file format. {0} is not a .nef file!'.format(filename))
        else:
            filename = filename + '.nef'

    if not isfile(filename):
        raise Exception('{0} does not exist!'.format(filename))

    #### This is the rawpy read call for the .nef files. The following options for postprocess will
    #### ensure that the resulting numpy array is a linear map from the raw sensor data from the
    #### camera. For more information see: https://pythonhosted.org/rawpy/api/rawpy.Params.html#rawpy.Params
    #### for all the parameters that are passed to the postprocess function.
    image = rawpy.imread(filename).postprocess(demosaic_algorithm=rawpy.DemosaicAlgorithm.LINEAR,
                                               output_color=rawpy.ColorSpace.raw,
                                               output_bps=16,
                                               no_auto_bright=True,
                                               adjust_maximum_thr=0.,
                                               gamma=(1, 1))
    #### This whole section is just to check for the color arg and make sure it is appropriate. cix is the color
    #### index for the 3rd dimension of the numpy array returned from the rawpy call above. 0=R, 1=G, and 2=B.
    if color is None:
        return image[:, :, :].astype('float64')
    elif type(color) is str:
        color = color.lower()
        allowedcolors = ['rR', 'g', 'b',
                         'red', 'green', 'blue']
        if color not in allowedcolors:
            raise Exception('{0} is not a valid color option. Pick from {1} (not case sensitive)'.format(color, allowedcolors))
        else:
            if color in ['rR', 'red']:
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
    #### Finally, the returned image, by adding the astype('float64') here, we're making this a float numpy array.
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
                                are 'rR', 'g', 'b', 'red', 'green', and 'blue' (not case
                                sensitive). Allowed ints are 0=Red, 1=Green, and 2=Blue.
            Returns:
                2D numpy.ndarray of image with background subtracted
                unless color is None, which returns 3D numpy.ndarray with last dimension being RGB.
    """
    ring_data = read_image(filename, color=color)
    bg_data = read_image(bgname, color=color)

    return ring_data - bg_data

def quick_plot(image, color=None, block=False):
    """ This function reads in a numpy.ndarray as an image and plots it with imshow quickly. If the numpy.ndarray
        is of dimension 3, the color argument may be used to select which color to plot.
        Dependency: matplotlib.pyplot as plt
        
        Args:
            image(np.ndarray): image data to be plotted, may be 2D or 3D with shape (x,y,3)
            color(str or int; optional): if image is a 3D array, color will select which color to plot.
                        Default is None, which will prompt the user to input a color or choose 'all', which
                        will plot all three colors as subplots. Allowed strings are
                        'rR', 'g', 'b', 'red', 'green', and 'blue' (not case sensitive). Allowed ints are
                        0=Red, 1=Green, and 2=Blue.
            block(bool; optional): Default is False, which will allow further input after this function finishes
                        without needing to close plot figure window. If executing in bash, however, you will want
                        to set block to True so that plot will show. (Note: in order to keep inputting commands 
                        without closing plot figure window, add "&" to the end of bash call.)
            
        Raises:
            ValueError: if image has greater than 3 dimensions
            ValueError: if 3rd dimension has shape greater than 3 (more colors than RGB)
            Exception: if color is string type and not in allowed strings (see above Args description)
            Exception: if color is int type and not 0, 1, or 2
    """
    #### This whole block is just to determine if we need to index the input array for a color and if so which color
    three_color = False
    if len(image.shape) > 3:
        raise ValueError("Too many dimensions! help!")
    elif len(image.shape) == 3 and color is None:
        if image.shape[2] != 3:
            raise ValueError("3D dimension has more colors than RGB!")
        prompt = True
        while prompt:
            c = raw_input("Which color do you want to see (R, G, B, all)?: ")
            if c.lower() not in ['rR', 'g', 'b', 'all']:
                print 'not a valid option, try again'
            else:
                prompt = False
        if c.lower() == 'rR':
            a = image[:, :, 0]
        if c.lower() == 'g':
            a = image[:, :, 1]
        if c.lower() == 'b':
            a = image[:, :, 2]
        if c.lower() == 'all':
            a = image
            three_color = True
    elif len(image.shape) == 3 and color is not None:
        if image.shape[2] != 3:
            raise ValueError("3D dimension has more colors than RGB!")
        if type(color) is str:
            color = color.lower()
            allowedcolors = ['rR', 'g', 'b',
                             'red', 'green', 'blue']
            if color not in allowedcolors:
                raise Exception(
                    '{0} is not a valid color option. Pick from {1} (not case sensitive)'.format(color, allowedcolors))
            else:
                if color in ['rR', 'red']:
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
        a = image[:, :, cix]
    else:
        a = image

    #### Plotting Time
    #### if three_color plots all three colors together
    #### else plots just a single color
    if three_color:
        f, axs = plt.subplots(1, 3, figsize=(18, 6))    #might want to tinker with figsize for your monitor
        names = ['Red', 'Green', 'Blue']
        im = []
        for i in range(3):
            im.append(axs[i].imshow(a[:, :, i], cmap='gray', origin='lower', vmin=a.min(), vmax=a.max()))
            axs[i].set_aspect('equal')
            axs[i].set_title(names[i])
        divider = make_axes_locatable(axs[1])
        f.colorbar(im[2], ax=np.transpose(axs).ravel().tolist(), orientation='horizontal', extend='both', aspect=50, fraction=0.06)
        plt.show(block=block)
    else:
        (x, y) = a.shape
        f, ax = plt.subplots(figsize=(8, float(x)/float(y) * 8.))   #this makes the figsize the same aspect as the image
        im = ax.imshow(a, cmap='gray', origin='lower')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.2)
        f.colorbar(im, extend='both', cax=cax)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.show(block=block)

if __name__ == "__main__":
    folder = "Images"
    shot_number = 15676
    fname = join(folder, "{0:07d}_000.nef".format(shot_number))
    bg_fname = join(folder, "{0:07d}_001.nef".format(shot_number))
    #fp_data = get_image_data(fname, bg_fname, color='b')
    fp_data = read_image(fname, color='b')
    quick_plot(fp_data, block=True)
