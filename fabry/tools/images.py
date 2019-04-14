from __future__ import print_function, division
import rawpy
import exifread
import numpy as np
import os.path as path
import collections
import skimage.io as io
from . import file_io
import matplotlib.pyplot as plt


image_readers = []

def retrieve_color_section(image, color):
    if color is None or len(image.shape) == 2:
        a = image
    elif color.lower() in ['r', 'red']:
        a = image[:, :, 0]
    elif color.lower() in ['g', 'green']:
        a = image[:, :, 1]
    elif color.lower() in ['b', 'blue']:
        a = image[:, :, 2]
    else:
        raise ValueError('not a valid color choice')
 
    return a


def register_reader(reader_func):
    """Decorator for registering functions as image readers

    Args:
        reader_func (collections.Callable): callable function for reading in image files
    """
    image_readers.append(reader_func)


def check_nef(filename):
    if filename[-4:].lower() != '.nef':
        if path.isfile(filename + '.NEF'):
            filename += '.NEF'
        elif path.isfile(filename + '.nef'):
            filename += '.nef'
        else:
            raise Exception('{0} does not exist!'.format(filename))
    return filename


@register_reader
def read_tiff(fname, **kwargs):
    """Reads .tif file

    Args:
        fname (str): filename to read
        image_idx (int): image idx to read from tiff stack
        return_mean (bool): returns the mean over the tiff stack, overrides image_idx
    Returns:
        np.ndarray: 2d image data, 3d if stack, first dimension being the stack
    """
    if path.splitext(fname)[-1].lower() == '.tif':
        image = io.imread(fname, plugin='tifffile')
        if len(image.shape) == 3:
            idx = kwargs.get('image_index', None)
            return_mean = kwargs.get('return_mean', None)

            if return_mean:
                image = np.mean(image, axis=0)
            elif idx is not None:
                image = image[idx, : , :]

            #print('temporarily averaging over the stack')
            #sigma = np.std(image, axis=0)
            #image = np.mean(image, axis=0)
            #print(sigma.shape)
            #fig, ax = plt.subplots()
            #im = ax.imshow(sigma)# / image * 100)
            #fig.colorbar(im)
            #plt.show()
            #print('temporarily taking the first image in the stack')
            #image = image[0, :, :]
            #print('temporarily taking the third image in the stack')
            #image = image[2, :, :]
        return image
    return None


@register_reader
def read_nef(fname, **kwargs):
    """Reads .nef image files

    Args:
        fname (str): filename to read

    Returns:
        np.ndarray: 2d image data
    """
    color = kwargs.get('color', None)
    if path.splitext(fname)[-1].lower() == '.nef':
        image = rawpy.imread(fname)
        data = image.postprocess(demosaic_algorithm=rawpy.DemosaicAlgorithm.LINEAR,
                                 output_color=rawpy.ColorSpace.raw, output_bps=16, no_auto_bright=True,
                                 adjust_maximum_thr=0., gamma=(1, 1)).astype('float64')
        if color is not None:
            data = retrieve_color_section(data, color)
        return data
    else:
        return None


@register_reader
def read_npy(fname, **kwargs):
    """Reads numpy binary files

    Args:
        fname (str): filename to read
    
    Returns: 
        np.ndarray: 2d image data
    """
    if path.splitext(fname)[-1].lower() == ".npy":
        data = np.load(fname)
        color = kwargs.get('color', None)
        if color is not None:
            data = retrieve_color_section(data, color)
        return data
    else:
        return None


@register_reader
def read_h5(fname, **kwargs):
    if path.splitext(fname)[-1].lower() in [".h5", ".hdf5"]:
        data = file_io.h5_2_dict(fname)
        data = data.get('image', None)
        data.astype(np.float64)
        color = kwargs.get('color', None)
        if color is not None:
            data = retrieve_color_section(data, color)
        return data
    else:
        return None


def get_data(filename, color=None, image_index=None, return_mean=False):
    """Reads image data from filename

    Args:
        filename (str): filename to read
        color (Union[int, str]): [0,2] for rgb, or a letter from rgb

    Returns: 
        np.ndarray: 2d image data
    """
    for reader in image_readers:
        image = reader(filename, color=color, image_index=image_index, 
                return_mean=return_mean)
        if image is not None:
            break
    else:
        raise ValueError('Image not found or not supported')

    if path.splitext(filename)[-1].lower() == '.tif':
        # Ignore color propery for tif files
        return image

    return image


def get_metadata(filename):
    """Returns metadata from image filename

    Args:
        filename (str): file to get metadata from

    Returns:
        dict: dictionary containing metadata
    """
    filename = check_nef(filename)
    with open(filename, 'rb') as f:
        tags = exifread.process_file(f, details=False)
    date = str(tags['Image DateTime'].values.replace(':', '_').split(' ')[0])
    time = str(tags['Image DateTime'].values.split(' ')[1])
    model = str(tags['Image Model'].values)
    return {'date': date, 'time': time, 'q': model}
