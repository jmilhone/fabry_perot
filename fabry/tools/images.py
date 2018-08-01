from __future__ import print_function, division
import rawpy
import exifread
import numpy as np
import os.path as path
import collections

image_readers = []


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
def read_nef(fname):
    """Reads .nef image files

    Args:
        fname (str): filename to read

    Returns:
        np.ndarray: 2d image data
    """
    if path.splitext(fname)[-1].lower() == '.nef':
        image = rawpy.imread(fname)
        data = image.postprocess(demosaic_algorithm=rawpy.DemosaicAlgorithm.LINEAR,
                                 output_color=rawpy.ColorSpace.raw, output_bps=16, no_auto_bright=True,
                                 adjust_maximum_thr=0., gamma=(1, 1)).astype('float64')
        return data
    else:
        return None


@register_reader
def read_npy(fname):
    """Reads numpy binary files

    Args:
        fname (str): filename to read
    
    Returns: 
        np.ndarray: 2d image data
    """
    if path.splitext(fname)[-1].lower() == ".npy":
        data = np.load(fname)
        return data
    else:
        return None


def get_data(filename, color=None):
    """Reads image data from filename

    Args:
        filename (str): filename to read
        color (Union[int, str]): [0,2] for rgb, or a letter from rgb

    Returns: 
        np.ndarray: 2d image data
    """
    for reader in image_readers:
        image = reader(filename)
        if image is not None:
            break
    else:
        raise ValueError('Image not found or not supported')

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
    return {'date': date, 'time': time, 'model': model}
