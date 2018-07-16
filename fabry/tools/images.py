import rawpy
import exifread
import numpy as np
import os.path as path

image_readers = []

def register_reader(reader_func):
    """
    Decorator for registering functions as image readers
    """
    image_readers.append(reader_func)

def check_nef(filename):
    if filename[-4:].lower() != '.nef':
        if path.isfile(filename+'.NEF'):
            filename += '.NEF'
        elif path.isfile(filename+'.nef'):
            filename += '.nef'
        else:
            raise Exception('{0} does not exist!'.format(filename))
    return filename

@register_reader
def read_nef(fname):
    """
    Reads .nef image files
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
    """
    Reads numpy binary files
    """
    if path.splitext(fname)[-1].lower() == ".npy":
        data = np.load(fname)
        return data
    else:
        return None

def get_data(filename, color=None):

    for reader in image_readers:
        image = reader(filename)
        if image is not None:
            break
    else:
        raise ValueError('Image not found or not supported')

    if color is None or len(image.shape) == 2:
        a = image
    elif color.lower() in ['r', 'red']:
        a = image[:,:,0]
    elif color.lower() in ['g', 'green']:
        a = image[:,:,1]
    elif color.lower() in ['b', 'blue']:
        a = image[:,:,2]
    else:
        raise ValueError('not a valid color choice')

    return a

def get_metadata(filename):
    filename = check_nef(filename)
    with open(filename, 'rb') as f:
        tags = exifread.process_file(f, details=False)
    date = str(tags['Image DateTime'].values.replace(':','_').split(' ')[0])
    time = str(tags['Image DateTime'].values.split(' ')[1])
    model = str(tags['Image Model'].values)
    return {'date':date, 'time':time, 'model':model}
