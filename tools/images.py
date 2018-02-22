import rawpy
import exifread
import numpy as np
from os.path import isfile

def check_nef(filename):
    if filename[-4:].lower() != '.nef':
        if isfile(filename+'.NEF'):
            filename += '.NEF'
        elif isfile(filename+'.nef'):
            filename += '.nef'
        else:
            raise Exception('{0} does not exist!'.format(filename))
    return filename

def get_data(filename, color=None):
    filename = check_nef(filename)
    image = rawpy.imread(filename).postprocess(demosaic_algorithm=rawpy.DemosaicAlgorithm.LINEAR,
            output_color=rawpy.ColorSpace.raw, output_bps=16, no_auto_bright=True, 
            adjust_maximum_thr=0., gamma=(1, 1)).astype('float64')

    if color is None:
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
