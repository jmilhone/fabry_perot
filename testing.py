import rawpy
import seaborn as sb
import matplotlib.pyplot as plt

image = rawpy.imread('Images/0015676_000.nef')

proc_image = image.postprocess(demosaic_algorithm=rawpy.DemosaicAlgorithm.LINEAR,
                                               output_color=rawpy.ColorSpace.raw,
                                               output_bps=16,
                                               no_auto_bright=True,
                                               adjust_maximum_thr=0.,
                                               gamma=(1,1))
blue_proc = proc_image[:, :, 2]
blue_raw = image.raw_image[image.raw_colors==2]
blue_proc = blue_proc[image.raw_colors==2]

sb.jointplot(blue_proc.flatten(), blue_raw.flatten())
plt.show()