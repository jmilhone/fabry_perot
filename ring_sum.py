import numpy as np
import matplotlib.pyplot as plt
from image_helpers import get_image_data
from os.path import join
stepsize = 2
folder = "Images"
shot_number = 15676
fname = join(folder, "{0:07d}_000.nef".format(shot_number))
bg_fname = join(folder, "{0:07d}_001.nef".format(shot_number))
print "getting image"
data = get_image_data(fname, bg_fname)
print "done"

ny, nx = data.shape

x0, y0 = (3066.44, 2036.44)
# x0, y0 = (3000., 2100.0)

x = np.arange(1, nx+1, 1)
y = np.arange(1, ny+1, 1)

xx, yy = np.meshgrid(x-x0, y-y0)
print xx


xmin = xx[0, 0]
xmax = xx[0, -1]
ymin = yy[0, 0]
ymax = yy[-1, 0]

rmax = np.min(np.abs([xmin, xmax, ymin, ymax]))
print rmax
ri = int(np.floor(rmax))

R = np.sqrt(xx**2 + yy**2)
binsize = 0.1

imax = int((rmax ** 2 - 2 * rmax - 1) / (1 + 2 * rmax) / binsize)
binarr = np.fromfunction(lambda i: np.sqrt(2. * (i + 1.) * rmax * binsize + (i + 1.) * binsize ** 2.), (imax,),
                         dtype='float64')

xi0 = int(round(x0)) # np.abs(x - x0).argmin()
yi0 = int(round(y0)) # np.abs(y - y0).argmin()

# chop data matrix into a quadrant for now
# ULdata = data[0:yi0+1, 0:xi0+1]
# URdata = data[0:yi0+1, xi0:]
# BLdata = data[yi0:, 0:xi0+1]
# BRdata = data[yi0:, xi0:]

ULdata = data[yi0-ri:yi0+1, xi0-ri:xi0+1]
URdata = data[yi0-ri:yi0+1, xi0:xi0+ri+1]
BLdata = data[yi0:yi0+ri+1, xi0-ri:xi0+1]
BRdata = data[yi0:yi0+ri+1, xi0:xi0+ri+1]
print "UL", ULdata.shape
print "UR", URdata.shape
print "BL", BLdata.shape
print "BR", BRdata.shape

ULsigarr, _ = np.histogram(R[yi0-ri:yi0+1, xi0-ri:xi0+1], bins=np.concatenate((np.array([0.]), binarr)), weights=ULdata)
URsigarr, _ = np.histogram(R[yi0-ri:yi0+1, xi0:xi0+ri+1], bins=np.concatenate((np.array([0.]), binarr)), weights=URdata)
BLsigarr, _ = np.histogram(R[yi0:yi0+ri+1, xi0-ri:xi0+1], bins=np.concatenate((np.array([0.]), binarr)), weights=BLdata)
BRsigarr, _ = np.histogram(R[yi0:yi0+ri+1, xi0:xi0+ri+1], bins=np.concatenate((np.array([0.]), binarr)), weights=BRdata)

# print ULsigarr.shape
# print URsigarr.shape
# print BLsigarr.shape
# print BRsigarr.shape
thres = 0.3 * np.max(ULsigarr + URsigarr)
i = np.where(ULsigarr + URsigarr > thres)[0]
ni = len(i)

# j = np.arange(-25, 26, 1)
# sarr = stepsize * j
# UB = np.zeros(len(j))
# RL = np.zeros(len(j))
ns = 25
sarr = stepsize * np.arange(-ns, ns+1, 1)
UB = np.zeros(len(sarr))
RL = np.zeros(len(sarr))
# for idx, ix in enumerate(j):
for ix in range(ns):
    UB[ns+ix] = np.sum(((ULsigarr[i-ix]+URsigarr[i-ix]) - (BLsigarr[i+ix]+BRsigarr[i+ix]))**2) / (ni-ix)
    RL[ns+ix] = np.sum(((URsigarr[i-ix]+BRsigarr[i-ix]) - (ULsigarr[i+ix]+BLsigarr[i+ix]))**2) / (ni-ix)
    UB[ns-ix] = np.sum(((ULsigarr[i+ix]+URsigarr[i+ix]) - (BLsigarr[i-ix]+BRsigarr[i-ix]))**2) / (ni-ix)
    RL[ns-ix] = np.sum(((URsigarr[i+ix]+BRsigarr[i+ix]) - (ULsigarr[i-ix]+BLsigarr[i-ix]))**2) / (ni-ix)
plt.plot(sarr, UB, 'or', ms=1)
plt.plot(sarr, RL, 'og', ms=1)
# # ufit = np.polyfit(uB, j, )
plt.show(block=False)

# fig, ax = plt.subplots()
# ax.imshow(ULdata, cmap='gray')
# ax.set_aspect(1.0)
# plt.plot(x0, y0, '.r', ms=1)
# plt.show(block=False)
# #
# fig, ax = plt.subplots()
# ax.plot(binarr, ULsigarr, label='UL', lw=1)
# ax.plot(binarr, URsigarr, label='UR', lw=1)
# ax.plot(binarr, BLsigarr, label='BL', lw=1)
# ax.plot(binarr, BRsigarr, label='BR', lw=1)
# plt.legend(loc='upper left')


# fig, ax = plt.subplots(2,2)
# ax[0][0].imshow(ULdata)
# ax[0][1].imshow(URdata)
# ax[1][0].imshow(BLdata)
# ax[1][1].imshow(BRdata)


# for axis in ax:
#     for a in axis:
#         a.set_aspect(1.0)
plt.show()


