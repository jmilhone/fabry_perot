import numpy as np
import matplotlib.pyplot as plt
from image_helpers import get_image_data
from os.path import join
stepsize = 1
folder = "Images"
shot_number = 15676
fname = join(folder, "{0:07d}_000.nef".format(shot_number))
bg_fname = join(folder, "{0:07d}_001.nef".format(shot_number))
print "getting image"
data = get_image_data(fname, bg_fname, color='b')
print "done"
print data.shape

def quick_ringsum(dat, x0, y0, binsize=0.1, quadrants=True):

    ny, nx = dat.shape

    x = np.arange(1, nx+1, 1)
    y = np.arange(1, ny+1, 1)
    xx, yy = np.meshgrid(x-x0, y-y0)
    R = np.sqrt(xx**2 + yy**2)
    xmin = xx[0, 0]
    xmax = xx[0, -1]
    ymin = yy[0, 0]
    ymax = yy[-1, 0]

    ri = int(np.min(np.abs([xmin, xmax, ymin, ymax])) - 1)
    ri = int(.5 * 4020)
    imax = int((ri ** 2. - 2. * ri - 1.) / (1. + 2. * ri) / binsize)
    binarr = np.fromfunction(lambda i: np.sqrt(2. * (i + 1.) * ri * binsize + (i + 1.) * binsize ** 2.), (imax,),
                             dtype='float64')

    xi0 = int(round(x0))  # np.abs(x - x0).argmin()
    yi0 = int(round(y0))  # np.abs(y - y0).argmin()
    print binarr
    if quadrants:
        ULdata = dat[yi0 - ri:yi0 + 1, xi0 - ri:xi0 + 1]
        URdata = dat[yi0 - ri:yi0 + 1, xi0:xi0 + ri + 1]
        BLdata = dat[yi0:yi0 + ri + 1, xi0 - ri:xi0 + 1]
        BRdata = dat[yi0:yi0 + ri + 1, xi0:xi0 + ri + 1]

        RUL = R[yi0-ri:yi0+1, xi0-ri:xi0+1]
        RUR = R[yi0-ri:yi0+1, xi0:xi0+ri+1]
        RBL = R[yi0:yi0+ri+1, xi0-ri:xi0+1]
        RBR = R[yi0:yi0+ri+1, xi0:xi0+ri+1]
        print RUL[-1, -1]
        print RUR[-1, 0]
        print RBL[0, -1]
        print RBR[0, 0]
        ULsigarr, _ = np.histogram(R[yi0-ri:yi0+1, xi0-ri:xi0+1], bins=np.concatenate((np.array([0.]), binarr)), weights=ULdata)
        URsigarr, _ = np.histogram(R[yi0-ri:yi0+1, xi0:xi0+ri+1], bins=np.concatenate((np.array([0.]), binarr)), weights=URdata)
        BLsigarr, _ = np.histogram(R[yi0:yi0+ri+1, xi0-ri:xi0+1], bins=np.concatenate((np.array([0.]), binarr)), weights=BLdata)
        BRsigarr, _ = np.histogram(R[yi0:yi0+ri+1, xi0:xi0+ri+1], bins=np.concatenate((np.array([0.]), binarr)), weights=BRdata)
        return binarr, ULsigarr, URsigarr, BLsigarr, BRsigarr
    else:
        sigarr, _ = np.histogram(R, bins=np.concatenate((np.array([0.]), binarr)), weights=dat)
        return binarr, sigarr

# ny, nx = data.shape

x0, y0 = (3066.44,#+0.7+.31+.14+.069,
          2036.44)#+.8+.27+.107+.0416)
binarr, ULsigarr, URsigarr, BLsigarr, BRsigarr = quick_ringsum(data, x0, y0)

fig, ax = plt.subplots()
# ax.plot(binarr, ULsigarr, label='UL', lw=1)
# ax.plot(binarr, URsigarr, label='UR', lw=1)
# ax.plot(binarr, BLsigarr, label='BL', lw=1)
# ax.plot(binarr, BRsigarr, label='BR', lw=1)
ax.plot(binarr, ULsigarr + URsigarr, label="U")
ax.plot(binarr, BLsigarr + BRsigarr, label="B")
plt.legend(loc='upper left')
plt.show(block=False)

thres = 0.5 * np.max(ULsigarr + URsigarr)
i = np.where(ULsigarr + URsigarr > thres)[0]
ni = len(i)
# print ni
# j = np.arange(-25, 26, 1)
# sarr = stepsize * j
# UB = np.zeros(len(j))
# RL = np.zeros(len(j))
ns = 25
sarr = stepsize * np.arange(-ns, ns+1, 1)
# print sarr
UB = np.zeros(len(sarr))
RL = np.zeros(len(sarr))

ULsigarr *= 1.0
URsigarr *= 1.0
ULsigarr *= 1.0
URsigarr *= 1.0


# for idx, ix in enumerate(j):
# for ix in range(ns):
#     UB[ns+ix] = np.sum(((ULsigarr[i-ix]+URsigarr[i-ix]) - (BLsigarr[i+ix]+BRsigarr[i+ix]))**2) / (1.*ni-ix)
#     RL[ns+ix] = np.sum(((URsigarr[i-ix]+BRsigarr[i-ix]) - (ULsigarr[i+ix]+BLsigarr[i+ix]))**2) / (1.*ni-ix)
#     UB[ns-ix] = np.sum(((ULsigarr[i+ix]+URsigarr[i+ix]) - (BLsigarr[i-ix]+BRsigarr[i-ix]))**2) / (1.*ni-ix)
#     RL[ns-ix] = np.sum(((URsigarr[i+ix]+BRsigarr[i+ix]) - (ULsigarr[i-ix]+BLsigarr[i-ix]))**2) / (1.*ni-ix)
fig, ax = plt.subplots()
plt.plot(ULsigarr[i])
plt.show()
for idx, ix in enumerate(sarr):
    UB[idx] = np.sum((ULsigarr[i-ix]+URsigarr[i-ix] - BLsigarr[i+ix]-BRsigarr[i+ix])**2) / (1.*ni-np.abs(ix))
    RL[idx] = np.sum((URsigarr[i-ix]+BRsigarr[i-ix] - ULsigarr[i+ix]-BLsigarr[i+ix])**2) / (1.*ni-np.abs(ix))
fig, ax = plt.subplots()
ax.plot(sarr, UB, 'or', ms=1)
ax.plot(sarr, RL, 'og', ms=1)
RLfit = np.polyfit(sarr, RL, 2)
print RLfit, RLfit[1] / (2. * RLfit[0])
UBfit = np.polyfit(sarr, UB, 2)
print UBfit, UBfit[1] / (2. * UBfit[0])
plt.show(block=False)





# # x0, y0 = (3068.39, 2031.85)
# # x0, y0 = (3000., 2100.0)

# x = np.arange(1, nx+1, 1)
# y = np.arange(1, ny+1, 1)
# xx, yy = np.meshgrid(x-x0, y-y0)
# print xx


# xmin = xx[0, 0]
# xmax = xx[0, -1]
# ymin = yy[0, 0]
# ymax = yy[-1, 0]

# rmax = np.min(np.abs([xmin, xmax, ymin, ymax])) - 1
# print rmax
# ri = int(np.floor(rmax))

# R = np.sqrt(xx**2 + yy**2)
# binsize = 0.1

# imax = int((rmax ** 2. - 2. * rmax - 1.) / (1. + 2. * rmax) / binsize)
# binarr = np.fromfunction(lambda i: np.sqrt(2. * (i + 1.) * rmax * binsize + (i + 1.) * binsize ** 2.), (imax,),
#                          dtype='float64')

# xi0 = int(round(x0))  # np.abs(x - x0).argmin()
# yi0 = int(round(y0))  # np.abs(y - y0).argmin()

# # chop data matrix into a quadrant for now
# # ULdata = data[0:yi0+1, 0:xi0+1]
# # URdata = data[0:yi0+1, xi0:]
# # BLdata = data[yi0:, 0:xi0+1]
# # BRdata = data[yi0:, xi0:]

# ULdata = data[yi0-ri:yi0+1, xi0-ri:xi0+1]
# URdata = data[yi0-ri:yi0+1, xi0:xi0+ri+1]
# BLdata = data[yi0:yi0+ri+1, xi0-ri:xi0+1]
# BRdata = data[yi0:yi0+ri+1, xi0:xi0+ri+1]
# # fig, ax = plt.subplots()
# # ax.plot(data[yi0, 0:xi0+1], lw=1)
# # plt.show()
# print "UL", ULdata.shape
# print "UR", URdata.shape
# print "BL", BLdata.shape
# print "BR", BRdata.shape

# ULsigarr, _ = np.histogram(R[yi0-ri:yi0+1, xi0-ri:xi0+1], bins=np.concatenate((np.array([0.]), binarr)), weights=ULdata)
# URsigarr, _ = np.histogram(R[yi0-ri:yi0+1, xi0:xi0+ri+1], bins=np.concatenate((np.array([0.]), binarr)), weights=URdata)
# BLsigarr, _ = np.histogram(R[yi0:yi0+ri+1, xi0-ri:xi0+1], bins=np.concatenate((np.array([0.]), binarr)), weights=BLdata)
# BRsigarr, _ = np.histogram(R[yi0:yi0+ri+1, xi0:xi0+ri+1], bins=np.concatenate((np.array([0.]), binarr)), weights=BRdata)
# # plt.plot(np.concatenate((np.array([0.]), binarr))**2)
# # plt.show()
# print type(ULsigarr)
# ULsigarr *= 1.0
# URsigarr *= 1.0
# BLsigarr *= 1.0
# BRsigarr *= 1.0

# # print ULsigarr.shape
# # print URsigarr.shape
# # print BLsigarr.shape
# # print BRsigarr.shape
# thres = 0.3 * np.max(ULsigarr + URsigarr)
# i = np.where(ULsigarr + URsigarr > thres)[0]
# ni = len(i)
# print ni
# # j = np.arange(-25, 26, 1)
# # sarr = stepsize * j
# # UB = np.zeros(len(j))
# # RL = np.zeros(len(j))
# ns = 25
# sarr = stepsize * np.arange(-ns, ns+1, 1)
# print sarr
# UB = np.zeros(len(sarr))
# RL = np.zeros(len(sarr))
# # for idx, ix in enumerate(j):
# for ix in range(ns):
#     UB[ns+ix] = np.sum(((ULsigarr[i-ix]+URsigarr[i-ix]) - (BLsigarr[i+ix]+BRsigarr[i+ix]))**2) / (ni-ix)
#     RL[ns+ix] = np.sum(((URsigarr[i-ix]+BRsigarr[i-ix]) - (ULsigarr[i+ix]+BLsigarr[i+ix]))**2) / (ni-ix)
#     UB[ns-ix] = np.sum(((ULsigarr[i+ix]+URsigarr[i+ix]) - (BLsigarr[i-ix]+BRsigarr[i-ix]))**2) / (ni-ix)
#     RL[ns-ix] = np.sum(((URsigarr[i+ix]+BRsigarr[i+ix]) - (ULsigarr[i-ix]+BLsigarr[i-ix]))**2) / (ni-ix)
# plt.plot(sarr, UB, 'or', ms=1)
# plt.plot(sarr, RL, 'og', ms=1)
# # RLfit = np.polyfit(RL, j, 2)
# # print RLfit
# plt.show(block=False)

# fig, ax = plt.subplots()
# i = ax.imshow(data, cmap='gray')
# fig.colorbar(i)
# plt.show(block=False)
# # fig, ax = plt.subplots()
# # ax.imshow(ULdata, cmap='gray')
# # ax.set_aspect(1.0)
# # plt.plot(x0, y0, '.r', ms=1)
# # plt.show(block=False)
# # #
# fig, ax = plt.subplots()
# ax.plot(binarr, ULsigarr, label='UL', lw=1)
# # ax.plot(binarr, URsigarr, label='UR', lw=1)
# # ax.plot(binarr, BLsigarr, label='BL', lw=1)
# # ax.plot(binarr, BRsigarr, label='BR', lw=1)
# plt.legend(loc='upper left')
# plt.show(block=False)

# # fig, ax = plt.subplots(2,2)
# # ax[0][0].imshow(ULdata)
# # ax[0][1].imshow(URdata)
# # ax[1][0].imshow(BLdata)
# # ax[1][1].imshow(BRdata)


# # for axis in ax:
# #     for a in axis:
# #         a.set_aspect(1.0)

plt.show()


