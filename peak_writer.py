import numpy as np
import json

n_wavelengths = int(raw_input("Number of wavelengths? "))
n_orders = int(raw_input("Number of orders for each wavelength? "))

w = []
LL = {}
RR = {} 

for i in range(n_wavelengths):
    wavelength = np.float64(raw_input("\nEnter wavelength (nm): "))
    w.append(wavelength)
    LL[wavelength] = []
    RR[wavelength] = []

    for j in range(n_orders):
        L = raw_input("\nLeft point of peak for order {0:d} for {1:f} nm: ".format(j, wavelength))
        R = raw_input("Right point of peak for order {0:d} for {1:f} nm: ".format(j, wavelength))

        LL[wavelength].append(np.float64(L))
        RR[wavelength].append(np.float64(R))

print LL
print RR


with open("newcalib_peak_locations.json", 'w') as outfile:
    locations = {'left pos': LL, 'right pos': RR}
    json.dump(locations, outfile, indent=4, separators=(',', ': '))


