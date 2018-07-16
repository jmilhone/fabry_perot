from __future__ import print_function, division
from future.builtins import input
import json
import argparse
from os.path import abspath


def get_F_limits():
    input_string = "Please enter finesse (F) limits: "
    F_lim = get_and_parse_user_input(input_string, float, n_elements=2)
    return list(F_lim)


def get_Arel_limits():
    input_string = "Please enter Th I relative amplitude (Arel) limits: "
    Arel_lim = get_and_parse_user_input(input_string, float, n_elements=2)
    return list(Arel_lim)


def get_Ti_limits():
    input_string = "Please enter Ar II Ti limits: "
    Ti_lim = get_and_parse_user_input(input_string, float, n_elements=2)
    return list(Ti_lim)


def get_wavelengths():
    input_string = "Please enter all of extra Th I wavelengths to consider: "
    wavelengths = get_and_parse_user_input(input_string, float)
    return list(wavelengths)


def get_extra_rel_amplitudes(wavelengths):
    input_string = "Please enter the relative amplitude limits for {0:f} nm: "
    amps = []
    for w in wavelengths:
        rel_amp = get_and_parse_user_input(input_string.format(w), float, n_elements=2)
        amps.append(list(rel_amp))
    return amps


def get_and_parse_user_input(input_string, type_to_cast, n_elements=None):
    user_input = input(input_string)
    user_input = user_input.split(",")
    if len(user_input[0]) == 0:
        print('im here')
        return []
    if n_elements:
        # we know how many elements we need
        return (type_to_cast(x) for _, x in zip(range(n_elements), user_input))
    else:
        return (type_to_cast(x) for x in user_input)


def main(fname):
    print("*** All user inputs should be comma separated ***")
    F_lim = get_F_limits()
    Arel_lim = get_Arel_limits()
    Ti_lim = get_Ti_limits()
    w0 = get_wavelengths()
    if w0:
        amps = get_extra_rel_amplitudes(w0)
    else:
        amps = []
    config = {'F_lim': F_lim,
              'Arel_lim': Arel_lim,
              'Ti_lim': Ti_lim,
              'w_extra': w0,
              'Arel_extra': amps,
             }

    with open(fname, 'w') as f:
        json.dump(config, f, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Writes MultiNest finesse solver prior information to a json file for the finesse solver.')
    parser.add_argument('filename', type=str, help='file to write json data to')
    args = parser.parse_args()

    fname = abspath(args.filename)
    main(fname)


