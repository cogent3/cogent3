#!/usr/bin/env python

# Simple script to run another command a certain number of times,
# recording how long each run took, and writing the results out to a file.

import os
import os.path
import re
import string
import sys
import time


__author__ = "Peter Maxwell and  Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Edward Lang"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

# Values that affect the running of the program.
minimum_accepted_time = 2

iterations = int(sys.argv[1])
args = sys.argv[2:]

script_re = re.compile(".py$")
type = ""

for i in range(len(args)):
    if script_re.search(args[i], 1):
        script = args[i][0 : string.index(args[i], ".")]
    if args[i] == "mpirun":
        type = "parallel"

if not type:
    type = "serial"

output = "timing/" + script + "-" + str(int(time.time())) + "-" + type


def usage():
    pass


def standard_dev(numbers=None, mean=1):
    numbers = numbers or []
    import math

    sum = 0

    for i in range(len(numbers)):
        sum = sum + math.pow(numbers[i] - mean, 2)

    sigma = math.sqrt(sum / (len(numbers) - 1))

    return sigma


def main():
    if args:
        command = " ".join(map(str, args))
    else:
        usage()
        sys.exit()

    total_time = 0.0
    times = []

    print('Running "%s" %d times...' % (command, iterations))
    i = 0
    attempt = 0
    while i < iterations:
        start_time = time.time()
        os.system(command + " > " + output + "." + str(i))
        end_time = time.time() - start_time

        if end_time > minimum_accepted_time:
            times.append(end_time)
            total_time = total_time + end_time
            print("Time for run %d: %.3f seconds" % (i, end_time))
            i = i + 1
            attempt = 0
        else:
            print(f"Discarding probably bogus time: {end_time:.3f} seconds")
            attempt = attempt + 1
        if attempt == 5:
            print("Aborting early due to multiple errors")
            sys.exit(3)

    times.sort()
    mean = total_time / len(times)
    sd = standard_dev(times, mean)
    print("")
    print(f"Fastest time   : {times[0]:.3f}")
    print(f"Slowest time   : {times[len(times) - 1]:.3f}")
    print(f"Mean           : {mean:.3f}")
    print(f"Standard dev   : {sd:.3f}")
    print(f"Total time     : {total_time:.3f}")

    print("")

    corrected_total = 0.0
    corrected_times = []

    for i in range(len(times)):
        if abs(mean - times[i]) < sd:
            corrected_times.append(times[i])
            corrected_total = corrected_total + times[i]
        else:
            print(f"Discarding value '{times[i]:.3f}'")

    if len(times) != len(corrected_times):
        corrected_mean = corrected_total / len(corrected_times)
        corrected_sd = standard_dev(corrected_times, corrected_mean)

        print("")
        print("CORRECTED RESULTS")
        print(f"Fastest time   : {corrected_times[0]:.3f}")
        print(f"Slowest time   : {corrected_times[len(corrected_times) - 1]:.3f}")
        print(f"Mean           : {corrected_mean:.3f}")
        print(f"Standard dev   : {corrected_sd:.3f}")


if __name__ == "__main__":
    main()
