import math
import os
import time

from collections import Counter

from cogent3.util import parallel


# the following environment variable is created by PBS on job execution
PBS_NCPUS = os.environ.get("PBS_NCPUS", None)
if PBS_NCPUS is None:
    raise RuntimeError("did not get cpu number from environment")

PBS_NCPUS = int(PBS_NCPUS)


def is_prime(n):
    # Postprocess the processor MPI rank to check your job got the resources
    # you requested
    r = parallel.get_rank()

    if n % 2 == 0:
        return False

    sqrt_n = int(math.floor(math.sqrt(n)))
    for i in range(3, sqrt_n + 1, 2):
        if n % i == 0:
            return False

    return r


PRIMES = (
    [
        112272535095293,
        112582705942171,
        112272535095293,
        115280095190773,
        115797848077099,
        117450548693743,
        993960000099397,
    ]
    * PBS_NCPUS
    * 20
)


def main():
    print(f"World size: {parallel.SIZE}\n")

    # Each worker will evaluate 20 prime numbers. This is just to slow the
    # script down!

    start = time.time()

    result = Counter(
        parallel.as_completed(is_prime, PRIMES, use_mpi=True, max_workers=PBS_NCPUS)
    )
    if sum(result.values()) != len(PRIMES):
        print(f" failed: {len(result)} != {len(PRIMES)} : {result=}")
    else:
        print(
            f"Time taken = {time.time() - start:.2f} seconds",
            f"CPU rank by number of jobs: {result}",
            sep="\n",
        )


if __name__ == "__main__":
    # This block is crucial! See
    # https://mpi4py.readthedocs.io/en/stable/mpi4py.futures.html
    # for why it needs to be done
    main()
