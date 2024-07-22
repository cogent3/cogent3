import math
import time
from collections import Counter

from cogent3.util import parallel


def is_prime(n):
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
    * 4
    * 20
)  # multiplying just to increase the amount of data to calculate


def main():
    print(f"World size: {parallel.SIZE}\n")

    start = time.time()
    result = Counter(parallel.as_completed(is_prime, PRIMES, max_workers=4))

    if sum(result.values()) != len(PRIMES):
        print(f" failed: {len(result)} != {len(PRIMES)} : {result=}")
    else:
        print(
            f"Time taken = {time.time() - start:.2f} seconds",
            f"CPU rank by number of jobs: {result}",
            sep="\n",
        )


if __name__ == "__main__":
    main()
