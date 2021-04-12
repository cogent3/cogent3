import math
import os
import time

from cogent3.util import parallel


def is_prime(n):
    r = parallel.get_rank()
    print(f"Rank: {r}")

    if n % 2 == 0:
        return False

    sqrt_n = int(math.floor(math.sqrt(n)))
    for i in range(3, sqrt_n + 1, 2):
        if n % i == 0:
            return False

    return True


PRIMES = [
    112272535095293,
    112582705942171,
    112272535095293,
    115280095190773,
    115797848077099,
    117450548693743,
    993960000099397,
] * 4 # multiplying just to increase the amount of data to calculate 


print(f"World size: {parallel.size}\n")

start = time.time()
result = parallel.map(is_prime, PRIMES, max_workers=4)

print(f"{time.time() - start:.2f} seconds\n")
