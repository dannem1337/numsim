import random
import numpy as np


def throwDices(N):
    arr = []
    for i in range(N):
        rand = random.randint(1,6)
        np.append(arr, rand)
        print(arr)
    print(np.mean(arr))

throwDices(5)
