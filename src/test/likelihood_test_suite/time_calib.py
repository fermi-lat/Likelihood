import time
from math import cos, sin, sqrt

NOMINAL = 0.45

def benchmark(count = 10):
    start = time.time()
    for _ in range(count):
        for i in range(1000000):
            x = pow(sin(i)**3.1415,cos(i**2))
        stop = time.time()
    return (stop - start)/10

if __name__ == "__main__":
    total = 0
    count = 10
    data = []
    for i in range(count):
        data.append(benchmark())
    ave = sum(data)/count

    std_dev = 0
    for i in range(count):
        std_dev += pow(ave-data[i],2)
    std_dev = sqrt(std_dev/count)
    print(f"average: {ave:.4f} +/- {std_dev:.4f}")