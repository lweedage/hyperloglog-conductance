import hashlib
import math
import struct
import numpy as np
import matplotlib.pyplot as plt
import time
import mmh3


def add(counter, hashed_bytes, b):
    i = left_most_bits(hashed_bytes, b)
    rho = leading_zeros(hashed_bytes, b)
    counter[i] = max(counter[i], rho)

def left_most_bits(hashed_bytes, b):
    return hashed_bytes >> (32-b)


def leading_zeros(hashed_bytes, b): #we take the leading zeros from the right instead of the left
    for n in range(32 - b):
        if hashed_bytes & (1 << n) != 0:
            return n + 1
    return 32 - b + 1


def union(old_counter, new_counter):
    np.maximum(old_counter, new_counter, new_counter)

def generic_hash(node):
    bytes = struct.pack('>Q', node)
    hashed_bytes = mmh3.hash(bytes, signed=False)
    return hashed_bytes

def combine_hash(hash1, hash2):
    bytes = struct.pack('>QQ', hash1, hash2)
    hashed_bytes = mmh3.hash(bytes, signed=False)
    return hashed_bytes

def node_hash(node):
    hashed_bytes = generic_hash(node)
    hashed_bytes_marker = generic_hash(3)
    hashed_bytes = combine_hash(hashed_bytes, hashed_bytes_marker)
    return hashed_bytes

def edge_hash_undirected(edge):
    edge = sorted(edge)
    hashed_bytes_1 = generic_hash(edge[0])
    hashed_bytes_2 = generic_hash(edge[1])
    hashed_bytes = combine_hash(hashed_bytes_1, hashed_bytes_2)
    hashed_bytes_marker = generic_hash(0)
    hashed_bytes = combine_hash(hashed_bytes, hashed_bytes_marker)
    return hashed_bytes

def edge_hash_directed(edge):
    edge = edge
    hashed_bytes_1 = generic_hash(edge[0])
    hashed_bytes_2 = generic_hash(edge[1])
    hashed_bytes = combine_hash(hashed_bytes_1, hashed_bytes_2)
    hashed_bytes_marker = generic_hash(1)
    hashed_bytes = combine_hash(hashed_bytes, hashed_bytes_marker)
    return hashed_bytes

def triangle_hash(triangle):
    triangle = sorted(triangle)
    hashed_bytes_1 = generic_hash(triangle[0])
    hashed_bytes_2 = generic_hash(triangle[1])
    hashed_bytes = combine_hash(hashed_bytes_1, hashed_bytes_2)
    hashed_bytes_3 = generic_hash(triangle[2])
    hashed_bytes = combine_hash(hashed_bytes, hashed_bytes_3)
    hashed_bytes_marker = generic_hash(2)
    hashed_bytes = combine_hash(hashed_bytes, hashed_bytes_marker)
    return hashed_bytes

def wedge_hash(wedge):
    hashed_bytes_1 = generic_hash(wedge[0])
    hashed_bytes_2 = generic_hash(wedge[1])
    hashed_bytes = combine_hash(hashed_bytes_1, hashed_bytes_2)
    hashed_bytes_3 = generic_hash(wedge[2])
    hashed_bytes = combine_hash(hashed_bytes, hashed_bytes_3)
    hashed_bytes_marker = generic_hash(3)
    hashed_bytes = combine_hash(hashed_bytes, hashed_bytes_marker)
    return hashed_bytes

mvlog = []

def size(counter, b, alpha):
    m = 2 ** b
    Z = np.sum(2.0**-counter)
    Z = 1 / Z
    E = alpha * m * m * Z
    # Small range correction:
    if E <= 5 / 2 * m:
        V = m - np.count_nonzero(counter)
        if V > 0:
            return mvlog[V]
        else:
            return E
    # Intermediate range - no correction
    elif E <= 1 / 30 * 2 ** 32:
        return E
    # Large range corrections
    elif E > 1 / 30 * 2 ** 32:
        print("Using large range correction")
        return -2 ** 32 * math.log(1 - E / (2 ** 32), 2)

def pre_compute_log(b):
    m = 2**b
    for V in range(m + 1):
        if V == 0:
            mvlog.append(None)
        else:
            mvlog.append(m * math.log(m/V))

def find_alpha(m):  # wikipedia
    if m < 17:
        return 0.673
    elif m < 33:
        return 0.697
    elif m < 65:
        return 0.709
    else:
        print(0.7213 / (1 + 1.079 / m))
        return 0.7213 / (1 + 1.079 / m)