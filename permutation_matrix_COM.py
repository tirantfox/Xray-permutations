# test the center of mass of all permutation matrices.

import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

def get_matrix_COM(A):
    # get the center of mass of a matrix
    # A: a matrix
    # return: the center of mass of the matrix 
    # coordinates are 0, 0 at top left entry, n, n at bottom right
    n = A.shape[0]
    m = A.shape[1]
    moment = 0
    total_mass = 0
    for i in range(n):
        for j in range(m):
            total_mass += A[i, j]
            moment += A[i, j] * np.array([i, j])
    return moment / total_mass

def expected_MOI(n):
    return (n**2 - 1)/6

def get_matrix_MOI(A):
    # get the moment of inertia of a matrix
    # A: a matrix
    # return: the moment of inertia of the matrix 
    # coordinates are 0, 0 at top left entry, n, n at bottom right
    n = A.shape[0]
    m = A.shape[1]
    centroid = np.array([n-1, m-1]) / 2
    second_moment = 0
    total_mass = 0
    for i in range(n):
        for j in range(m):
            total_mass += A[i, j]
            second_moment += A[i, j] * np.linalg.norm((np.array([i, j]) - centroid))**2
    return second_moment / total_mass

def get_matrix_MOI_around_diagonal(A):
    # get the moment of inertia of a matrix around the main diagonal
    # A: a matrix
    # return: the moment of inertia of the matrix 
    # coordinates are 0, 0 at top left entry, n, n at bottom right
    n = A.shape[0]
    m = A.shape[1]
    centroid = np.array([n-1, m-1]) / 2
    second_moment = 0
    total_mass = 0
    for i in range(n):
        for j in range(m):
            total_mass += A[i, j]
            # for main diagonal: take dot product of the position relative to centroid with [1, 1]
            position = np.array([i, j]) - centroid
            second_moment += A[i, j] * np.dot(position, [1, 1])**2
    return second_moment / total_mass

def get_permutation_matrices(n):
    # a generator for all permutation matrices of size n
    # pick a position for the first 1 in the row, then recurse on the rest.
    if n == 0:
        yield np.zeros((0, 0))
    else:
        for i in range(n):
            A = np.zeros((n, n))
            A[0, i] = 1
            for B in get_permutation_matrices(n - 1):
                A[1:, :i] = B[:, :i]
                A[1:, i+1:] = B[:, i:]
                yield A



for n in range(1, 8):
    # print all permutation matrices of size n.
    print("\nn = ", n)
    print("Expected MOI: ", expected_MOI(n))

    counts = Counter()
    for A in get_permutation_matrices(n):
        # print(A)
        counts[np.round(get_matrix_MOI_around_diagonal(A), 4)] += 1

    # plot the counts as a histogram
    plt.bar(counts.keys(), counts.values())
    plt.show()



