#!/usr/bin/python

import numpy as np
from sys import maxsize

class TreeNode:
    def __init__(self, label):
        self.label = label
        self.children = []
        self.distance_to_parent = 0.0

def calculateQ(d):
    r = d.shape[0]
    q = np.zeros((r, r))
    for i in range(r):
        for j in range(r):
            if i == j:
                q[i][j] = 0
            else:
                sumI = np.sum(d[i, :])
                sumJ = np.sum(d[j, :])
                q[i][j] = (r - 2) * d[i][j] - sumI - sumJ

    return q

def findLowestPair(q):
    r = q.shape[0]
    minVal = maxsize
    for i in range(0, r):
        for j in range(i, r):
            if (q[i][j] < minVal):
                minVal = q[i][j]
                minIndex = (i, j)
    return minIndex

def doDistOfPairMembersToNewNode(i, j, d):
    r = d.shape[0]
    sumI = np.sum(d[i, :])
    sumJ = np.sum(d[j, :])

    if r - 2 == 0:
        print("Warning: Division by zero avoided.")
        return (0, 0)  # Avoid division by zero

    dfu = (1. / (2. * (r - 2.))) * ((r - 2.) * d[i][j] + sumI - sumJ)
    dgu = (1. / (2. * (r - 2.))) * ((r - 2.) * d[i][j] - sumI + sumJ)

    return (dfu, dgu)

def calculateNewDistanceMatrix(f, g, d):
    r = d.shape[0]
    nd = np.zeros((r - 1, r - 1))

    ii = jj = 1
    for i in range(0, r):
        if i == f or i == g:
            continue
        for j in range(0, r):
            if j == f or j == g:
                continue
            nd[ii][jj] = d[i][j]
            jj += 1
        ii += 1
        jj = 1

    ii = 1
    for i in range(0, r):
        if i == f or i == g:
            continue
        nd[0][ii] = (d[f][i] + d[g][i] - d[f][g]) / 2.
        nd[ii][0] = (d[f][i] + d[g][i] - d[f][g]) / 2.
        ii += 1

    return nd

def build_tree_structure(d, labels):
    nodes = [TreeNode(label) for label in labels]

    while len(nodes) > 1:
        q = calculateQ(d)
        lowestPair = findLowestPair(q)
        i = lowestPair[0]
        j = lowestPair[1]

        pairDist = doDistOfPairMembersToNewNode(i, j, d)
        
        new_node = TreeNode(labels[i] + labels[j])
        new_node.children.append(nodes[i])
        new_node.children.append(nodes[j])
        new_node.distance_to_parent = pairDist[0]

        nodes.pop(j)
        nodes.pop(i)
        nodes.append(new_node)

        d = calculateNewDistanceMatrix(i, j, d)

    return nodes[0]

def print_tree(node, indent=0):
    print(" " * indent + node.label, ":", node.distance_to_parent)
    for child in node.children:
        print_tree(child, indent + 2)

def doNeighbourJoining(d, labels):
    tree_structure = build_tree_structure(d, labels)
    print("Final tree structure:")
    print_tree(tree_structure)

if __name__ == "__main__":
    distMatrix = np.array([
        [0.0000, 0.1890, 0.1100, 0.1130, 0.2150],
        [0.1890, 0.0000, 0.1790, 0.1920, 0.2110],
        [0.1100, 0.1790, 0.0000, 0.0941, 0.2050],
        [0.1130, 0.1920, 0.0940, 0.0000, 0.2140],
        [0.2150, 0.2110, 0.2050, 0.2140, 0.0000]
    ])

    labels = ["A", "B", "C", "D", "E"]

    doNeighbourJoining(distMatrix, labels)
