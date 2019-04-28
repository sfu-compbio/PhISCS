import numpy as np
import sys
import os
from math import *

def readMatrixIntoHash(pathInputFile):
	assert os.path.exists(pathInputFile), "There does not exist file " + pathInputFile
	inputFile  = open(pathInputFile, "r")
	inputLines = inputFile.readlines()
	inputFile.close()
	assert len(inputLines) > 0, "ERROR. Input file " + pathInputFile + " is empty."

	headerEntries	= inputLines[0].strip().split()
	columnIDs	= headerEntries[1:len(headerEntries)]
	numColumns 	= len(columnIDs)
	assert numColumns > 0, "ERROR. First (header) line in " + pathSCFile + " is empty. Exiting!!!"
	inputLinesWithoutHeader = inputLines[1:len(inputLines)]
	D = {}
	for line in inputLinesWithoutHeader:
		lineColumns = line.strip().split()
		rowID = lineColumns[0].strip()
		assert rowID not in D.keys(), "ERROR in function readMatrixIntoHash. " + rowID + " is already in keys."
		D[rowID] = {}
		for i in range(numColumns):
			D[rowID][columnIDs[i]] = int(lineColumns[1+i])
	return D


def get_liklihood(inputSCMatrixFile, outputCFMatrixFile, fn, fp):
	D = readMatrixIntoHash(inputSCMatrixFile)
	E = readMatrixIntoHash(outputCFMatrixFile)
	alpha = float(fp)
	beta  = float(fn)
	missingEntryCharacter = 2 

	objectiveValueFromCFMatrix = 0.0
	cellIDs = E.keys()
	mutIDs  = E[cellIDs[0]].keys()
	dummyVariable = 1
	for i in cellIDs:
		for j in mutIDs:
			if D[i][j] == 1:
				if E[i][j] == 0:
					objectiveValueFromCFMatrix += log(alpha)
				elif E[i][j] == 1:
					objectiveValueFromCFMatrix += log(1-alpha)
			elif D[i][j] == 0:
				if E[i][j] == 1:
					objectiveValueFromCFMatrix += log(beta)
				elif E[i][j] == 0:
					objectiveValueFromCFMatrix += log(1-beta)
			elif D[i][j] == missingEntryCharacter:
				dummyVariable = 1
	return objectiveValueFromCFMatrix


def draw_tree(filename, addBulk, bulkfile):
    import pandas as pd
    import pygraphviz as pyg

    graph = pyg.AGraph(strict=False, directed=True)
    font_name = 'Avenir'

    class Node:
        def __init__(self, name, parent):
            self.name = name
            self.parent = parent
            self.children = []
            if parent:
                parent.children.append(self)

    def print_tree(node):
        graph.add_node(node.name, label=node.name, fontname=font_name, color='black', penwidth=3.5)
        for child in node.children:
            graph.add_edge(node.name, child.name)
            print_tree(child)

    def contains(col1, col2):
        for i in range(len(col1)):
            if not col1[i] >= col2[i]:
                return False
        return True

    def write_tree(matrix, names):
        i = 0
        while i < matrix.shape[1]:
            j = i + 1
            while j < matrix.shape[1]:
                if np.array_equal(matrix[:,i], matrix[:,j]):
                    matrix = np.delete(matrix, j, 1)
                    x = names.pop(j)
                    names[i] += '<br/><br/>' + x
                    j -= 1
                j += 1
            names[i] = '<'+names[i]+'>'
            i += 1

        rows = len(matrix)
        cols = len(matrix[0])
        dimensions = np.sum(matrix, axis=0)
        # ordered indeces
        indeces = np.argsort(dimensions)
        dimensions = np.sort(dimensions)
        mutations_name = []
        for i in range(cols):
            mutations_name.append(names[indeces[i]])

        root = Node(mutations_name[-1], None)
        mut_nod = {}
        mut_nod[mutations_name[cols-1]] = root

        i = cols - 2
        while i >=0:
            if dimensions[i] == 0:
                break
            attached = False
            for j in range(i+1, cols):
                if contains(matrix[:, indeces[j]], matrix[:, indeces[i]]):
                    node = Node(mutations_name[i], mut_nod[mutations_name[j]])
                    mut_nod[mutations_name[i]] = node
                    attached = True
                    break
            if not attached:
                node = Node(mutations_name[i], root)
                mut_nod[mutations_name[i]] = node
            i -=1
        print_tree(root)

    if addBulk:
        vafs = {}
        bulkMutations = readMutationsFromBulkFile(bulkfile)
        sampleIDs = bulkMutations[0].getSampleIDs()
        for mut in bulkMutations:
            temp_vaf = []
            for sample in sampleIDs:
                temp_vaf.append('<font color="blue">' + str(mut.getVAF(sampleID=sample)) + '</font>')
            vafs[mut.getID()] = '{} ({})'.format(mut.getID(), ','.join(temp_vaf))

    inp = np.genfromtxt(filename, skip_header=1, delimiter='\t')
    with open(filename, 'r') as fin:
        if addBulk:
            mutation_names = [vafs[x] for x in fin.readline().strip().split('\t')[1:]]
        else:
            mutation_names = fin.readline().strip().split('\t')[1:]
    sol_matrix = np.delete(inp, 0, 1)
    write_tree(sol_matrix, mutation_names)
    graph.layout(prog='dot')
    outputpath = filename[:-len('.CFMatrix')]
    graph.draw('{}.png'.format(outputpath))
