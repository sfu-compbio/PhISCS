#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==============================================================================
# Written by : Salem Malikic
# Modified by: Farid Rashidi
# Last Update: Apr 25, 2019
# ==============================================================================


from helperFunctions import *
import argparse
import errno
import subprocess
import sys


# COMMAND LINE ARGUMENTS PARSING
parser = argparse.ArgumentParser(description='PhISCS-B by Z3, aspino, Maxino, Open-WBO, QMaxSAT and MSCG solvers', add_help=True)
# Required arguments:
parser.add_argument('-SCFile', '--SCFile', required=True, type=str,
                    help='Path to single cell data matrix file')
parser.add_argument('-fn', '--fnProbability', required=True, type=float,
                    help='Probablity of false negative')
parser.add_argument('-fp', '--fpProbability', required=True, type=float,
                    help='Probablity of false positive')

# Optional:
parser.add_argument('-o', '--outDir', default='.', type=str,
                    help='Output directory')
parser.add_argument('-kmax', '--maxMutationsToEliminate', default=0, type=int,
                    help='Max number of mutations to be eliminated [default value is 0]')
parser.add_argument('-bulkFile', '--bulkFile', default=None, type=str,
                    help='Path to bulk data file')
parser.add_argument('-delta', '--delta', default=0.20, type=float,
                    help='Delta parameter accounting for VAF variance [default value is 0.20]')
parser.add_argument('-time', '--time', type=int, default=86400,
                    help='Max time (in seconds) allowed for the computation (supported only for Z3 solver) [default value is 24 hours]')
parser.add_argument('--drawTree', action='store_true',
                    help='Draw output tree by Graphviz')

parser.add_argument('-w', '--colEliminationWeight', default=0, type=float,
                    help='Weight of column elimination [default value is 0]')
parser.add_argument('-threads', '--threads', default=1, type=int,
                    help='Number of threads [default value is 1]')
parser.add_argument('--solver', type=str, default='Z3',
                    help='Solver can be Z3, aspino, Maxino, Open-WBO, QMaxSAT, MSCG or other (provided the path to exe file of solver)')


args = parser.parse_args()


# assert os.path.exists(args.outDir) == False, "ERROR!!! There already exists file or folder with name " + args.outDir + ". Exiting."
try:
    os.makedirs(args.outDir)
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(args.outDir):
        pass
    else:
        raise

usingBulk = False
if args.bulkFile:
    usingBulk = True

filename = os.path.splitext(os.path.basename(args.SCFile))[0]
outfile = os.path.join(args.outDir, filename)

if args.solver.lower() == 'z3':
    if usingBulk:
        cmds = ['python', 'PhISCS-B/csp_z3.py', '-f', args.SCFile, '-n', args.fnProbability, '-p', args.fpProbability, 
                '-w', args.colEliminationWeight, '-o', args.outDir, '-t 1', '--timeout', args.time, '-e', args.delta,
                '-b', args.bulkFile, '-m', args.maxMutationsToEliminate]
    else:
        cmds = ['python', 'PhISCS-B/csp_z3.py', '-f', args.SCFile, '-n', args.fnProbability, '-p', args.fpProbability, 
                '-w', args.colEliminationWeight, '-o', args.outDir, '-t 1', '--timeout', args.time, '-m', args.maxMutationsToEliminate]
else:
    cmds = ['PhISCS-B/csp_maxsat', '-f', args.SCFile, '-n', args.fnProbability, '-p', args.fpProbability, '-o', args.outDir, '-i']

cmd = ' '.join(str(v) for v in cmds)
# print(cmd)
subprocess.check_output(cmd, shell=True)

if args.drawTree:
    draw_tree("{}.CFMatrix".format(outfile), usingBulk, args.bulkFile)
