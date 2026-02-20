"""CLI for invoking solve_qp."""

import argparse
import sys

import numpy as np
import scipy.io as sio
from qpsolvers import solve_qp

parser = argparse.ArgumentParser(description="Solve QP problem.")
parser.add_argument("infile", help="Input file with QP problem.")
parser.add_argument("outfile", help="Output file with QP solution.")

# Pass the specific args as the script may be called within an
# interpreter so can have variable number of args.
args = parser.parse_args(sys.argv[-2:])
infile = args.infile
outfile = args.outfile

qp = sio.loadmat(infile)["qp"][0][0]
QP = {name: qp[name].astype(np.float64) for name in qp.dtype.names}

X = solve_qp(
    P=QP["H"],
    q=QP["f"],
    G=QP["Aineq"],
    h=QP["bineq"],
    lb=None,
    ub=None,
    solver="cvxopt",
    initvals=None,  # do not initialize, in testing, initializing with a known solution is actually slower (why???)
)

if X is None:
    raise ValueError("Python CVXOPT QP solver failed to find a solution.")

sio.savemat(outfile, {"X": X})
