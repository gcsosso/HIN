import argparse
import warnings
import sys
import os
warnings.filterwarnings("ignore")
pathToFunctions=os.getcwd()+"/FUNCTIONS"
sys.path.append(os.getcwd())
sys.path.append(pathToFunctions)

from FUNCTIONS.qt import *
from FUNCTIONS.f3_f4 import *

# To do:
# - Add new plotting type: 'distance', with options for plotting a pair correlation function or probability density (e.g. for Tom's work)
# - Inputs for x and y labels/ other style options (if additional flag -s is included?)
# - Additional informaion for each plotting type: #parser.add_argument("-i","--info", action='store_true', help="Show more information about the plotting parameters for selected plot")

parser = argparse.ArgumentParser(description="Post-processing and plotting tool for HIN order parameters")
parser.add_argument("parameter", type=str, help="Specify the order parameter to plot", choices=['qt','f3','f4'])

args = parser.parse_args()

if args.parameter == "qt":
    print("\n 0   Radial/Z distribution\n 1   Probability density\n")
    plotType = int(input("Select plot: "))

    if plotType == 0:
        qtData = qtDist()
        qtData.read()
        qtData.plot()

    elif plotType == 1:
        qtData = qtDens()
        qtData.read()
        qtData.plot()

    else:
        print("Must select either: Radial/Z distribution [enter 0] ; Probability density [enter 1]")
        exit()

elif (args.parameter == "f3" or "f4"):
    if args.parameter == "f3":
        paramLoc=3
    elif args.parameter == "f4":
        paramLoc=4
    print("\n 0   Radial/Z distribution\n 1   Probability density\n")
    plotType = int(input("Select plot: "))

    if plotType == 0:
        fData = fDist()
        fData.read(paramLoc,4)
        fData.plot(paramLoc)

    elif plotType == 1:
        qtData = fDens()
        qtData.read(paramLoc,4)
        qtData.plot(paramLoc)

    else:
        print("Must select either: Radial/Z distribution [enter 0] ; Probability density [enter 1]")
        exit()
