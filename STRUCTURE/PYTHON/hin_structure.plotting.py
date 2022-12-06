import argparse
import warnings
import sys
import os
warnings.filterwarnings("ignore")
pathToScript=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
pathToFunctions=pathToScript+"/FUNCTIONS"
sys.path.append(pathToScript)
sys.path.append(pathToFunctions)

from FUNCTIONS.qt import *
from FUNCTIONS.f3_f4 import *
from FUNCTIONS.dist import *

# To do:
# - Add new plotting type: 'distance', with options for plotting a pair correlation function or probability density (e.g. for Tom's work)
# - Inputs for x and y labels/ other style options (if additional flag -s is included?)
# - Additional informaion for each plotting type: #parser.add_argument("-i","--info", action='store_true', help="Show more information about the plotting parameters for selected plot")

parser = argparse.ArgumentParser(description="Post-processing and plotting tool for HIN order parameters")
parser.add_argument("parameter", type=str, help="Specify the order parameter to plot", choices=['qt','q3','f3','f4','distance'])

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

elif (args.parameter in ["f3", "f4", "q3"]):
    if args.parameter == "q3":
        paramLoc=3
        paramMax=3
    elif args.parameter == "f3":
        paramLoc=3
        paramMax=4
    elif args.parameter == "f4":
        paramLoc=4
        paramMax=4
    print("\n 0   Radial/Z distribution\n 1   Probability density\n")
    plotType = int(input("Select plot: "))

    if plotType == 0:
        fData = fDist()
        fData.read(paramLoc,paramMax,args.parameter)
        fData.plot(paramLoc,args.parameter)

    elif plotType == 1:
        qtData = fDens()
        qtData.read(paramLoc,paramMax,args.parameter)
        qtData.plot(paramLoc,args.parameter)

    else:
        print("Must select either: Radial/Z distribution [enter 0] ; Probability density [enter 1]")
        exit()

elif args.parameter == "distance":
    print("\n 0   Probability density\n")
    plotType = int(input("Select plot: "))

    if plotType == 0:
        distData = distDens()
        distData.read()
        distData.plot()

    else:
        print("Must select: Probability density [enter 0]")
        exit()
