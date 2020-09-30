import argparse
import warnings
warnings.filterwarnings("ignore")

# To do:
# - Add default inputs
# - Inputs for y-axis min and max
# - Inputs for other style options (if additional flag is included?)
# - Additional informaion for each plotting type: #parser.add_argument("-i","--info", action='store_true', help="Show more information about the plotting parameters for selected plot")
# - Add DATA folder in PYTHON directory where output files are placed for plotting (without specifying input path)

parser = argparse.ArgumentParser(description="Post-processing and plotting tool for HIN_structure.")
parser.add_argument("parameter", type=str, help="Specify the order parameter to plot.", choices=['qt'])

args = parser.parse_args()

if args.parameter == "qt":
    from FUNCTIONS.qt import *
    plotType = int(input("Select plot:\n 0   Probability density, f(qT)\n 1   Radial distribution, qT(r)\n\n"))

    if plotType == 0:
        qtData = f_qt()
        qtData.read()
        qtData.plot()

    elif plotType == 1:
        qtData = qt_r()
        qtData.read()
        qtData.plot()
