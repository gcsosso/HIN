import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Read output file into nested dictionary, where each internal dictionary contains the data for each frame
def readData(inpath,paramLoc,paramMax):
    data = {}
    with open(inpath) as fp:
        counter=0
        frame = {}
        for i, line in enumerate(fp):
            counter+=1
            values=line.split()
            if counter>paramMax:
                data[(i/paramMax)-1]=frame
                frame = {}
                counter=1
            if counter==1:
                frame['index']=list(map(int, values))
            elif counter==2:
                frame['filter']=list(map(float, values))
            elif counter==paramLoc:
                frame['order']=list(map(float, values))
    fp.close()
    return(data)

# Binning
def binData(data,bins,minFilt,maxFilt):
    dr=(maxFilt-minFilt)/bins
    half_dr=dr/2
    mesh=[None]*bins
    for i in range(bins):
        mesh[i]=minFilt+((i+0.5)*dr)
    frames=len(data)
    count=[0]*bins
    sum=[0.0]*bins
    avg=[None]*bins
    for j in range(frames):
        current=data[j]
        filtered=len(current['index'])
        for k in range(filtered):
            dist=current['filter'][k]
            order=current['order'][k]
            for l in range(bins):
                if ((dist > mesh[l]-half_dr) and (dist <= mesh[l]+half_dr)):
                    count[l]+=1
                    sum[l]=sum[l]+order
    for r in range(bins):
        try:
            avg[r]=sum[r]/count[r]
        except ZeroDivisionError:
            avg[r]=0
    return(avg,mesh)

# Processing
def splitData(data,minFilt,maxFilt):
    nframes=len(data)
    order_out=[]
    filter_out=[]
    for i in range(nframes):
        current=data[i]
        nfiltered=len(current['index'])
        for j in range(nfiltered):
            dist=current['filter'][j]
            order=current['order'][j]
            if ((dist > minFilt) and (dist <= maxFilt)):
                order_out.append(order)
                filter_out.append(dist)
    return(order_out,filter_out)

# Distance (r or z) distribution of F3/4
class fDist:
    def __init__(self):
        self.minFilt = float(input("Minimum radius (nm): "))
        self.maxFilt = float(input("Maximum radius (nm): "))
        self.bins = int(input("Number of bins: "))
        inputPaths = []
        inputPaths.append(str(input("Path to input file: ")))
        inputLabels = []
        inputLabels.append(str(input("Data label: ")))
        subplots = int(input("Number of additional subplots: ") or 0)
        for i in range(subplots):
            path = str(input(("Path to input file for subplot %s: ") %(i+1)))
            inputPaths.append(path)
            label = str(input(("Data label for subplot %s: ") %(i+1)))
            inputLabels.append(label)
        self.plts = subplots+1
        self.paths = inputPaths
        self.labels = inputLabels
        self.outDir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))+"/PLOTS"

    def read(self,paramLoc,paramMax):
        avgData = []
        meshData = []
        for i in range(self.plts):
            dict = readData(self.paths[i],paramLoc,paramMax)
            pickle.dump(dict, open(('{}/{}.F{}.pkl').format(self.outDir, self.labels[i],paramLoc), "wb"))
            print("Data saved as dictionary to {}/{}.F{}.pkl".format(self.outdir, self.labels[i],paramLoc))
            avg, mesh = binData(dict,self.bins,self.minFilt,self.maxFilt)
            avgData.append(avg)
            meshData.append(mesh)
        self.avg = avgData
        self.mesh = meshData

    def plot(self,paramLoc):
        fig, ax = plt.subplots(num=None, figsize=(7,4), dpi=300, facecolor='w', edgecolor='k')
        for i in range(self.plts):
            sns.lineplot(y=self.avg[i],x=self.mesh[i],label=self.labels[i])
        plt.xlim(self.minFilt,self.maxFilt)
        #plt.ylim(0.6,0.7) # User input for these?
        plt.ylabel('F{}(d)'.format(paramLoc), labelpad=10, fontsize=11)
        plt.xlabel('Distance (nm)', labelpad=10, fontsize=11)
        plt.legend(prop={'size': 10}) #loc='center left', bbox_to_anchor=(1.05, 0.7))
        plt.tight_layout()
        plt.savefig('{}/{}.F{}Dist.png'.format(self.outDir,'_'.join(self.labels),paramLoc))
        print("Plot saved as {}/{}.F{}Dist.png".format(self.outdir,'_'.join(self.labels),paramLoc))

# Probability density distribution of qT
class fDens:
    def __init__(self):
        self.minFilt = float(input("Minimum radius (nm): "))
        self.maxFilt = float(input("Maximum radius (nm): "))
        inputPaths = []
        inputPaths.append(str(input("Path to input file: ")))
        inputLabels = []
        inputLabels.append(str(input("Data label: ")))
        subplots = int(input("Number of additional subplots: ") or 0)
        for i in range(subplots):
            path = str(input(("Path to input file for subplot %s: ") %(i+1)))
            inputPaths.append(path)
            label = str(input(("Data label for subplot %s: ") %(i+1)))
            inputLabels.append(label)
        self.plts = subplots+1
        self.paths = inputPaths
        self.labels = inputLabels
        self.outDir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))+"/PLOTS"

    def read(self,paramLoc,paramMax):
        paramData = []
        filtData = []
        for i in range(self.plts):
            dict = readData(self.paths[i],paramLoc,paramMax)
            pickle.dump(dict, open(('{}/{}.F{}.pkl').format(self.outDir, self.labels[i],paramLoc), "wb"))
            print("Data saved as dictionary to {}/{}.F{}.pkl".format(self.outdir, self.labels[i],paramLoc))
            param, filt = splitData(dict,self.minFilt,self.maxFilt)
            paramData.append(param)
            filtData.append(filt)
        self.data = paramData

    def plot(self,paramLoc):
        fig, ax = plt.subplots(num=None, figsize=(7,4), dpi=300, facecolor='w', edgecolor='k')
        for i in range(self.plts):
            sns.distplot(self.data[i], ax=ax, hist=False, kde=True, hist_kws={'edgecolor':'black'}, kde_kws={'shade': True, 'linewidth': 2}, label=self.labels[i])
        #plt.xlim(-0.25,1)
        plt.legend(prop={'size': 10})
        plt.ylabel('Probability density'.format(paramLoc), labelpad=10, fontsize=11)
        plt.xlabel('F{}'.format(paramLoc), labelpad=10, fontsize=11)
        plt.tight_layout()
        plt.savefig('{}/{}.F{}Dens.png'.format(self.outDir,'_'.join(self.labels),paramLoc))
        print("Plot saved as {}/{}.F{}Dens.png".format(self.outdir,'_'.join(self.labels),paramLoc))
