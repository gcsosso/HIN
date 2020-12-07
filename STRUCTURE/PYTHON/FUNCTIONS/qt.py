import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Read output file into nested dictionary, where each internal dictionary contains the data for each frame
def readData(inpath):
    data = {}
    with open(inpath) as fp:
        counter=0
        for i, line in enumerate(fp):
            frame = {'index':[],'filter':[],'order':[]}
            values=line.split()
            column=1
            for j in range(len(values)):
                if (column==1):
                    frame['index'].append(int(values[j]))
                    column=2
                elif(column==2):
                    frame['filter'].append(float(values[j]))
                    column=3
                elif(column==3):
                    frame['order'].append(float(values[j]))
                    column=1
            data[counter]=frame
            counter+=1
    fp.close()
    return(data)

# Binning (q(r))
def binData(data,nbins,rmin,rmax):
    dr=(rmax-rmin)/nbins
    half_dr=dr/2
    mesh=[None]*nbins
    for i in range(nbins):
        mesh[i]=rmin+((i+0.5)*dr)
    nframes=len(data)
    count=[0]*nbins
    sum=[0.0]*nbins
    avg=[None]*nbins
    for j in range(nframes):
        current=data[j]
        nfiltered=len(current['index'])
        for k in range(nfiltered):
            dist=current['filter'][k]
            order=current['order'][k]
            for l in range(nbins):
                if ((dist > mesh[l]-half_dr) and (dist <= mesh[l]+half_dr)):
                    count[l]+=1
                    sum[l]=sum[l]+order
    for r in range(nbins):
        try:
            avg[r]=sum[r]/count[r]
        except ZeroDivisionError:
            avg[r]=0
    return(avg,mesh)

# Processing (f(q))
def splitData(data,rmin,rmax):
    nframes=len(data)
    order_out=[]
    filter_out=[]
    for i in range(nframes):
        current=data[i]
        nfiltered=len(current['index'])
        for j in range(nfiltered):
            dist=current['filter'][j]
            order=current['order'][j]
            if ((dist > rmin) and (dist <= rmax)):
                order_out.append(order)
                filter_out.append(dist)
    return(order_out,filter_out)

# Distance (r or z) distribution of qT
class qtDist:
    def __init__(self):
        self.rmin = float(input("Minimum radius (nm): "))
        self.rmax = float(input("Maximum radius (nm): "))
        self.bins = int(input("Number of bins: "))
        inputPaths = []
        inputPaths.append(str(input("Path to input file (hin_structure.out.t_order): ")))
        inputLabels = []
        inputLabels.append(str(input("Data label: ")))
        subplots = int(input("Number of additional subplots: ") or 0)
        for i in range(subplots):
            path = str(input(("Path to input file (hin_structure.out.t_order) for subplot %s: ") %(i+1)))
            inputPaths.append(path)
            label = str(input(("Data label for subplot %s: ") %(i+1)))
            inputLabels.append(label)
        self.plts = subplots+1
        self.paths = inputPaths
        self.labels = inputLabels
        self.outdir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))+"/PLOTS"

    def read(self):
        avgData = []
        meshData = []
        for i in range(self.plts):
            dict = readData(self.paths[i])
            pickle.dump(dict, open(('{}/{}.qt.pkl').format(self.outdir, self.labels[i]), "wb"))
            print("Data saved as dictionary to {}/{}.qt.pkl".format(self.outdir, self.labels[i]))
            avg, mesh = binData(dict,self.bins,self.rmin,self.rmax)
            avgData.append(avg)
            meshData.append(mesh)
        self.avg = avgData
        self.mesh = meshData

    def plot(self):
        fig, ax = plt.subplots(num=None, figsize=(7,4), dpi=300, facecolor='w', edgecolor='k')
        for i in range(self.plts):
            sns.lineplot(y=self.avg[i],x=self.mesh[i],label=self.labels[i])
        plt.xlim(self.rmin,self.rmax)
        #plt.ylim(0.6,0.7) # User input for these?
        plt.ylabel('qT(d)', labelpad=10, fontsize=11)
        plt.xlabel('Distance (nm)', labelpad=10, fontsize=11)
        plt.legend(prop={'size': 10}) #loc='center left', bbox_to_anchor=(1.05, 0.7))
        plt.tight_layout()
        plt.savefig('{}/{}.qtDist.png'.format(self.outdir,'_'.join(self.labels)))
        print ("Plot saved as {}/{}.qtDist.png".format(self.outdir,'_'.join(self.labels)))

# Probability density distribution of qT
class qtDens:
    def __init__(self):
        self.rmin = float(input("Minimum radius (nm): "))
        self.rmax = float(input("Maximum radius (nm): "))
        inputPaths = []
        inputPaths.append(str(input("Path to input file (hin_structure.out.t_order): ")))
        inputLabels = []
        inputLabels.append(str(input("Data label: ")))
        subplots = int(input("Number of additional subplots: ") or 0)
        for i in range(subplots):
            path = str(input(("Path to input file (hin_structure.out.t_order) for subplot %s: ") %(i+1)))
            inputPaths.append(path)
            label = str(input(("Data label for subplot %s: ") %(i+1)))
            inputLabels.append(label)
        self.plts = subplots+1
        self.paths = inputPaths
        self.labels = inputLabels
        self.outdir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))+"/PLOTS"

    def read(self):
        paramData = []
        filtData = []
        for i in range(self.plts):
            dict = readData(self.paths[i])
            pickle.dump(dict, open(('{}/{}.qt.pkl').format(self.outdir, self.labels[i]), "wb"))
            print("Data saved as dictionary to {}/{}.qt.pkl".format(self.outdir, self.labels[i]))
            param, filt = splitData(dict,self.rmin,self.rmax)
            paramData.append(param)
            filtData.append(filt)
        self.data = paramData

    def plot(self):
        fig, ax = plt.subplots(num=None, figsize=(7,4), dpi=300, facecolor='w', edgecolor='k')
        for i in range(self.plts):
            sns.distplot(self.data[i], ax=ax, hist=False, kde=True, hist_kws={'edgecolor':'black'}, kde_kws={'shade': True, 'linewidth': 2}, label=self.labels[i])
        plt.xlim(-0.25,1)
        plt.legend(prop={'size': 10})
        plt.ylabel('Probability density', labelpad=10, fontsize=11)
        plt.xlabel('qT', labelpad=10, fontsize=11)
        plt.tight_layout()
        plt.savefig('{}/{}.qtDens.png'.format(self.outdir,'_'.join(self.labels)))
        print ("Plot saved as {}/{}.qtDens.png".format(self.outdir,'_'.join(self.labels)))
