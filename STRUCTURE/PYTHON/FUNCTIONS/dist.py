import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import os

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

class distDens:
    def __init__(self):
        self.rmin = float(input("Minimum radius (nm): "))
        self.rmax = float(input("Maximum radius (nm): "))
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
        self.outdir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))+"/PLOTS"

    def read(self):
        freqData = []
        filtData = []
        for i in range(self.plts):
            dict = readData(self.paths[i])
            #pickle.dump(dict, open(('{}/{}.qt.pkl').format(self.outdir, self.labels[i]), "wb"))
            #print("Data saved as dictionary to {}/{}.qt.pkl".format(self.outdir, self.labels[i]))
            param, filt = splitData(dict,self.rmin,self.rmax)
            filtData.append(filt)
        # filtData = [i * 10 for i in filtData] # Convert to Å
        self.data = filtData

    def plot(self):
        fig, ax = plt.subplots(num=None, figsize=(6,4), dpi=300, facecolor='w', edgecolor='k')
        for i in range(self.plts):
            sns.distplot(self.data[i], ax=ax, hist=False, kde=True, hist_kws={'edgecolor':'black'}, kde_kws={'shade': True, 'linewidth': 2}, label=self.labels[i])
        plt.xlim(self.rmin,self.rmax)
        plt.legend(prop={'size': 10})
        plt.ylabel('Probability density', labelpad=10, fontsize=11)
        plt.xlabel('r (nm)', labelpad=10, fontsize=11)
        plt.tight_layout()
        plt.savefig('{}/{}.distPDF.png'.format(self.outdir,'_'.join(self.labels)))
        print ("Plot saved as {}/{}.distPDF.png".format(self.outdir,'_'.join(self.labels)))
