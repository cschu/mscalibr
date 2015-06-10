import sys
import os
import math

import numpy as np
import matplotlib
try:
    fig = matplotlib.pyplot.figure(figsize=(10.0, 3.0))
except:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from scipy.stats import norm
from pyteomics import mgf

# functions to determine "bins" for plot-generation
def precMassBinFunc(spectrum, params={'winsize': 200}):
    # binning by precursor mass: bin = floor(precursor mass / window size)
    return int(spectrum['params']['pepmass'][0] // params['winsize'])
def precChargeBinFunc(spectrum, params={}):
    # binning by charge: bin = charge
    return spectrum['params']['charge'][0]
def precIntensityBinFunc(spectrum, params={}):
    # binning by precursor intensity: bin = log10(intensity) for intensity > 0
    # intensities less or equal than 0 are in the "None" bin
    intensity = spectrum['params']['pepmass'][1]
    return 'E%02i' % int(math.log(intensity, 10)) if intensity > 0 else 'NA'
def maxMS2IntensityBinFunc(spectrum, params={}):
    # binning by maximum intensity: bin = log10(maximum intensity)
    maxIntensity = max(spectrum['intensity array'])
    return int(math.log(maxIntensity, 10)) if maxIntensity > 0 else -1
def precIntensityAndMassBinFunc(spectrum, params={'winsize': 200}):
    # binning by precursor (intensity, mass)
    return precIntensityBinFunc(spectrum), precMassBinFunc(spectrum, params=params)
def precIntensityAndChargeBinFunc(spectrum, params={}):
    # binning by precursor (intensity, mass)
    return precIntensityBinFunc(spectrum), precChargeBinFunc(spectrum)
def maxMS2IntensityAndPrecMassBinFunc(spectrum, params={'winsize': 200}):
    # binning by max MS2 intensity and precursor mass
    return maxMS2IntensityBinFunc(spectrum), precMassBinFunc(spectrum, params=params)
def maxMS2IntensityAndChargeBinFunc(spectrum, params={}):
    # binning by max MS2 intensity and precursor charge (is that spectrum charge or precursor charge?)
    return maxMS2IntensityBinFunc(spectrum), precChargeBinFunc(spectrum)

def binSpectra(fi, binfunc, binparams):
    # this function reads an mgf file
    # and assigns bins to the spectra
    # according to the given bin-function and bin-parameters
    bins = {}
    with mgf.read(fi) as reader:
        for spectrum in reader:
            pmass = spectrum['params']['pepmass'][0]
            # here the bin-function is called
            bin_ = binfunc(spectrum, binparams)
            if bin_ not in bins:
                # raw_intensities, raw_masses, norm_masses
                bins[bin_] = [[], [], []]
            bins[bin_][0].extend(spectrum['intensity array'])
            # to save space we only store the normalised masses
            # bins[bin_][1].extend(spectrum['m/z array'])
            bins[bin_][2].extend(spectrum['m/z array'] - pmass)

    # sort the peaks in each bin by mass, then intensity
    for bin_ in bins:
        bins[bin_] = list(reversed(map(list, zip(*sorted(zip(bins[bin_][2], bins[bin_][0]))))))
        bins[bin_].insert(1, [])
    return bins

def makePlot(figure, where, title, x, y):
    pl = figure.add_subplot(where)
    pl.set_title(title, fontsize=12)
    pl.set_xlabel('m/z[Da]', fontsize=10)
    pl.set_ylabel('Intensity', fontsize=10)
    pl.scatter(x, y, s=1, c='b', marker='.')
    return pl

def getPlotLabels(fn, ptype, bin_, binparams):
    name = os.path.basename(fn)
    if ptype == "pmass":
        winsize = binparams['winsize']
        interval = (winsize * bin_, winsize * (bin_ + 1))
        title = '%s\n%i-%iDa' % ((name,) + interval)
        figname = '%s_%i-%i' % ((ptype,) + interval)
    elif ptype == "pcharge":
        title = '%s\ncharge=%i+' % (name, bin_)
        figname = '%s_%i+' % (ptype, bin_)
    elif ptype == "pintensity":
        title = '%s\nprecIntensity=%s' % (name, bin_)
        figname = '%s_%s' % (ptype, bin_)
    elif ptype == "maxms2intensity":
        title = '%s\nmaxMS2Intensity=%s' % (name, bin_)
        figname = '%s_%s' % (ptype, bin_)
    elif ptype == "pintensity+mass":
        winsize = binparams['winsize']
        interval = (winsize * bin_[1], winsize * (bin_[1] + 1)) # careful: bin_ is a tuple here
        title = '%s\nprecIntensity=%s,precMass=%i-%iDa' % ((name, bin_[0]) + interval)
        figname = '%s_%s+%i-%i' % ((ptype, bin_[0]) + interval)
    elif ptype == "pintensity+charge":
        title = '%s\nprecIntensity=%s,charge=%i+' % ((name,) + bin_)
        figname = '%s_%s+%i+' % ((ptype,) + bin_)
    elif ptype == "maxms2intensity+pmass":
        winsize = binparams['winsize']
        interval = (winsize * bin_[1], winsize * (bin_[1] + 1))
        title = '%s\nmaxMS2Intensity=%s,precMass=%i-%iDa' % ((name, bin_[0]) + interval)
        figname = '%s_%s+%i-%i' % ((ptype, bin_[0]) + interval)
    elif ptype == "maxms2intensity+pcharge":
        title = '%s\nmaxMS2Intensity=%s,charge=%i+' % ((name,) + bin_)
        figname = '%s_%s+%i+' % ((ptype,) + bin_)
    return title, figname

"""
def maxMS2IntensityAndPrecMassBinFunc(spectrum, params={'winsize': 200}):
    # binning by max MS2 intensity and precursor mass
    return maxMS2IntensityBinFunc(spectrum), precMassBinFunc(spectrum, params=params)
def maxMS2IntensityAndChargeBinFunc(spectrum, params={}):
    # binning by max MS2 intensity and precursor charge (is that spectrum charge or precursor charge?)
    return maxMS2IntensityBinFunc(spectrum), precChargeBinFunc(spectrum)
"""



def adjustAxes(plot1, plot2, axes=None):
    ax1, ax2 = plot1.axis(), plot2.axis()
    if axes is None:
	    axes = min(ax1[0], ax2[0]), max(ax1[1], ax2[1]), min(ax1[2], ax2[2]), max(ax1[3], ax2[3])
    plot1.axis(axes)
    plot2.axis(axes)
    pass


def getMassHistogram(fi, binsize=50):
    masses = []
    with mgf.read(fi) as reader:
        for spectrum in reader:
            pmass = spectrum['params']['pepmass'][0]
            masses.extend(spectrum['m/z array'] - pmass)
    # the histogram of the data with histtype='step'
    fig = plt.figure()
    plot = fig.add_subplot(111)

    n, bins, patches = plot.hist(masses, binsize, normed=1, histtype='stepfilled')
    plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)

    # add a line showing the expected distribution
    y = matplotlib.mlab.normpdf(bins, np.mean(masses), np.std(masses))
    l = plot.plot(bins, y, 'r--', linewidth=1.5)
    fig.tight_layout()
    fig.savefig(os.path.basename(fi) + '.massHist1.png', dpi=300)
    fig.savefig(os.path.basename(fi) + '.massHist1.svg', dpi=300)
    plt.close(fig)

    fig = plt.figure()
    plot = fig.add_subplot(111)

    n, bins, patches = plot.hist(masses, binsize, normed=1, histtype='bar')

    # add a line showing the expected distribution
    y = matplotlib.mlab.normpdf(bins, np.mean(masses), np.std(masses))
    l = plot.plot(bins, y, 'r--', linewidth=1.5)
    fig.tight_layout()
    fig.savefig(os.path.basename(fi) + '.massHist2.png', dpi=300)
    fig.savefig(os.path.basename(fi) + '.massHist2.svg', dpi=300)
    plt.close(fig)

    pass

def makePlots(fi1, fi2, ptype, binfunc, binparams):
    #bins1, bins2 = map(lambda x:binSpectra(x, binfunc, binparams), [fi1, fi2])
    bins1 = binSpectra(fi1, binfunc, binparams)
    bins2 = binSpectra(fi2, binfunc, binparams)
    for bin_ in sorted(set(bins1).intersection(set(bins2))):
        if bin_ is None:
            continue
        title1, figname = getPlotLabels(fi1, ptype, bin_, binparams)
        title2, figname = getPlotLabels(fi2, ptype, bin_, binparams)

        fig = plt.figure(figsize=(10.0, 3.0))
        plot1 = makePlot(fig, 121, title1, bins1[bin_][2], bins1[bin_][0])
        plot2 = makePlot(fig, 122, title2, bins2[bin_][2], bins2[bin_][0])
        xValues = bins1[bin_][2] + bins2[bin_][2]
        yValues = bins1[bin_][0] + bins2[bin_][0]
        adjustAxes(plot1, plot2, axes=(min(xValues), max(xValues), min(yValues), max(yValues)))

        # 0,0 lower-left, 1,1 upper-right
        plot1.text(0.95, 0.95, 'N=%i' % len(bins1[bin_][2]), verticalalignment='top', horizontalalignment='right',
                   style='italic', transform=plot1.transAxes, color='black', fontsize=10)
        plot2.text(0.95, 0.95, 'N=%i' % len(bins2[bin_][2]), verticalalignment='top', horizontalalignment='right',
                   style='italic', transform=plot2.transAxes, color='black', fontsize=10)

        fig.tight_layout()
        fig.savefig(figname + '.png', dpi=300)
        fig.savefig(figname + '.svg', dpi=300)
        plt.close(fig)

    pass

def getBinMembers(fi, binfunc, binparams):
    bins = {}
    with mgf.read(fi) as reader:
        for spectrum in reader:
            pmass = spectrum['params']['pepmass'][0]
            bin_ = binfunc(spectrum, binparams)
            if bin_ not in bins:
                # raw_intensities, raw_masses, norm_masses
                bins[bin_] = []
            bins[bin_].append(spectrum['params']['scans'])
    return bins

def writeBinMembers(fi, binfunc, binparams):
    members = getBinMembers(fi, binfunc, binparams)
    with open(os.path.basename(fi) + '.bmembers.txt', 'wb') as fo:
        for bin_ in sorted(members):
            fo.write('BIN %s\n' % bin_)
            for member in members[bin_]:
                fo.write('%s\n' % member)

def main():
    # uncomment to generate plot0})
    # makePlots(sys.argv[1], sys.argv[2], 'pcharge', precChargeBinFunc, None)
    # makePlots(sys.argv[1], sys.argv[2], 'pintensity', precIntensityBinFunc, None)
    # makePlots(sys.argv[1], sys.argv[2], 'maxms2intensity', maxMS2IntensityBinFunc, None)
    # makePlots(sys.argv[1], sys.argv[2], 'pintensity+mass', precIntensityAndMassBinFunc, binparams={'winsize': 200})
    # makePlots(sys.argv[1], sys.argv[2], 'pintensity+charge', precIntensityAndChargeBinFunc, None)
    # makePlots(sys.argv[1], sys.argv[2], 'maxms2intensity+pmass', maxMS2IntensityAndPrecMassBinFunc, binparams={'winsize': 200})
    # makePlots(sys.argv[1], sys.argv[2], 'maxms2intensity+pcharge', maxMS2IntensityAndChargeBinFunc, None)

    getMassHistogram(sys.argv[1], binsize=50)
    getMassHistogram(sys.argv[2], binsize=50)


    # uncomment to extract bin members
    # writeBinMembers(sys.argv[1], maxMS2IntensityBinFunc, None)
    # writeBinMembers(sys.argv[2], maxMS2IntensityBinFunc, None)
    pass


main()







"""
fi = '/mnt/slproj/MS_results_temp/Christian/HeLa_150506_20ng_r1.mgf'

raw_intensities = []
raw_masses, norm_masses = [], []

with mgf.read(fi) as reader:
    for spectrum in reader:
        print spectrum['params']['pepmass']
        raw_intensities.extend(spectrum['intensity array'])
        raw_masses.extend(spectrum['m/z array'])
        norm_masses.extend(spectrum['m/z array'] - spectrum['params']['pepmass'][0])
    print 'Done.'

raw.plot(raw_masses, raw_intensities)
norm.plot(norm_masses, raw_intensities)

plt.tight_layout()
plt.savefig(fi.replace('.mgf', '.png'), dpi=300)
"""
