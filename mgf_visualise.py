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
from pyteomics import mgf

def precMassBinFunc(spectrum, params={'winsize': 200}):
    return int(spectrum['params']['pepmass'][0] // params['winsize'])
def precChargeBinFunc(spectrum, params={}):
    return spectrum['params']['charge'][0]
def precIntensityBinFunc(spectrum, params={}):
    intensity = spectrum['params']['pepmass'][1]
    return int(math.log(intensity, 10)) if intensity > 0 else None
def maxMS2IntensityBinFunc(spectrum, params={}):
    maxIntensity = max(spectrum['intensity array'])
    return int(math.log(maxIntensity, 10)) if maxIntensity > 0 else -1

def binSpectra(fi, binfunc, binparams):
    bins = {}
    with mgf.read(fi) as reader:
        for spectrum in reader:
            pmass = spectrum['params']['pepmass'][0]
            bin_ = binfunc(spectrum, binparams)
            if bin_ not in bins:
                # raw_intensities, raw_masses, norm_masses
                bins[bin_] = [[], [], []]
            bins[bin_][0].extend(spectrum['intensity array'])
            # bins[bin_][1].extend(spectrum['m/z array'])
            bins[bin_][2].extend(spectrum['m/z array'] - pmass)

    for bin_ in bins:
        bins[bin_] = list(reversed(map(list, zip(*sorted(zip(bins[bin_][2], bins[bin_][0]))))))
        bins[bin_].insert(1, [])
    return bins

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

def makePlot(figure, where, title, x, y):
    pl = figure.add_subplot(where)
    pl.set_title(title, fontsize=12)
    pl.set_xlabel('m/z[D]', fontsize=10)
    pl.set_ylabel('Intensity', fontsize=10)
    pl.scatter(x, y, s=1)
    return pl

def getPlotLabels(fn, ptype, bin_, binparams):
    name = os.path.basename(fn)
    if ptype == "pmass":
        winsize = binparams['winsize']
        interval = (winsize * bin_, winsize * (bin_ + 1))
        title = '%s\n%i-%iD' % ((name,) + interval)
        figname = '%s_%i-%i.png' % ((ptype,) + interval)
    elif ptype == "pcharge":
        title = '%s\ncharge=%i+' % (name, bin_)
        figname = '%s_%i+.png' % (ptype, bin_)
    elif ptype == "pintensity":
        title = '%s\nprecIntensity=E%02i' % (name, bin_)
        figname = '%s_E%02i.png' % (ptype, bin_)
    elif ptype == "maxms2intensity":
        title = '%s\nmaxMS2Intensity=E%02i' % (name, bin_)
        figname = '%s_E%02i.png' % (ptype, bin_)
    return title, figname

def adjustAxes(plot1, plot2):
    ax1, ax2 = plot1.axis(), plot2.axis()
    sharedAxis = min(ax1[0], ax2[0]), max(ax1[1], ax2[1]), min(ax1[2], ax2[2]), max(ax1[3], ax2[3])
    plot1.axis(sharedAxis)
    plot2.axis(sharedAxis)
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
        plot1, plot2 = makePlot(fig, 121, title1, bins1[bin_][2], bins1[bin_][0]), makePlot(fig, 122, title2, bins2[bin_][2], bins2[bin_][0])
        adjustAxes(plot1, plot2)

        fig.tight_layout()
        fig.savefig(figname, dpi=300)
        plt.close(fig)
    pass

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
