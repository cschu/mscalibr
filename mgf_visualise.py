import sys
import os
import math

from pyteomics import mgf
import matplotlib.pyplot as plt
import numpy as np
import pylab

plt.figure(figsize=(10.0, 3.0))

raw = plt.subplot(1, 2, 1)
#raw.ylabel('raw')
# plt.ylabel('average')
# plt.plot(data.mean(axis=0))

norm = plt.subplot(1, 2, 2)
#norm.ylabel('normalised')
#plt.plot(data.max(axis=0))

def precMassBinFunc(spectrum, params={'winsize': 200}):
    return int(spectrum['params']['pepmass'][0] // params['winsize'])
def precChargeBinFunc(spectrum, params={}):
    return spectrum['params']['charge'][0]
def precIntensityBinFunc(spectrum, params={}):
    intensity = spectrum['params']['pepmass'][1]
    return int(math.log(intensity, 10)) if intensity > 0 else None

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

def makePlots(fi1, fi2, figfn, binfunc, binparams):
    bins1, bins2 = map(lambda x:binSpectra(x, binfunc, binparams), [fi1, fi2])
    for bin_ in sorted(set(bins1).intersection(set(bins2))):
        if bin_ is None:
            continue
        name1, name2 = map(os.path.basename, [fi1, fi2])
        if figfn == "pmass":
            winsize = binparams['winsize']
            interval = (winsize * bin_, winsize * (bin_ + 1))
            title1 = '%s\n%i-%iD' % ((name1,) + interval)
            title2 = '%s\n%i-%iD' % ((name2,) + interval)
            figname = '%s_%i-%i.png' % ((figfn,) + interval)
        elif figfn == "pcharge":
            title1 = '%s\ncharge=%i+' % (name1, bin_)
            title2 = '%s\ncharge=%i+' % (name2, bin_)
            figname = '%s_%i+.png' % (figfn, bin_)
        elif figfn == "pintensity":
            title1 = '%s\nprecIntensity=E%02i' % (name1, bin_)
            title2 = '%s\nprecIntensity=E%02i' % (name2, bin_)
            figname = '%s_E%02i.png' % (figfn, bin_)

        fig = plt.figure(figsize=(10.0, 3.0))
        plot1, plot2 = fig.add_subplot(1, 2, 1), fig.add_subplot(1, 2, 2)
        plot1.set_title(title1)
        plot2.set_title(title2)
        plot1.plot(bins1[bin_][2], bins1[bin_][0])
        plot2.plot(bins2[bin_][2], bins2[bin_][0])
        # CAREFUL: barplots are memory hogs and crash!!!
        #plot1.bar(bins1[bin_][2], bins1[bin_][0], width=0.1, linewidth=2, edgecolor='black')
        #plot2.bar(bins2[bin_][2], bins2[bin_][0], width=0.1, linewidth=2, edgecolor='black')
        fig.tight_layout()
        fig.savefig(figname, dpi=300)
    pass


makePlots(sys.argv[1], sys.argv[2], 'pmass', precMassBinFunc, {'winsize': 200})
makePlots(sys.argv[1], sys.argv[2], 'pcharge', precChargeBinFunc, None)
makePlots(sys.argv[1], sys.argv[2], 'pintensity', precIntensityBinFunc, None)







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
