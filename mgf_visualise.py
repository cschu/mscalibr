import sys
import os
import math
import itertools as it
from collections import Counter

import numpy as np
import matplotlib
try:
    fig = matplotlib.pyplot.figure(figsize=(10.0, 3.0))
except:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

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

"""
Before I forget, let me write what I think you could do with Jan's data
For each peptide plot density of M/Z and if possible plot 10 peptides of similar mass together.
And plot 10 of such plots together (so you have batches of 100peptides together)
Plot first 100/1000 and see what we get.
"""

def n_subplots(n):
    if n == 1: return 1, 1
    if n == 2: return 1, 2
    if n == 3: return 1, 3
    if n == 4: return 2, 2
    if n == 5: return 2, 3
    if n == 6: return 2, 3
    if n == 7: return 3, 3
    if n == 8: return 3, 3
    if n == 9: return 3, 3
    if n == 10: return 3, 4
    if n == 11: return 3, 4
    if n == 12: return 3, 4
    if n == 13: return 4, 4
    if n == 14: return 4, 4
    if n == 15: return 4, 4
    if n == 16: return 4, 4
    return -1, 5


def plotDensities(fi, precMassBinSize=100):


    binf = maxMS2IntensityBinFunc
    massBins = {}
    spectrum_bin_map = {}
    with mgf.read(fi) as reader:
        spectra = sorted(((spectrum['params']['pepmass'][0], spectrum['params']['title'], binf(spectrum, params={'winsize': precMassBinSize})) for spectrum in reader), key=lambda x:(x[2], x[0]))
        for i, spectrum in enumerate(spectra):
            spectrum_bin_map[i] = spectrum[2]

    grid = n_subplots(len(set(spectrum_bin_map.values())))
    f, axes = plt.subplots(grid[0], grid[1], figsize=(16,16), sharey=True)
    sns.despine(left=True)
    sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 1.0})
    pal = sns.color_palette("Reds", n_colors=101)

    with mgf.read(fi) as reader:
        for i, spectrum in enumerate(reader):
            massBins[spectrum_bin_map[i]] = massBins.get(spectrum_bin_map[i], []) + [spectrum]

    row, col = 0, 0
    for mbin in sorted(massBins):
        print row, col, mbin, len(massBins[mbin])
        maxIntensity = max(max(spectrum['intensity array']) for spectrum in massBins[mbin])

        spectra = sorted(((max(spectrum['intensity array']), spectrum['params']['pepmass'][0], i) for i, spectrum in enumerate(massBins[mbin])), key=lambda x:x[0])

        for intensity, precMass, i in spectra:
            dp = sns.distplot(massBins[mbin][i]['m/z array'] - precMass, hist=False, color=pal[int(intensity/maxIntensity * 100 + 0.5)], ax=axes[row,col])
            axes[row, col].set_ylim(0, 0.004)
        if col == grid[1] - 1:
            row += 1
            col = 0
        else:
            col += 1

    plt.tight_layout()
    plt.savefig("test.png", dpi=300)





def plotHeatmap(fi, nSpectra=4000, fragBinsize=10, precursorBinsize=5):
    binf = maxMS2IntensityBinFunc

    # filter spectra so that at most 4000 spectra with log10(maxMS2Intensity) == 2 are used
    # spectra are sorted by log10(maxMS2Intensity), precursor mass
    with mgf.read(fi) as reader:
        spectra = sorted([(spectrum['params']['pepmass'][0], spectrum['params']['title'], binf(spectrum)) for spectrum in reader],
                         key=lambda x:(x[2], x[0]))
        spectra = [spectrum for spectrum in spectra if spectrum[2] == 2][:nSpectra]

    for spectrum in spectra:
        print spectrum

    # build dict spectrum-id => index (might have to be precursor mass)
    usedSpectra = dict([(spectrum[1], i) for i, spectrum in enumerate(spectra)])

    # extract spectrum masses/intensities from file according to usedSpectra list
    grid_d = {}

    X,Y = [], []

    with mgf.read(fi) as reader:
        for spectrum in reader:
            if spectrum['params']['title'] in usedSpectra:
                # find precursor bin (y-axis)
                binPrecursor = usedSpectra[spectrum['params']['title']] // precursorBinsize

                maxIntensity, maxIntensityFM = 0, 0
                for fragMass, fragIntensity in it.izip(spectrum['m/z array'], spectrum['intensity array']):
                    # find fragment mass bin (x-axis) for each fragment
                    binFrag = int(fragMass // fragBinsize)
                    if fragIntensity > maxIntensity:
                        maxIntensity = fragIntensity
                        maxIntensityFM = binFrag
                        maxIntensityBP = binPrecursor

                    # key = (binFrag, binPrecursor)
                    key = (binPrecursor, binFrag)

                    # grid_d holds 2D bins (precursor mass, fragment mass)
                    if key not in grid_d:
                        grid_d[key] = []
                    grid_d[key].append(fragIntensity)
                for i in xrange(int(maxIntensity + 0.5)):
                    X.append(int(maxIntensityFM + 0.5))
                    Y.append(int(maxIntensityBP + 0.5))




    with sns.axes_style("white"):
        # f, axes = plt.subplots(1, 1, figsize=(16,16), sharey=True)
        #fig = plt.figure()
        #sp = fig.add_subplot(1,1,1)
        sns.despine(left=True)
        sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 1.0})
        pl = sns.jointplot(x=np.array(X), y=np.array(Y), kind="hex", color="k")
        pl.savefig('hex.png', dpi=300)

        #plt.tight_layout()
        #plt.savefig("hex.png", dpi=300)
        del X
        del Y

    def mean(L):
        return sum(L) / len(L)
    nrows, ncols = 2000 // precursorBinsize, 2000 // fragBinsize

    for key in sorted(grid_d):
        print key, grid_d[key]

    # build colormesh grid - each cell is either 0.0 or the mean intensity
    grid = []
    for x in xrange(ncols):
        for y in xrange(nrows):
            if (x, y) in grid_d:
                grid.append(mean(grid_d[(x,y)]))
                del grid_d[(x,y)]
            else:
                grid.append(0.0)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolormesh(np.array(grid).reshape((nrows, ncols)), cmap=matplotlib.cm.seismic, vmin=min(grid), vmax=max(grid))
    """
    grid = []
    for y in xrange(nrows):
        row = []
        for x in xrange(ncols):
            if (x, y) in grid_d:
                row.append(mean(grid_d[(x,y)]))
                del grid_d[(x,y)]
            else:
                row.append(0.0)
        grid.append(row)

    grid = grid[::-1]
    grid = [item for sublist in grid for item in sublist]

    grid = np.array(grid).reshape((nrows, ncols))
    plt.imshow(grid,
               extent=(0, ncols, 0, nrows),
               interpolation='nearest',
               cmap=matplotlib.cm.seismic)
    """

    plt.tight_layout()
    plt.savefig(os.path.basename(fi) + '.colormap.png', dpi=300)
    plt.savefig(os.path.basename(fi) + '.colormap.svg', dpi=300)
    plt.close()

    pass

def plotHeatmapOld(fi, limit=1000):
    # X, Y, Z = [], [], []
    grid_d = {}
    xmin, xmax, ymin, ymax = None, None, None, None
    processedSpectra = 0
    with mgf.read(fi) as reader:
        for spectrum in reader:
            if processedSpectra > limit:
                break
            processedSpectra += 1
            x = int(spectrum['params']['pepmass'][0] * 100000)
            xmin, xmax = (min(x, xmin), max(x, xmax)) if xmin is not None else (x, x)
            for y, z in it.izip(spectrum['m/z array'], spectrum['intensity array']):
                y = int(y + 0.5)
                ymin, ymax = (min(y, ymin), max(y, ymax)) if ymin is not None else (y, y)
                grid_d[(x,y)] = z

            # X.extend([spectrum['params']['pepmass'][0] for i in spectrum['intensity array']])
            # Y.extend(spectrum['m/z array'])
            # Z.append(spectrum['intensity array'])
            # Z.extend(spectrum['intensity array'])

    grid = []
    for x in xrange(xmin, xmax + 1):
        for y in xrange(ymin, ymax + 1):
            if (x, y) in grid_d:
                grid.append(grid_d[(x,y)])
                del grid_d[(x,y)]
            else:
                grid.append(0.0)

    nrows, ncols = xmax - xmin + 1, ymax - ymin + 1
    grid2 = np.array(grid).reshape((nrows, ncols))
    grid = grid2

    #fig = plt.figure()
                                                    #plot = fig.add_subplot(111)

    plt.imshow(grid,
               extent=(xmin, xmax, ymin, ymax),
               interpolation='nearest',
               cmap=matplotlib.cm.seismic)

    # X = np.array(X)
    # Y = np.array(Y)
    # Z = np.matrix(Z)
    #range_ = min(Z), max(Z)

    #Z = matplotlib.cm.rainbow(map(lambda x:x/range_[1], Z))
    #plot.scatter(X, Y, color=Z) # pcolormesh(X, Y, Z)

    plt.tight_layout()
    plt.savefig(os.path.basename(fi) + '.colormap.png', dpi=300)
    plt.savefig(os.path.basename(fi) + '.colormap.svg', dpi=300)
    plt.close()







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
                bins[bin_] = [[], 0, []]
            bins[bin_][0].extend(spectrum['intensity array'])
            # to save space we only store the normalised masses
            # bins[bin_][1].extend(spectrum['m/z array'])
            bins[bin_][1] += 1
            bins[bin_][2].extend(spectrum['m/z array'] - pmass)

    # sort the peaks in each bin by mass, then intensity
    for bin_ in bins:
        nSpectra = bins[bin_][1]
        bins[bin_] = list(reversed(map(list, zip(*sorted(zip(bins[bin_][2], bins[bin_][0]))))))
        bins[bin_].insert(1, nSpectra)
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
        plot1.text(0.95, 0.95, 'fragments=%i\nprecursors=%i' % (len(bins1[bin_][2]), bins1[bin_][1]), verticalalignment='top', horizontalalignment='right',
                   style='italic', transform=plot1.transAxes, color='black', fontsize=8)
        plot2.text(0.95, 0.95, 'fragments=%i\nprecursors=%i' % (len(bins2[bin_][2]), bins2[bin_][1]), verticalalignment='top', horizontalalignment='right',
                   style='italic', transform=plot2.transAxes, color='black', fontsize=8)

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

    #getMassHistogram(sys.argv[1], binsize=50)
    #getMassHistogram(sys.argv[2], binsize=50)
    # plotHeatmap(sys.argv[1])
    plotHeatmap(sys.argv[2])
    #plotDensities(sys.argv[1])

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
