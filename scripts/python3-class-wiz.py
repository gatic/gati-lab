#### Script to investigate convergence behaviour in relion 3D classification
#### Copyright Cornelius Gati 2020 - SLAC Natl Acc Lab - cgati@stanford.edu

import os
import sys
import numpy as np
import scipy
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import operator
import collections
import matplotlib.backends.backend_pdf

################INPUT

print('Please specify:')
print('--f 		folder path 							(default: current folder)')
print('--root 		root name 							(default: \'run\')')
print('--o		output pdf name 						(default: output.pdf)')
print('--filt 		\'true\' or \'false\' obtain filtered.star file 			(default: false)')
print('--sigmafac 	cutoff for \'filt\', how many sigma above mean 			(default: 1)')
print('--mic		minimum cutoff for CTFFIND/Gctf resolution estimate 		(default: none)')

folder = '.'
rootname = 'run'
output = 'output.pdf'
plottype = 'bar'
micfilt = ''
filtstar = 'false'
sigmafac = 1

for (si, s) in enumerate(sys.argv):
    if s == '--f':
        folder = sys.argv[si + 1]

    if s == '--root':
        rootname = sys.argv[si + 1]

    if s == '--o':
        output = sys.argv[si + 1]

    if s == '--plot':
        plottype = sys.argv[si + 1]

    if s == '--filt':
        filtstar = sys.argv[si + 1]

    if s == '--sigmafac':
        sigmafac = sys.argv[si + 1]

    if s == '--mic':
        micfilt = sys.argv[si + 1]

#### List all files in folder and sort by name

filesdir = sorted(os.listdir(folder))
unwanted = []

##Check number of iterations

iterationlist = []
iterationtemp = []
iterations = []
filesdir = sorted(filesdir)

for datafile in filesdir:
    if 'data.star' in datafile and 'sub' not in datafile:
        if 'it001' in datafile and len(folder) == 0:  # Get rootname by looking at first iteration if not specified
            rootname = datafile.split('_it')[0]
        if 'ct' in datafile.split('_')[-3]:
            if int((datafile.split('_')[-2])[2:]) \
                != int((datafile.split('_')[-3])[2:]):
                iterationtemp.append(datafile)
        if 'ct' not in datafile.split('_')[-3]:
            iterationtemp.append(datafile)

        iterations.append(int((datafile.split('_')[-2])[2:]))

### Open PDF for output of figures

pdf = matplotlib.backends.backend_pdf.PdfPages('%s' % output)

## List of all input files used

iterationlist = sorted(iterationtemp, key=lambda x: x.split('_')[-2])
if len(iterations) == 0:
    print('')
    print('I cannot find any data.star files in the provided folder!')
    sys.exit()

iterations = max(iterations) + 1

##Print used input data.star files

print('')
for files in iterationlist:
    print('Using %s as input' % files)

##Check number of particles, number of classes, number of micrographs

part = 0
classes = []
micrographlist = []
checklist = []
checklistcol = []
anglerot = []
rescol = []

with open('%s/%s_it016_data.star' % (folder, rootname), 'r') as f:
    for l in f:
        if l[0] == '_' and 10 < len(l) < 50:  # check header
            if l.split()[0] == '_rlnClassNumber':  # check header for ClassNumber column
                classcolumn = int((l.split()[1])[1:]) - 1
            if l.split()[0] == '_rlnMicrographName':  # check header for MicrographName column
                miccolumn = int((l.split()[1])[1:]) - 1
            if l.split()[0] == '_rlnImageName':  # check header for ParticleName column
                particolumn = int((l.split()[1])[1:]) - 1
            if l.split()[0] == '_rlnAngleRot':  # check header for ParticleName column
                anglerot = int((l.split()[1])[1:]) - 1
            if l.split()[0] == '_rlnCtfMaxResolution':  # check header for CtfMaxResolution column
                rescol = int((l.split()[1])[1:]) - 1

            if '_rln' in l.split()[0]:
                checklist.append(int((l.split()[1])[1:]) - 1)
                checklistcol.append(l.split()[0])

        if '@' in l and l.split()[0] != 'opticsGroup1':
            part += 1
            classes.append(int(l.split()[classcolumn]))
            micrographlist.append(l.split()[miccolumn])

classes = max(classes)

print('')
print(('Plots will be generated for the following columns:',
       checklistcol))

## Initial colorbar

fig = plt.figure(figsize=(9, 11))
ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
cmap = plt.get_cmap('jet', int(classes) + 1)
norm = matplotlib.colors.Normalize(vmin=1, vmax=int(classes) + 2)

cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
        orientation='horizontal', label='Color for each class')

# colorbar stuff

labels = np.arange(0, classes + 2)

loc = labels + 0.5
cb1.set_ticks(loc)

a = ax1.get_xticks().tolist()
a[1] = 'no class'
a[2:] = labels[2:] - 1
result = ax1.set_xticklabels(a)

pdf.savefig()
plt.close()
###### Go into each iteration_data.star file and read in information, such as particle class assignments etc.

groupnumarray = np.zeros((part, iterations), dtype=np.double)  # Class assignments over all iterations
checkarray = np.empty((part, len(checklist)), dtype=np.double)  # Stats of final iteration
check = np.empty((part, 2, iterations), dtype=np.double)  # Stats of final iteration
changes = []

micnumdict = collections.defaultdict(list)

print('')
for datafile in iterationlist:

    miciterdict = collections.defaultdict(list)
    if 'ct' in datafile:
        iteration = int((datafile.split('_')[-2])[2:])
    if 'ct' not in datafile:
        iteration = int((datafile.split('_')[-2])[2:])

    particle = 0
    changesum = 0
    if int(iteration) > 1:
        with open('%s/%s' % (folder, datafile), 'r') as f:
            for l in f:
                if '@' in l:
                    groupnum = l.split()[classcolumn]  # # Class number
                    if int(groupnum) > int(classes):
                        groupnum = 0
                    micrograph = l.split()[miccolumn]  # # Micrograph name
                    particlename = l.split()[particolumn]  # # Particle Name

                    groupnumarray[particle, iteration] = groupnum
                    miciterdict[micrograph].append(groupnum)
                    if int(iteration) == iterations - 1:  # last iteration: all columns of star file to checkarray
                        for i in range(0, int(checklist[-1] + 1)):
                            if 'mrc' not in l.split()[i]:  # No filenames in checklist
                                colval = l.split()[i]
                                if 'group' in l.split()[i]:
                                    colval = int((l.split()[i])[-2:])
                                checkarray[particle, i] = colval
                            else:
                                checkarray[particle, i] = 0

                    if int(iteration) >= iterations - 2:
                        for i in range(0, len(checklist) - 1):
                            check[particle, 0, iteration] = \
                                l.split()[anglerot]
                            check[particle, 1, iteration] = \
                                l.split()[classcolumn]

                    if groupnumarray[particle, iteration] \
                        == groupnumarray[particle, int(iteration) - 1]:
                        change = 'nan'
                        change1 = 0
                    if groupnumarray[particle, iteration] \
                        != groupnumarray[particle, int(iteration) - 1]:
                        change = groupnum
                        change1 = 1
                        changesum += 1

                    particle += 1
    changes.append(changesum)

    print('Iteration %s: %s particles changed class assignments' \
        % (iteration, changesum))

    ## After each iteration create Histogram of class assignments for each micrograph

    miciterdict = \
        collections.OrderedDict(sorted(list(miciterdict.items()),
                                key=operator.itemgetter(0)))
    for (key, value) in list(miciterdict.items()):
        michisto = []
        for numb in value:
            michisto.append(float(numb))
        classes = int(classes)
        michisto = np.histogram(michisto, bins=classes, range=(1,
                                classes + 1))
        micnumdict[key].append(michisto[0].tolist())

######## Plot rotational and translational accuracy over each iteration

rotationcol = 2
translationcol = 3
rotation = np.zeros((int(classes) + 1, iterations), dtype=np.double)
translation = np.zeros((int(classes) + 1, iterations), dtype=np.double)

for datafile in iterationlist:
    iteration = int((datafile.split('_')[-2])[2:])
    with open('%s/%s_model.star' % (folder, datafile[:-10]), 'r') as f:
        for l in f:
            if '_rlnAccuracyRotations' in l:
                rotationcol = int(l.split('#')[-1]) - 1
            if '_rlnAccuracyTranslations' in l:
                translationcol = int(l.split('#')[-1]) - 1
            if 'class' in l and 'mrc' in l and 'mrcs' not in l:
                classnum = int((l.split('.mrc')[0])[-3:])
                rotation[classnum, iteration] = l.split()[rotationcol]
                translation[classnum, iteration] = \
                    l.split()[translationcol]

rotation = np.array(rotation[1:])
translation = np.array(translation[1:])

if len(set(rotation[0])) == 1:
    print('You did not perform image alignment during classification - skipping these two plots!')

if len(set(rotation[0])) > 1:

# Rotational

    cmap = plt.get_cmap('jet', int(classes) + 1)
    plt.figure(num=None, dpi=80, facecolor='white',figsize=(9, 11))
    plt.title('RotationalAccuracy', fontsize=16, fontweight='bold')
    plt.xlabel('Iteration #', fontsize=13)
    plt.ylabel('RotationalAccuracy', fontsize=13)
    plt.grid()
    colors = np.arange(1, int(classes) + 1)
    d = 0
    for (c, r) in zip(colors, rotation):
        d = c
        plt.plot(r[:], linewidth=3, color=cmap(c), label='Class %s' % d)
    ticks = np.arange(2, iterations)
    plt.xlim(2, iterations)
    plt.legend(loc='best')
    pdf.savefig()
    plt.close()
# Translational

    cmap = plt.get_cmap('jet', int(classes) + 1)
    plt.figure(num=None, dpi=80, facecolor='white', figsize=(9, 11))
    plt.title('TranslationalAccuracy', fontsize=16, fontweight='bold')
    plt.xlabel('Iteration #', fontsize=13)
    plt.ylabel('TranslationalAccuracy', fontsize=13)
    plt.grid()
    colors = np.arange(1, int(classes) + 1)
    d = 0
    for (c, t) in zip(colors, translation):
        d = c
        plt.plot(t[:], linewidth=3, color=cmap(c), label='Class %s' % d)
    ticks = np.arange(2, iterations)
    plt.xlim(2, iterations)
    plt.legend(loc='best')
    pdf.savefig()
plt.close()
###########################################################################

### Sort group assignment array column by column

sortindices = np.lexsort(groupnumarray[:, 1:].T)
groupnumarraysorted = groupnumarray[sortindices]

### Heat map of group sizes

H = groupnumarraysorted[:, :]
cmap = plt.get_cmap('jet', int(classes) + 1)
norm = matplotlib.colors.Normalize(vmin=0, vmax=int(classes) + 1)
plt.figure(num=None, dpi=120, facecolor='white', figsize=(9, 11))
plt.title('Class assignments of each particle', fontsize=16,
          fontweight='bold')
plt.xlabel('Iteration #', fontsize=13)
plt.ylabel('Particle #', fontsize=13)
mat = plt.imshow(H, aspect='auto', interpolation='nearest', cmap=cmap,
                 norm=norm)

# colorbar stuff

labels = np.arange(0, classes + 2)

# labels[-1] = 'unassigned'

cb1 = plt.colorbar(mat, ticks=labels)
loc = labels + .5

# print loc

cb1.set_ticks(loc)
a[0] = 'no class'
a[1:-1] = labels[1:-1]
cb1.set_ticklabels(a)
cb1.ax.tick_params(labelsize=16)
cb1.set_label('Class #')
plt.xlim(2, iterations - .5)

pdf.savefig()
plt.close()

# plt.show()

### Plot heat map of the last 5 iterations (close-up)

H = groupnumarraysorted[:, :]
cmap = plt.get_cmap('jet', int(classes) + 1)
norm = matplotlib.colors.Normalize(vmin=0, vmax=int(classes) + 1)
plt.figure(num=None, dpi=120, facecolor='white', figsize=(9, 11))
plt.title('Class assignments of each particle - last 5 iterations',
          fontsize=16, fontweight='bold')
plt.xlabel('Iteration #', fontsize=13)
plt.ylabel('Particle #', fontsize=13)
mat = plt.imshow(H, aspect='auto', interpolation='nearest', cmap=cmap,
                 norm=norm)

# colorbar stuff

labels = np.arange(0, classes + 2)
cb1 = plt.colorbar(mat, ticks=labels)
loc = labels + .5
cb1.set_ticks(loc)
a[0] = 'no class'
a[1:-1] = labels[1:-1]
cb1.set_ticklabels(a)
cb1.ax.tick_params(labelsize=16)
cb1.set_label('Class #')
plt.xlim(iterations - 6.5, iterations - .5)
pdf.savefig()
plt.close()
# plt.show()

#########################################################################################################################################################
### Sort group assignment array column by column

check = check[sortindices]

#### Jumper analysis

checkdict = collections.defaultdict(list)
checktest = []
labelsY = []

for c in check:
    checkdict[c[:][1][-1]].append(c[:][1][-2])
for (key, value) in list(checkdict.items()):
    (hist, bins) = np.histogram(value, bins=np.arange(1, int(classes)
                                + 2))  # , normed=True)
    checktest.append(hist)
    labelsY.append(int(key))

fig = plt.figure(num=None, dpi=80, facecolor='white', figsize=(9, 11))
ax = fig.add_subplot(111)
plt.title('Class assignment of each particle - last iteration',
          fontsize=16, fontweight='bold')
plt.xlabel('Class assignment iteration %s' % (int(iterations) - 2),
           fontsize=13)
plt.ylabel('Class assignment iteration %s' % (int(iterations) - 1),
           fontsize=13)
plt.grid()
ticks = np.arange(0, int(classes) + 1)
labelsX = np.arange(1, int(classes) + 2)
plt.xticks(ticks, labelsX)
# plt.yticks(ticks, labelsY)
plt.imshow(checktest, aspect='auto', interpolation='nearest',
           origin='lower')  # FIXME values in box
cb2 = plt.colorbar(ticks=np.arange(0, 1, 0.1))
cb2.set_label('Fraction of particles went into group #')
plt.figtext(0, 0,
            'Particle class assignment changes in the last iteration')
pdf.savefig()
plt.close()
# plt.show()

######################################################################################################################
###########################################################################
#### Class assignments per micrograph of the last iteration

### Convert dictionary to 2D array

micval = []
micfirst = []
micticks = []
for mics in micnumdict:
    (max_index, max_value) = max(enumerate(micnumdict[mics][-1]),
                                 key=operator.itemgetter(1))
    micval.append(micnumdict[mics][-1])
    micticks.append(mics)

### Sort last iter by most prominent group

micval = np.array(micval)

### Plot heat map last iteration

cmap = plt.get_cmap('jet', np.max(micval) - np.min(micval) + 1)
plt.figure(num=None, dpi=80, facecolor='white', figsize=(9, 11))
plt.title('Class assignments of each micrograph - last iteration',
          fontsize=16, fontweight='bold')
plt.xlabel('Class #', fontsize=13)
plt.ylabel('Micrograph #', fontsize=13)
ticks = np.arange(0, int(classes) + 1)
labels = np.arange(1, int(classes) + 2)
plt.xticks(ticks, labels)
plt.grid()
plt.imshow(micval, aspect='auto', interpolation='nearest', cmap=cmap)
cb3 = plt.colorbar()
cb3.set_label('Total number of particles in class #')
plt.figtext(0, 0,
            'Micrographs contributing to certain classes (e.g. important when merging datasets)'
            )
pdf.savefig()
plt.close()
# plt.show()

###### Find out how often particles are jumping

scorelist = []
for g in groupnumarraysorted[:, 2:]:

    count = 1
    stayed = 0
    g = g[::-1]  # # reverse order
    for (i, ig) in enumerate(g):
        if i > 0 and count > 0:

            if len(set(g)) == 1:
                stayed = 0
                count = 0
            if int(ig) == int(g[i - 1]):
                count += 1
            if int(ig) != int(g[i - 1]):
                stayed = count
                count = 0

    g = g.tolist()
    counter = [[x, g.count(x)] for x in set(g)]

    # score3 = (iterations-2)/float(len(counter))

    score2 = float(len(counter))
    score1 = stayed / score2 / (iterations - 2)
    scorelist.append(score1)

plt.figure(num=None, dpi=80, facecolor='white', figsize=(9, 11))
plt.title('Particle jump score', fontsize=16, fontweight='bold')
plt.xlabel('Score = (# class assignments/# of iterations)', fontsize=13)
plt.ylabel('# of particles with score normalized', fontsize=13)
plt.grid()

# histbins = sorted(set(scorelist))

histbins = np.arange(0, .5, 0.05)
plt.hist(scorelist, bins=histbins, density=True)
mean1 = np.mean(scorelist)
variance1 = np.var(scorelist)
sigma1 = np.sqrt(variance1)

# x1 = np.linspace(min(scorelist), max(scorelist), 100)

x1 = np.linspace(0, .5, 100)
plt.figtext(0, 0, 'Gaussian sigma: %s, variance: %s, mean: %s'
            % (sigma1, variance1, mean1))
plt.plot(x1, scipy.stats.norm.pdf(x1, mean1, sigma1))
pdf.savefig()
plt.close()
### Number of assignment changes per iteration

plt.figure(num=None, dpi=80, facecolor='white', figsize=(9, 11))
plt.title('Total assignment changes per iteration', fontsize=16,
          fontweight='bold')
plt.xlabel('Iteration #', fontsize=13)
plt.ylabel('# of changed assignments', fontsize=13)
plt.grid()
plt.plot(changes[2:])
pdf.savefig()
plt.close()
# plt.show()

#########################################################################################################################################################

checkarraysorted = checkarray[sortindices]

### Plot histogram of each column in data.star sorted by the class assignments of the last iteration

histocolarray = collections.defaultdict(list)

for column in checkarraysorted.transpose():
    lastitgrouparray = collections.defaultdict(list)
    for (ci, col) in enumerate(column):
        lastitgrouparray[groupnumarraysorted[:, -1][ci] - 1].append(col)

        # print groupnumarraysorted[:,-1][ci]-1, ci, col

    for (key, value) in list(lastitgrouparray.items()):
        histocolarray[key].append(value)

        # print key....#STILL has 0,3,4!!!

fincolarray = collections.defaultdict(list)

for (key, value) in list(histocolarray.items()):
    for (vi, v) in enumerate(value):

        # print key, vi, v....................########### LOSING KEY

        fincolarray[vi].append(v)

for (key2, value2) in list(fincolarray.items()):
    rangeval = []  # temp2 = [];
    temp = []

    for (vi2, v2) in enumerate(value2):
        for v3 in v2:
            rangeval.append(v3)

    if len(set(rangeval)) > 1:
        for (vi2, v2) in enumerate(value2):
            temp.append(v2)

    if len(temp) > 0:
        plt.figure(num=None, dpi=120, facecolor='white', figsize=(9, 11))

        if plottype == 'bar':
            (n, bins, patches) = plt.hist(temp, range=(min(rangeval),
                    max(rangeval)), histtype='bar')
            cmap = plt.get_cmap('jet', classes + 1)
            colors = np.arange(0, classes + 1)
            for (c, p) in zip(colors, patches):

                # d = c+1

                d = int(list(histocolarray.keys())[c]) + 1  # recolor based on class HACK
                plt.setp(p, 'facecolor', cmap(d))

            # plt.colorbar(ticks=np.arange(1, int(classes)+1))

        plt.title('Histogram Column %s' % checklistcol[int(key2)],
                  fontsize=16, fontweight='bold')
        plt.xlabel('%s' % checklistcol[int(key2)], fontsize=13)
        plt.ylabel('# particles per bin', fontsize=13)
        plt.grid()
        pdf.savefig()
        plt.close()

################################################################### DELETE UNWANTED PARTICLES FROM INITIAL STAR FILE ##FIXME

particcount = 0
initstarfile = '%s_it001_data.star' % rootname
a1 = open('%s_filtered.star' % rootname, 'w')
if filtstar != 'false' or micfilt != '':

 # #######################

    print('The mean jump score is: 					%s' % mean1)
    if sigmafac == 1:
        print('You did not specify a cutoff, I will use a sigma of 1 above mean: %s' \
            % (float(mean1) + float(sigmafac) * float(sigma1)))
    if sigmafac != 1:
        print('I will use a cutoff of: 					%s' % (float(mean1)
                + float(sigmafac) * float(sigma1)))
    cutoff = float(mean1) + float(sigmafac) * float(sigma1)

    for (i, score) in enumerate(scorelist):
        if score > cutoff:
            unwanted.append(i)

 # #######################

    with open(initstarfile, 'r') as g:
        for m in g:
            if '@' not in m:
                a1.write(m)
            if '@' in m:
                if micfilt != '':
                    if float(m.split()[rescol]) <= float(micfilt):
                        a1.write(m)

                # print('micfilt')

                if filtstar != 'false':
                    if particcount not in unwanted:
                        a1.write(m)

                # print('nomicfilt')

                if filtstar != 'false' and micfilt != '':
                    if float(m.split()[rescol]) <= float(micfilt) \
                        and particcount not in unwanted:
                        a1.write(m)

                # print('both')

                particcount += 1
    print('Saved %s_filtered.star file ommitting %s out of %s particles that changed classes too often' \
        % (rootname, len(unwanted), particcount))
a1.close()

print('Saved all plots in %s' % output)
pdf.close()
