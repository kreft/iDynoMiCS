import numpy
import matplotlib.pyplot as plt
import os
import csv

path = '/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/stats/proportional/time_course/'
os.chdir(path)
AR_DS, AR_AR, DS_DS = 'SAR_ASNR_all.csv', 'control_SAR_SAR_all.csv', 'control_ASNR_ASNR_all.csv'

#Biomass
sn = 'Biomass_time_courses.png'
adding = 1

#Population number
#sn = 'Population_size_time_courses.png'
#adding = 2


fig = plt.figure(figsize=(8.27, 5.51))
ax1 = plt.subplot2grid((10, 4), (0, 0), rowspan=4)
ax2, ax3, ax4 = plt.subplot2grid((10, 4), (0, 1), sharey=ax1, rowspan=4), plt.subplot2grid((10, 4), (0, 2), sharey=ax1, rowspan=4), plt.subplot2grid((10, 4), (0, 3), sharey=ax1, rowspan=4)
ax5 = plt.subplot2grid((10, 10), (6, 1), colspan=2, rowspan=4)
ax6, ax7 = plt.subplot2grid((10, 10), (6, 4), sharey=ax5, colspan=2, rowspan=4), plt.subplot2grid((10, 10), (6, 7), sharey=ax5, colspan=2, rowspan=4)
all_ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]
t_ax = [ax1, ax2, ax3, ax4]
id_ax = [ax5, ax6, ax7]
remy = [ax2, ax3, ax4]
for a in t_ax:
    a.set_xlabel('Time (h)')
for a in id_ax:
    a.set_xlabel('Initial cell density')
for b in remy:
    plt.setp(b.get_yticklabels(), visible=False)

ax1.set_title('4 cells'), ax2.set_title('8 cells'), ax3.set_title('16 cells'), ax4.set_title('32 cells')
ax5.set_title('AR/DS'), ax6.set_title('AR/AR'), ax7.set_title('DS/DS')
ax1.set_ylabel('log(AR/DS) biomass')
ax5.set_ylabel('log(AR/DS) biomass')
ax6.set_ylabel('log(AR/AR) biomass')
ax7.set_ylabel('log(DS/DS) biomass')

with open(AR_DS, 'rU') as f:
    AR_DS = []
    for row in csv.reader(f):
        AR_DS.append(row)

with open(AR_AR, 'rU') as f:
    AR_AR = []
    for row in csv.reader(f):
        AR_AR.append(row)

with open(DS_DS, 'rU') as f:
    DS_DS = []
    for row in csv.reader(f):
        DS_DS.append(row)

t250 = [[[], [], [], []], [[], [], [], []], [[], [], [], []]]
t250_means = [[], [], []]
t250_sd = [[], [], []]
files = [AR_DS, AR_AR, DS_DS]
for a in range(len(files)):
    rows = files[a]
    b = 0
    num, ind = ['4', '8', '16', '32'], [0, 1, 2, 3]
    while b < 800:
        prop = float(rows[b+adding][len(rows[b+adding])-1])
        if rows[b][0][1] == 'c':
            rb0 = rows[b][0][0]
        else:
            rb0 = rows[b][0][:2]
        for d in range(len(num)):
            if num[d] == rb0:
                i = ind[d]
        t250[a][i].append(prop)
        b += 4

for c in range(len(t250)):
    for d in range(len(t250[c])):
        t250_means[c].append(numpy.mean(t250[c][d]))
        t250_sd[c].append(numpy.std(t250[c][d]))

x = [4, 8, 16, 32]
mean_plots = [ax5, ax6, ax7]
for e in range(len(t250_means)):
    for f in range(len(t250_means[e])):
        mean_plots[e].errorbar(x[f], t250_means[e][f], yerr=t250_sd[e][f], marker='o', capsize=2)
    mean_plots[e].plot([0, 35], [0, 0], 'k--')
    mean_plots[e].set_xlim([0, 35])
    plt.setp(mean_plots[e], xticks=[4, 8, 16, 32], xticklabels=num)

times, biomass = [[], [], [], []], [[], [], [], []]
h = 0
while h < 800:
    this_time, this_biom = [], []
    for i in range(len(AR_DS[h])):
        if i > 1:
            this_time.append(float(AR_DS[h][i]))
            this_biom.append(float(AR_DS[h+adding][i]))
    if AR_DS[h][0][1] == 'c':
        rb0 = AR_DS[h][0][0]
    else:
        rb0 = AR_DS[h][0][:2]
    for j in range(len(num)):
            if num[j] == rb0:
                ni = ind[j]
    times[ni].append(this_time)
    biomass[ni].append(this_biom)
    h += 4
biom_plots = [ax1, ax2, ax3, ax4]
for k in range(len(times)):
    for l in range(len(times[k])):
        biom_plots[k].plot(times[k][l], biomass[k][l], lw=0.8)
    biom_plots[k].plot([0, 300], [0, 0], 'k--')
    plt.setp(biom_plots[k], xticks=[0, 100, 200, 300], xticklabels=[0, 100, 200, 300])
    biom_plots[k].set_xlim([0, 300])

#plt.tight_layout()
os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/proportional_biofilms')
plt.savefig(sn, bbox_inches='tight', dpi=600)