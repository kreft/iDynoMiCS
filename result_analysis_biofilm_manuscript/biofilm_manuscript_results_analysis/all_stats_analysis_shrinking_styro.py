#do stats with ratio files
import csv
import scipy.stats as stats
import os
import ready_create_ratio_files_for_stats as cr
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

def open_file_get_prop(fn, m='b'):
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    if m == 'b': ind = 2
    elif m == 'p': ind = 3
    elif m == 'gr': ind = 4
    ratios = [0, 0, 0, 0]
    competition, ratio, number = [], [], [0, 0, 0, 0]
    for a in range(len(rows)):
        if a > 0:
            competition.append(rows[a][0])
            ratio.append(rows[a][ind])
    for b in range(len(competition)):
        name = competition[b]
        ratio[b] = float(ratio[b])
        if name[0] == '4':
            number[0] += 1
            if ratio[b] > 0:
                ratios[0] += 1
        elif name[0] == '8':
            number[1] += 1
            if ratio[b] > 0:
                ratios[1] += 1
        elif name[0] == '1':
            number[2] += 1
            if ratio[b] > 0:
                ratios[2] += 1
        elif name[0] == '3':
            number[3] += 1
            if ratio[b] > 0:
                ratios[3] += 1
    return ratios, number

def open_file_get_ratio(fn, count, all_ax, m='b', ):
    ax1, ax2, ax3, ax4 = all_ax[0], all_ax[1], all_ax[2], all_ax[3]
    sp, coun = ['', ''], 0
    a = 0
    while a < len(fn)-8:
        if fn[a] == 'c':
            a += 8
        if fn[a] == '_':
            coun += 1
            a += 1
        if coun < 2:
            sp[coun] += fn[a]
            a += 1
    for b in range(len(sp)):
        if sp[b] == 'ASNR':
            sp[b] = 'DS'
        elif sp[b] == 'SAR':
            sp[b] = 'AR'
        elif sp[b] == 'SFR':
            sp[b] = 'FR'
    name = sp[0]+'/'+sp[1]
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    a = 0
    ml = []
    if m == 'b': add = 1
    elif m == 'p': add = 2
    elif m == 'gr': add = 3
    while a < len(rows):
        if rows[a][0][0] == '4': ax = ax1
        elif rows[a][0][0] == '8': ax = ax2
        elif rows[a][0][0] == '1': ax = ax3
        elif rows[a][0][0] == '3': ax = ax4
        time, ratios = [], []
        ma = 0
        for b in range(len(rows[a])):
            if b > 1:
                time.append(float(rows[a][b]))
                if base_folder == '/Volumes/Robyn_W_2/july_2018/paper/styrofoam_competitions_new/' and sp == ['DS', 'FR']:
                    ratios.append(float(rows[a+add][b])*(-1))
                else:
                    ratios.append(float(rows[a+add][b]))
        ax.plot(time, ratios)
        a += 4
        if time[-1] > ma:
            ma = time[-1]+24
            ax.set_xlim([0, ma])
        if ma > 250:
            ml.append(MultipleLocator(100))
        else:
            ml.append(MultipleLocator(50))
    majorLocatory = MultipleLocator(1)
    for a in range(len(all_ax)):
        b = a
        a = all_ax[a]
        a.plot([0, ma], [0, 0], 'k--')
        if a != ax1:
            plt.setp(a.get_yticklabels(), visible=False)
        a.xaxis.set_major_locator(ml[b])
        a.yaxis.set_major_locator(majorLocatory)
    ax1.set_ylabel('log('+name+') biomass')
    return
    
def get_biomass_plots(save_stats, base_folder, m):
    files = os.listdir(save_stats+'time_course/')
    os.chdir(save_stats+'time_course/')
    new_files = []
    fig = plt.figure(figsize=(10, 12))
    all_ax = []
    ax1 = plt.subplot(6,4,1)
    ax2, ax3, ax4 = plt.subplot(6,4,2, sharey=ax1), plt.subplot(6,4,3, sharey=ax1), plt.subplot(6,4,4, sharey=ax1)
    ax1.set_title('4 cells'), ax2.set_title('8 cells'), ax3.set_title('16 cells'), ax4.set_title('32 cells')
    all_ax.append([ax1, ax2, ax3, ax4])
    ax1 = plt.subplot(6,4,5)
    ax2, ax3, ax4 = plt.subplot(6,4,6, sharey=ax1), plt.subplot(6,4,7, sharey=ax1), plt.subplot(6,4,8, sharey=ax1)
    all_ax.append([ax1, ax2, ax3, ax4])
    ax1 = plt.subplot(6,4,9)
    ax2, ax3, ax4 = plt.subplot(6,4,10, sharey=ax1), plt.subplot(6,4,11, sharey=ax1), plt.subplot(6,4,12, sharey=ax1)
    all_ax.append([ax1, ax2, ax3, ax4])
    ax1 = plt.subplot(6,4,13)
    ax2, ax3, ax4 = plt.subplot(6,4,14, sharey=ax1), plt.subplot(6,4,15, sharey=ax1), plt.subplot(6,4,16, sharey=ax1)
    all_ax.append([ax1, ax2, ax3, ax4])
    ax1 = plt.subplot(6,4,17)
    ax2, ax3, ax4 = plt.subplot(6,4,18, sharey=ax1), plt.subplot(6,4,19, sharey=ax1), plt.subplot(6,4,20, sharey=ax1)
    all_ax.append([ax1, ax2, ax3, ax4])
    ax1 = plt.subplot(6,4,21)
    ax2, ax3, ax4 = plt.subplot(6,4,22, sharey=ax1), plt.subplot(6,4,23, sharey=ax1), plt.subplot(6,4,24, sharey=ax1)
    all_ax.append([ax1, ax2, ax3, ax4])
    for a in range(len(files)):
        if files[a] != '.DS_Store' and files[a][-4:] != '.png':
            new_files.append(files[a])
    files = sorted(new_files)
    for a in range(len(files)):
        open_file_get_ratio(files[a], a, all_ax[a], m)
    #put this bit in if plotting proportional, and change below to all_ax[2] not all_ax[5]
    """
    for ax in range(len(all_ax)):
        if ax > 2:
            for a in all_ax[ax]:
                a.set_xticklabels(['']*10)
                a.set_yticklabels(['']*10)
                for spine in ['right', 'left', 'top', 'bottom']:
                    a.spines[spine].set_color('none')
                a.tick_params(top="off", bottom="off", right="off", left="off")
    """
    for a in all_ax[5]:
        a.set_xlabel('Time (h)')
    fig.subplots_adjust(wspace=0.15, hspace=0.4)
    plt.savefig('All_ratios '+m+'.png', bbox_inches='tight', dpi=600)
    plt.close()
    return

def get_stats(save_stats, base_folder, m):
    files = os.listdir(save_stats+'end_only/')
    os.chdir(save_stats+'end_only/')
    names, all_ratios, all_stats = [], [], []
    for a in files:
        if a != '.DS_Store':
            ratios, number = open_file_get_prop(a, m)
            names.append(a[:-4])
            all_ratios.append(ratios)
            stats_p = []
            for b in range(len(ratios)):
                stats_p.append(stats.binom_test(ratios[b], number[b], 0.5))
                ratios[b] = str(ratios[b])+' of '+str(number[b])
            all_stats.append(stats_p)
    top_line = ['Competition', '4 cells', '8 cells', '16 cells', '32 cells', '4 cells stats', '8 cells stats', '16 cells stats', '32 cells stats']
    os.chdir(save_stats)
    with open('Competition winners '+m+'.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(top_line)
        for a in range(len(names)):
            writer.writerow([names[a]]+all_ratios[a]+all_stats[a])
    return


#run with proportional, shrinking and styrofoam to have all files
#change save_stats to where the stats should be saved to, and base_folder to where the simulations are


#proportional
#save_stats = '/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/stats/proportional/'
#base_folder = '/Volumes/Robyn_W_2/july_2018/paper/biofilms_proportional/'

#shrinking
save_stats = '/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/stats/shrinking/'
base_folder = '/Volumes/Robyn_W_2/july_2018/paper/shrinking_biofilms/'

#styrofoam
#save_stats = '/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/stats/styrofoam/'
#base_folder = '/Volumes/Robyn_W_2/july_2018/paper/styrofoam_competitions_new/'


cr.run_functions(save_stats, base_folder, end_only=False)
cr.run_functions(save_stats, base_folder, end_only=True)
ms = ['b', 'p', 'gr']
#ms = ['b']
for m in ms:
    get_stats(save_stats, base_folder, m)
    get_biomass_plots(save_stats, base_folder, m)