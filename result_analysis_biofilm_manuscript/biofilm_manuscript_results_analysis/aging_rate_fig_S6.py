from __future__ import division
from __future__ import with_statement
import matplotlib.pyplot as plt
import numpy
import os
import toolbox_idynomics
import geoMean as gm

#main_folder = '/Volumes/Robyn_W_2/july_2018/paper/aging_rate/tests_3/'
main_folder = 'F:\\paper\\aging_rate\\tests_3'
save_folder = os.getcwd()
os.chdir(main_folder)

sar = [['dynamic_SAR_01(20180702_1051)', 'dynamic_proportional_SAR_0000(20180630_1512)', 'dynamic_proportional_SAR_0500(20180630_1754)', 'dynamic_proportional_SAR_1000(20180630_2018)', 'dynamic_proportional_SAR_1500(20180630_2225)', 'dynamic_proportional_SAR_2000(20180701_0013)', 'dynamic_proportional_SAR_2100(20180701_0141)', 'dynamic_proportional_SAR_2200(20180701_0302)', 'dynamic_proportional_SAR_2300(20180701_0416)', 'dynamic_proportional_SAR_2400(20180701_0530)', 'dynamic_proportional_SAR_2500(20180701_1040)', 'dynamic_proportional_SAR_3000(20180701_1134)'],
       ['dynamic_SAR_01(20180702_1054)', 'dynamic_proportional_SAR_0000(20180630_1519)', 'dynamic_proportional_SAR_0500(20180630_1756)', 'dynamic_proportional_SAR_1000(20180630_2021)', 'dynamic_proportional_SAR_1500(20180630_2226)', 'dynamic_proportional_SAR_2000(20180701_0014)', 'dynamic_proportional_SAR_2100(20180701_0142)', 'dynamic_proportional_SAR_2200(20180701_0303)', 'dynamic_proportional_SAR_2300(20180701_0419)', 'dynamic_proportional_SAR_2400(20180701_0531)', 'dynamic_proportional_SAR_2500(20180706_1226)', 'dynamic_proportional_SAR_3000(20180701_1135)'], 
        ['dynamic_SAR_01(20180702_1820)', 'dynamic_proportional_SAR_0000(20180708_2123)', 'dynamic_proportional_SAR_0500(20180708_2358)', 'dynamic_proportional_SAR_1000(20180709_0216)', 'dynamic_proportional_SAR_1500(20180709_0727)', 'dynamic_proportional_SAR_2000(20180709_0914)', 'dynamic_proportional_SAR_2100(20180709_1142)', 'dynamic_proportional_SAR_2200(20180709_1315)', 'dynamic_proportional_SAR_2300(20180709_1429)', 'dynamic_proportional_SAR_2400(20180709_1603)', 'dynamic_proportional_SAR_2500(20180709_1710)', 'dynamic_proportional_SAR_3000(20180709_2208)']]
sfr = [['dynamic_SFR_01(20180702_1227)', 'dynamic_proportional_SFR_0000(20180701_1144)', 'dynamic_proportional_SFR_0500(20180701_1414)', 'dynamic_proportional_SFR_1000(20180701_1635)', 'dynamic_proportional_SFR_1500(20180701_2052)', 'dynamic_proportional_SFR_2000(20180701_2243)', 'dynamic_proportional_SFR_2100(20180702_0019)', 'dynamic_proportional_SFR_2200(20180702_0148)', 'dynamic_proportional_SFR_2300(20180702_0303)', 'dynamic_proportional_SFR_2400(20180702_0411)', 'dynamic_proportional_SFR_2500(20180702_0511)', 'dynamic_proportional_SFR_3000(20180702_0604)'], 
       ['dynamic_SFR_01(20180702_1230)', 'dynamic_proportional_SFR_0000(20180701_1145)', 'dynamic_proportional_SFR_0500(20180701_1419)', 'dynamic_proportional_SFR_1000(20180706_1328)', 'dynamic_proportional_SFR_1500(20180701_2053)', 'dynamic_proportional_SFR_2000(20180701_2245)', 'dynamic_proportional_SFR_2100(20180702_0027)', 'dynamic_proportional_SFR_2200(20180702_0158)', 'dynamic_proportional_SFR_2300(20180702_0313)', 'dynamic_proportional_SFR_2400(20180702_0420)', 'dynamic_proportional_SFR_2500(20180702_0520)', 'dynamic_proportional_SFR_3000(20180702_0610)'], 
        ['dynamic_SFR_01(20180703_1034)', 'dynamic_proportional_SFR_0000(20180710_0903)', 'dynamic_proportional_SFR_0500(20180710_1152)', 'dynamic_proportional_SFR_1000(20180710_1619)', 'dynamic_proportional_SFR_1500(20180710_1909)', 'dynamic_proportional_SFR_2000(20180711_0351)', 'dynamic_proportional_SFR_2100(20180711_0524)', 'dynamic_proportional_SFR_2200(20180711_0647)', 'dynamic_proportional_SFR_2300(20180711_0803)', 'dynamic_proportional_SFR_2400(20180711_0912)', 'dynamic_proportional_SFR_2500(20180711_1012)', 'dynamic_proportional_SFR_3000(20180711_1101)']]
asnr = [['dynamic_ASNR_01(20180702_0803)', 'dynamic_proportional_ASNR_0000(20180630_0111)', 'dynamic_proportional_ASNR_0500(20180630_0353)', 'dynamic_proportional_ASNR_1000(20180630_0626)', 'dynamic_proportional_ASNR_1500(20180630_0832)', 'dynamic_proportional_ASNR_2000(20180630_1013)', 'dynamic_proportional_ASNR_2100(20180630_1120)', 'dynamic_proportional_ASNR_2200(20180630_1221)', 'dynamic_proportional_ASNR_2300(20180630_1315)', 'dynamic_proportional_ASNR_2400(20180630_1400)', 'dynamic_proportional_ASNR_2500(20180630_1437)', 'dynamic_proportional_ASNR_3000(20180630_1505)'],
        ['dynamic_ASNR_01(20180702_1652)', 'dynamic_proportional_ASNR_0000(20180706_0920)', 'dynamic_proportional_ASNR_0500(20180630_0357)', 'dynamic_proportional_ASNR_1000(20180630_0632)', 'dynamic_proportional_ASNR_1500(20180630_0836)', 'dynamic_proportional_ASNR_2000(20180630_1017)', 'dynamic_proportional_ASNR_2100(20180630_1125)', 'dynamic_proportional_ASNR_2200(20180630_1227)', 'dynamic_proportional_ASNR_2300(20180630_1320)', 'dynamic_proportional_ASNR_2400(20180630_1406)', 'dynamic_proportional_ASNR_2500(20180630_1443)', 'dynamic_proportional_ASNR_3000(20180630_1512)'], 
        ['dynamic_ASNR_01(20180702_1653)', 'dynamic_proportional_ASNR_0000(20180707_2233)', 'dynamic_proportional_ASNR_0500(20180708_0911)', 'dynamic_proportional_ASNR_1000(20180708_1202)', 'dynamic_proportional_ASNR_1500(20180708_1420)', 'dynamic_proportional_ASNR_2000(20180708_1603)', 'dynamic_proportional_ASNR_2100(20180708_1714)', 'dynamic_proportional_ASNR_2200(20180708_1821)', 'dynamic_proportional_ASNR_2300(20180708_1919)', 'dynamic_proportional_ASNR_2400(20180708_2006)', 'dynamic_proportional_ASNR_2500(20180708_2045)', 'dynamic_proportional_ASNR_3000(20180708_2115)']]
labels = ['Not proportional\n0.1', '0.0', '0.05', '0.1', '0.15', '0.2', '0.21', '0.22', '0.23', '0.24', '0.25', '0.3']
colors = ['k', '#9900FF', '#0000FF', '#8899FF', '#99FFFF', '#99FFCC', '#00CC00', '#FF0000', '#009966', '#336600', '#999900', '#FFFF00']
l1, l2, l3 = '-', '--', '--'
linestyles = [l1, l3, l3, l3, l3, l3, l3, l2, l3, l3, l3, l3]
alpha_rest = 0.8
alpha = [1, alpha_rest, alpha_rest, alpha_rest, alpha_rest, 1, alpha_rest, 1, alpha_rest, alpha_rest, alpha_rest, alpha_rest, alpha_rest]

ax1 = plt.subplot(331)
ax2, ax3, ax4 = plt.subplot(332, sharey=ax1), plt.subplot(333, sharey=ax1), plt.subplot(334)
ax5, ax6, ax7 = plt.subplot(335, sharey=ax4), plt.subplot(336, sharey=ax4), plt.subplot(337)
ax8, ax9 = plt.subplot(338, sharey=ax7), plt.subplot(339, sharey=ax7)
ax1.set_title('Adaptive repair'), ax2.set_title('Fixed repair'), ax3.set_title('Asymmetric \nsegregation')
ax1.set_ylabel('Growth rate'), ax4.set_ylabel('Population size'), ax7.set_ylabel('Population age')
plt.setp(ax1.get_xaxis(), visible=False), plt.setp(ax2.get_xaxis(), visible=False), plt.setp(ax3.get_xaxis(), visible=False), plt.setp(ax4.get_xaxis(), visible=False), plt.setp(ax5.get_xaxis(), visible=False), plt.setp(ax6.get_xaxis(), visible=False)
plt.setp(ax2.get_yaxis(), visible=False), plt.setp(ax3.get_yaxis(), visible=False), plt.setp(ax5.get_yaxis(), visible=False), plt.setp(ax6.get_yaxis(), visible=False), plt.setp(ax8.get_yaxis(), visible=False), plt.setp(ax9.get_yaxis(), visible=False)
ax8.set_xlabel('Time (days)')

def get_params(agent_output):
    biomass_names=['activeBiomassGrowth', 'activeBiomassRepair', 'inactiveBiomassGrowth', 'inactiveBiomassRepair']
    growth, ages, total_biomass = [], [], []
    pop = 0
    for species in agent_output.species_outputs:
        for cell in species.members:
            pop += 1
            growth.append(float(cell.get_specific_growth_rate(biomass_names)))
            ages.append(float(cell.vars['age']))
            ar = float(cell.vars['activeBiomassRepair'])
            ir = float(cell.vars['inactiveBiomassRepair'])
            ag = float(cell.vars['activeBiomassGrowth'])
            ig = float(cell.vars['inactiveBiomassGrowth'])
            total_biomass.append(ar+ir+ag+ig)
    if max(growth) > 0:
        growthrate, delta = gm.geomMeanExtension(growth)
    else:
        growthrate = 0
    sd = numpy.std(growth)
    age = numpy.mean(ages)
    age_sd = numpy.std(ages)
    biom = sum(total_biomass)
    return growthrate, sd, pop, age, age_sd, biom

iters = [0, 240, 480, 720, 960, 1200, 1440, 1680, 1920, 2160, 2400]
x = []
for z in range(len(iters)):
    x.append(iters[z]/240)

for a in range(3):
    os.chdir(main_folder)
    fig = plt.figurefigsize=(9, 9)
    ax1 = plt.subplot(331)
    ax2, ax3, ax4 = plt.subplot(332, sharey=ax1), plt.subplot(333, sharey=ax1), plt.subplot(334)
    ax5, ax6, ax7 = plt.subplot(335, sharey=ax4), plt.subplot(336, sharey=ax4), plt.subplot(337)
    ax8, ax9 = plt.subplot(338, sharey=ax7), plt.subplot(339, sharey=ax7)
    ax1.set_title('Adaptive repair'), ax2.set_title('Fixed repair'), ax3.set_title('Asymmetric \nsegregation')
    ax1.set_ylabel('Growth rate'), ax4.set_ylabel('Population size'), ax7.set_ylabel('Population age')
    plt.setp(ax1.get_xaxis(), visible=False), plt.setp(ax2.get_xaxis(), visible=False), plt.setp(ax3.get_xaxis(), visible=False), plt.setp(ax4.get_xaxis(), visible=False), plt.setp(ax5.get_xaxis(), visible=False), plt.setp(ax6.get_xaxis(), visible=False)
    plt.setp(ax2.get_yaxis(), visible=False), plt.setp(ax3.get_yaxis(), visible=False), plt.setp(ax5.get_yaxis(), visible=False), plt.setp(ax6.get_yaxis(), visible=False), plt.setp(ax8.get_yaxis(), visible=False), plt.setp(ax9.get_yaxis(), visible=False)
    ax8.set_xlabel('Time (days)')
    for b in range(3):
        if b == 0:
            ax = [ax1, ax4, ax7]
            files = sar
        elif b == 1:
            ax = [ax2, ax5, ax8]
            files = sfr
        elif b == 2:
            ax = [ax3, ax6, ax9]
            files = asnr
        for c in range(12):
            if c != 2 and c != 4 and c != 11:
                f = files[a][c]
                sim = toolbox_idynomics.SimulationDirectory(f)
                params = [[], [], [], [], [], []]
                for d in range(len(iters)):
                    iter_info = sim.get_single_iterate(iters[d])
                    gr, gr_sd, p, ag, a_sd, biom = get_params(iter_info.agent_output)
                    params[0].append(gr), params[1].append(gr_sd), params[2].append(p), params[3].append(ag), params[4].append(a_sd), params[5].append(biom)
                sim.clean_up()
                ax[0].errorbar(x, params[0], linestyle=linestyles[c], linewidth=1, color=colors[c], label=labels[c], markersize=3, alpha=alpha[c])
                ax[1].errorbar(x, params[2], linestyle=linestyles[c], color=colors[c], linewidth=1, markersize=3, alpha=alpha[c])
                ax[2].errorbar(x, params[3], linestyle=linestyles[c], linewidth=1, color=colors[c], markersize=3, alpha=alpha[c])
                if c == 7:
                    all_params = params
        if b == 2:
            ax[0].legend(bbox_to_anchor=(2.5, 1.1))
        ax[0].errorbar(x, all_params[0], linestyle=linestyles[7], linewidth=1, color=colors[7], markersize=3, alpha=alpha[7])
        ax[1].errorbar(x, all_params[2], linestyle=linestyles[7], linewidth=1, color=colors[7], markersize=3, alpha=alpha[7])
        ax[2].errorbar(x, all_params[3], linestyle=linestyles[7], linewidth=1, color=colors[7], markersize=3, alpha=alpha[7])
    plt.sca(ax1)
    plt.yticks([0, 0.1, 0.2, 0.3])
    plt.sca(ax4)
    plt.yticks([800, 1200, 1600, 2000, 2400])
    plt.sca(ax7)
    plt.yticks([0, 0.1, 0.2, 0.3, 0.4])
    #plt.tight_layout()
    print('Saving figure ', str(a))
    os.chdir(save_folder)
    plt.savefig('Replicate '+str(a)+'.png', bbox_inches='tight', dpi=600)
    plt.close()
            
