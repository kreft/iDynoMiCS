axisB = fig.add_subplot('B', 222)
paths = ['an10_const_steady_results.xml', 'ao10_const_steady_results.xml', sys.argv[1]]
attributes = {'name':'specific growth rate population structure',
                    'starting_time':'2400', 'header':'bin,frequency'}

xlim, ylim = 0.082, 0.6

for i in range(3):
    results_path = os.path.join(input_path, paths[i])
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        bin_edges = result.vars['bin'].split('-')
        result.vars['mid_point'] = (float(bin_edges[0])+float(bin_edges[1]))/2
    result_set.members.sort(key=lambda r: r.vars['mid_point'])
    x = [float(r.vars['frequency']) for r in result_set.members]
    y = [r.vars['mid_point'] for r in result_set.members]
    axisB.plot(x, y, color=colors[i])
axisB.text(0.021, 0.575,'0', color=colors[0], fontsize=fs)
axisB.text(0.055, 0.36, '1', color=colors[0], fontsize=fs)
axisB.text(0.019, 0.245, '2', color=colors[0], fontsize=fs)
axisB.text(0.009, 0.08, '3', color=colors[0], fontsize=fs)
axisB.text(0.009, 0.0, '4', color=colors[0], fontsize=fs)
axisB.text(0.03, 0.515,'0', color=colors[1], fontsize=fs)
axisB.text(0.075, 0.38, '1', color=colors[1], fontsize=fs)
axisB.text(0.012, 0.35, '2', color=colors[1], fontsize=fs)
axisB.text(0.0155, 0.31, '3', color=colors[1], fontsize=fs)
axisB.text(0.006, 0.24, '4', color=colors[1], fontsize=fs)
axisB.text(0.003, 0.16, '5', color=colors[1], fontsize=fs)
aging_extras.draw_cell(axisB, 0.6*xlim, 0.7*xlim, 0.9*ylim, 0.05*xlim,
                               y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisB, 0.8*xlim, 0.8*ylim, xlim, ylim,
                                                            y2xscale=ylim/xlim)
axisB.set_xlim([0.00, xlim])
axisB.set_xticks([0.00, 0.02, 0.04, 0.06, 0.08])
axisB.set_ylim([0, ylim])
axisB.set_xlabel('Frequency in population')
axisB.set_title('Growth Rate Distribution')
setp( axisB.get_yticklabels(), visible=False)
