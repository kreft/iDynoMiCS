#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import matplotlib
from matplotlib import rcParams
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
#from PIL import Image
#import Image
from pylab import *
import matplotlib.cbook as cbook
import numpy
import os
import random
from scipy.spatial import KDTree
from scipy.stats import linregress
import toolbox_basic


def mm2inch(mm):
    return mm*0.0393700787


class JournalFigure:
    def __init__(self):
        self.subplots = {}
    def save(self, path, dpi=300, clear=True):
        print('Saving as '+path)
        matplotlib.pyplot.savefig(path, dpi=dpi)
        if clear:
            matplotlib.pyplot.clf()
    def set_font_size(self, fontsize):
        matplotlib.rc('font', **{'size':fontsize})
    def bottom_axis_only(self, axis):
        axis.axes.get_yaxis().set_visible(False)
        axis.spines['right'].set_color('none')
        axis.spines['left'].set_color('none')
        axis.spines['top'].set_color('none')
    def add_subplot(self, label, position,
                    frameon=True, aspect='auto', axisbg='w', projection=None, sharey=True):
        axis = self.fig.add_subplot(position, aspect=aspect,
                                     axisbg=axisbg, projection=projection)
        # Moved frameon here as the keyword was causing problems in 3D plots
        if not frameon:
            axis.axis('off')
        setp( axis.get_xticklabels(), visible=frameon)
        setp( axis.get_yticklabels(), visible=frameon)
        if not frameon:
            axis.tick_params(top="off", bottom="off", right="off", left="off")
        self.subplots[axis] = label
        return axis
    def find_axis_from_label(self, label):
        for ax, lab in self.subplots.iteritems():
            if lab == label:
                return ax
        toolbox_basic.error_message('Could not find subplot with',
                                                        'label '+str(label))
        return None
    def use_image(self, axis, image_path):
        image_path = toolbox_basic.check_path(image_path)
        datafile = cbook.get_sample_data(image_path)
        image = Image.open(datafile)
        axis.set_ylim(0,image.size[0])
        axis.set_ylim(image.size[1],0)
        return imshow(image)
    def process_subplots(self):
        pass
    def subplots_adjust(self, left=None, right=None, top=None, bottom=None,
                                   hspace=None, wspace=None):
        self.fig.subplots_adjust(left=left, right=right, top=top,
                                  bottom=bottom, hspace=hspace, wspace=wspace)
    def inset_axes(self, padding=0.02):
        for axis, label in self.subplots.iteritems():
            xlim = axis.get_xlim()
            xdiff = padding*(xlim[1] - xlim[0])
            ylim = axis.get_ylim()
            ydiff = padding*(ylim[1] - ylim[0])
            axis.set_xlim([xlim[0]-xdiff, xlim[1]+xdiff])
            axis.set_ylim([ylim[0]-ydiff, ylim[1]+ydiff])
    def process_lines(self):
        for axis, label in self.subplots.iteritems():
            for line in axis.spines.itervalues():
                line.set_linewidth(1)
            axis.tick_params(top='off', right='off')
            for line in axis.get_xticklines() + axis.get_yticklines():
                line.set_markeredgewidth(1)
    def set_yaxis_label_positions(self, x, y, labels=[]):
        if labels == []:
            labels = self.subplots.values()
        for (axis, label) in self.subplots.iteritems():
            if label in labels:
                axis.yaxis.set_label_coords(x, y, transform=axis.transAxes)


class BmcFigure(JournalFigure):
    def __init__(self, double_column=False, height=0.0):
        JournalFigure.__init__(self)
        print("Using figure dimensions for BMC:")
        # width and height
        single, double, maximum = mm2inch(85), mm2inch(170), mm2inch(225)
        if double_column: self.width = double
        else:             self.width = single
        if height == 0.0: self.height = self.width
        elif height == 'single': self.height = single
        elif height == 'double': self.height = double
        elif height == 'max':    self.height = maximum
        else:                    self.height = min(height, maximum)
        print("width x height = "+str(self.width)+" x "+str(self.height))
        # default settings
        matplotlib.rc('font',**{'family':'Arial','weight':'normal','size':10})
        matplotlib.rc('mathtext', fontset='stixsans', default='regular')
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        # the figure!
        self.fig = matplotlib.pyplot.figure(figsize=(self.width, self.height),
                                                             facecolor='white')
    def process_subplots(self, label_pos=(1, 3), padding=0.1):
        for axis, label in self.subplots.iteritems():
            xlim = axis.get_xlim()
            xdiff = 0.02*(xlim[1] - xlim[0])
            ylim = axis.get_ylim()
            ydiff = 0.02*(ylim[1] - ylim[0])
            axis.set_xlim([xlim[0]-xdiff, xlim[1]+xdiff])
            axis.set_ylim([ylim[0]-ydiff, ylim[1]+ydiff])
            axis.text(label_pos[0], label_pos[1], label,
                         transform=axis.transAxes,
                         va='top', fontsize=10, fontweight='normal')
            for line in axis.spines.itervalues():
                line.set_linewidth(1)
            axis.tick_params(bottom='on', top='off', left='on', right='off')
            for line in axis.get_xticklines() + axis.get_yticklines():
                line.set_markeredgewidth(1)


class PlosFigure(JournalFigure):
    def __init__(self, double_column=False, height=0.0):
        JournalFigure.__init__(self)
        print("Using figure dimensions for PLoS:")
        # width and height
        single, double, maximum = 3.27, 6.83, 9.19
        if double_column: self.width = double
        else:             self.width = single
        if height == 0.0: self.height = self.width
        elif height == 'single': self.height = single
        elif height == 'double': self.height = double
        elif height == 'max': self.height = maximum
        else:             self.height = min(height, maximum)
        print("width x height = "+str(self.width)+" x "+str(self.height))
        # default settings
        matplotlib.rc('font',**{'family':'Arial','weight':'normal','size':10})
        matplotlib.rc('mathtext', fontset='stixsans', default='regular')
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        # the figure!
        self.fig = matplotlib.pyplot.figure(figsize=(self.width, self.height),
                                                             facecolor='white')
    def add_labels(self, label_pos=(-0.18, 1.01)):
        for axis, label in self.subplots.iteritems():
            axis.text(label_pos[0], label_pos[1], label, fontsize=12, va='top',
                                   transform=axis.transAxes, fontweight='bold')
    def process_subplots(self, label_pos=(-0.18, 1.01), padding=0.02):
        self.inset_axes(padding=padding)
        self.add_labels(label_pos=label_pos)
        self.process_lines()


class NatCommsFigure(JournalFigure):
    def __init__(self, double_column=False, height=0.0):
        JournalFigure.__init__(self)
        print("Using figure dimensions for Nature Communications:")
        # width and height
        single, double, maximum = mm2inch(89), mm2inch(183), mm2inch(247)
        if double_column: self.width = double
        else:             self.width = single
        if height == 0.0: self.height = self.width
        elif height == 'single': self.height = single
        elif height == 'double': self.height = double
        elif height == 'max': self.height = maximum
        else:             self.height = min(height, maximum)
        print("width x height = "+str(self.width)+" x "+str(self.height))
        # default settings
        matplotlib.rc('font',**{'family':'Arial','weight':'normal','size':10})
        matplotlib.rc('mathtext', fontset='stixsans', default='regular')
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        # the figure!
        self.fig = matplotlib.pyplot.figure(figsize=(self.width, self.height),
                                                             facecolor='white')
    def add_labels(self, label_pos=(-0.18, 1.01)):
        for axis, label in self.subplots.iteritems():
            axis.text(label_pos[0], label_pos[1], label.lower(), va='top',
                      transform=axis.transAxes, fontsize=12, fontweight='bold')
    def process_subplots(self, label_pos=(-0.18, 1.01), padding=0.02):
        self.inset_axes(padding=padding)
        self.add_labels(label_pos=label_pos)
        self.process_lines()


# University of Birmingham (UK) thesis guidelines state:
# - inner margin = 30mm, outer margin = 20mm, so double column = 160mm
# - top margin = bottom margin = 30mm, so maximum height = 237mm
class ThesisFigure(JournalFigure):
    def __init__(self, double_column=False, height=0.0):
        JournalFigure.__init__(self)
        print("Using figure dimensions for a thesis (University of Birmingham):")
        # width and height
        single, double, maximum = mm2inch(80), mm2inch(160), mm2inch(237)
        if double_column: self.width = double
        else:             self.width = single
        if height == 0.0: self.height = self.width
        elif height == 'single': self.height = single
        elif height == 'double': self.height = double
        elif height == 'max': self.height = maximum
        else:             self.height = min(height, maximum)
        print("width x height = "+str(self.width)+" x "+str(self.height))
        # default settings
        matplotlib.rc('font',**{'family':'Arial','weight':'normal','size':12})
        matplotlib.rc('mathtext', fontset='stixsans', default='regular')
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        # the figure!
        self.fig = matplotlib.pyplot.figure(figsize=(self.width, self.height),
                                                             facecolor='white')
    def add_labels(self, label_pos=(-0.18, 1.01)):
        for axis, label in self.subplots.iteritems():
            axis.text(label_pos[0], label_pos[1], label, fontsize=12, va='top',
                                   transform=axis.transAxes, fontweight='bold')
    def process_subplots(self, label_pos=(-0.18, 1.01), padding=0.02):
        self.inset_axes(padding=padding)
        self.add_labels(label_pos=label_pos)
        self.process_lines()


class SlideFigure(JournalFigure):
    def __init__(self, height=7.5, width=10):
        JournalFigure.__init__(self)
        self.height = height
        self.width = width
        matplotlib.rc('font',**{'family':'sans serif',
                                'weight':'normal',
                                'size':10})
        matplotlib.rc('mathtext', fontset='stixsans', default='regular')
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        self.fig = matplotlib.pyplot.figure(
                    figsize=(self.width, self.height))
    def process_subplots(self, padding=0.02):
        self.inset_axes(padding=padding)
        for axis, label in self.subplots.iteritems():
            axis.text(-0.25, 1.02, label, transform=axis.transAxes,
                                    va='top', fontsize=12, fontweight='bold')
            for line in axis.spines.itervalues():
                line.set_linewidth(1)
            axis.tick_params(bottom='on', top='off', left='on', right='off')
            for line in axis.get_xticklines() + axis.get_yticklines():
                line.set_markeredgewidth(1)


class PovrayFile:
    def __init__(self, agent_output, save_path):
        self.agent_output = agent_output
        self.save_path = os.path.join( \
                        toolbox_basic.check_path(os.path.dirname(save_path)),
                        os.path.basename(save_path))
        nI, nJ = agent_output.grid_nI, agent_output.grid_nJ
        nK, res = agent_output.grid_nK, agent_output.grid_res
        max_dimension = max(nI, nJ, nK)
        x_scale = nI / max_dimension
        y_scale = nJ / max_dimension
        self.scaling = res * max_dimension
        ### Set the header
        self.script = '#declare white = color rgb < 1.0, 1.0, 1.0 >;\n\n'
        self.script += 'background { white }\n'
        # The camera
        self.script += 'camera {\n'
        self.script += '\t location < 0.5 %f 1 >\n'%(y_scale/2)
        self.script += '\t look_at < 0.5, %f -1 >\n'%(y_scale/2)
        self.script += '\t up < 0.0, %f, 0.0 >\n'%(1.2*y_scale)
        self.script += '\t right < %f, 0.0, 0.0 >\n}\n\n'%(1.2*x_scale)
        # The lights
        for i in [0, 1]:
            self.script += 'light_source {\n\t < %d, 1, 1 >\n'%(i)
            self.script += '\t white\n\t shadowless}\n\n'
        # The substratum
        self.script += 'box {\n\t < 0.0, 0.0, 0.0 >\n'
        self.script += '\t< 1.0, -0.01, -0.001 >\n'
        self.script += '\t pigment { color rgb < 0.25, 0.25, 0.25 > }\n}\n\n'
    def append_cell(self, cell_output, total_radius=True):
        # You MUST have already set cell_output.vars['color']
        x = float(cell_output.vars['locationX']) / self.scaling
        y = float(cell_output.vars['locationY']) / self.scaling
        z = float(cell_output.vars['locationZ']) / self.scaling
        r_name = 'totalRadius' if total_radius else 'radius'
        r = float(cell_output.vars[r_name]) / self.scaling
        self.script += 'sphere {\n'
        self.script += '\t < %f, %f, %f >\n\t %f\n'%(x, y, z, r)
        self.script += '\t pigment { %s }\n}\n\n'%(cell_output.vars['colorRGB'])
    def save(self):
        print('Saving at %s'%(self.save_path))
        with open(self.save_path, 'wb') as f:
            f.write(self.script)
    def call(self, x_pixels=1600, y_pixels=1000):
        POVcmd = '''povray '%s' +A -D -geometry %dx%d''' \
                        %(self.save_path, x_pixels, y_pixels)
        os.system(POVcmd)
    def delete(self):
        os.remove(self.save_path)


class PlotPoint:
    def __init__(self, x, y,
                edgecolor='k', facecolor='none', size=5, style='o'):
        self.x          = float(x)
        self.y          = float(y)
        self.edgecolor  = edgecolor
        self.facecolor  = facecolor
        self.markersize = float(size)
        self.style      = style
    def plot(self, axis):
        axis.plot([self.x], [self.y], self.style, color=self.facecolor,
                    markeredgecolor=self.edgecolor, markersize=self.markersize)


class ScatterColumn:
    def __init__(self, x, name):
        self.x          = float(x)
        self.name       = name
        self.max_width  = 0.5
        self.edgecolor  = 'k'
        self.facecolor  = 'white'
        self.markersize = 5
        self.min_dist   = 0.01
        self.style      = 'o'
        self.points     = []
    def set_defaults(self, edgecolor=None, facecolor=None, max_width=None,
                                   markersize=None, min_dist=None, style=None):
        self.edgecolor  = self.edgecolor  if (edgecolor  == None) else edgecolor
        self.facecolor  = self.facecolor  if (facecolor  == None) else facecolor
        self.max_width  = self.max_width  if (max_width  == None) else max_width
        self.markersize = self.markersize if (markersize == None) else markersize
        self.min_dist   = self.min_dist   if (min_dist   == None) else min_dist
        self.style      = self.style      if (style      == None) else style
    def add_point(self, y, edgecolor=None, facecolor=None, style=None):
        edgecolor = self.edgecolor if (edgecolor == None) else edgecolor
        facecolor = self.facecolor if (facecolor == None) else facecolor
        style     = self.style     if (style     == None) else style
        self.points.append(PlotPoint(self.x, y, edgecolor=edgecolor,
                             facecolor=facecolor, style=style, size=self.markersize))
    def jiggle_points(self):
        max_x, min_x = self.x + self.max_width/2, self.x - self.max_width/2
        more_to_do = (len(self.points) > 0)
        while more_to_do:
            for i in range(len(self.points)):
                current = self.points.pop(0)
                points_array = numpy.array([[p.x, p.y] for p in self.points])
                kdtree = KDTree(points_array)
                ids = kdtree.query_ball_point([current.x, current.y], self.min_dist)
                # If there are any neighbours too near
                counter = 0
                while not ids == []:
                    rand = random.uniform(-self.min_dist, self.min_dist)
                    current.x = min(max_x, max(min_x, current.x + rand))
                    ids = kdtree.query_ball_point([current.x, current.y], self.min_dist)
                    counter += 1
                    if counter > 100: break
                self.points.append(current)
            for i in range(len(self.points)):
                current = self.points.pop(0)
                points_array = numpy.array([[p.x, p.y] for p in self.points])
                kdtree = KDTree(points_array)
                ids = kdtree.query_ball_point([current.x, current.y], self.min_dist)
                self.points.append(current)
                more_to_do = (not ids == [])
                if more_to_do: break
    def plot(self, axis):
        for point in self.points:
            point.plot(axis)


class ScatterRow:
    def __init__(self, y, name):
        self.y          = float(y)
        self.name       = name
        self.max_height  = 0.5
        self.edgecolor  = 'k'
        self.facecolor  = 'white'
        self.markersize = 5
        self.min_dist   = 0.01
        self.style      = 'o'
        self.points     = []
    def set_defaults(self, edgecolor=None, facecolor=None, max_height=None,
                                   markersize=None, min_dist=None, style=None):
        self.edgecolor  = self.edgecolor  if (edgecolor  == None) else edgecolor
        self.facecolor  = self.facecolor  if (facecolor  == None) else facecolor
        self.max_height = self.max_height if (max_height == None) else max_height
        self.markersize = self.markersize if (markersize == None) else markersize
        self.min_dist   = self.min_dist   if (min_dist   == None) else min_dist
        self.style      = self.style      if (style      == None) else style
    def add_point(self, x, edgecolor=None, facecolor=None, style=None):
        edgecolor = self.edgecolor if (edgecolor == None) else edgecolor
        facecolor = self.facecolor if (facecolor == None) else facecolor
        style     = self.style     if (style     == None) else style
        self.points.append(PlotPoint(x, self.y, edgecolor=edgecolor,
                             facecolor=facecolor, style=style, size=self.markersize))
    def jiggle_points(self):
        max_y, min_y = self.y + self.max_height/2, self.y - self.max_height/2
        more_to_do = (len(self.points) > 0)
        while more_to_do:
            for i in range(len(self.points)):
                current = self.points.pop(0)
                points_array = numpy.array([[p.x, p.y] for p in self.points])
                kdtree = KDTree(points_array)
                ids = kdtree.query_ball_point([current.x, current.y], self.min_dist)
                # If there are any neighbours too near
                counter = 0
                while not ids == []:
                    rand = random.uniform(-self.min_dist, self.min_dist)
                    current.y = min(max_y, max(min_y, current.y + rand))
                    ids = kdtree.query_ball_point([current.x, current.y], self.min_dist)
                    counter += 1
                    if counter > 100: break
                self.points.append(current)
            for i in range(len(self.points)):
                current = self.points.pop(0)
                points_array = numpy.array([[p.x, p.y] for p in self.points])
                kdtree = KDTree(points_array)
                ids = kdtree.query_ball_point([current.x, current.y], self.min_dist)
                self.points.append(current)
                more_to_do = (not ids == [])
                if more_to_do: break
    def plot(self, axis):
        for point in self.points:
            point.plot(axis)


def padding_axis(axis, side="right", pad=0.04):
    divider = make_axes_locatable(axis)
    cax = divider.append_axes(side, size="5%", pad=pad, axisbg='none')
    return cax


def empty_padding_axis(axis, side="right"):
    cax = padding_axis(axis, side=side)
    cax.set_xticklabels(['']*10)
    cax.set_yticklabels(['']*10)
    for spine in ['right', 'left', 'top', 'bottom']:
        cax.spines[spine].set_color('none')
    cax.tick_params(top="off", bottom="off", right="off", left="off")
    return cax


def make_colorbar(axis, colorscheme, side="top", fontsize=10, pad=0.04, label='',ticks=[0, 0.5, 1],ticklab=[0,0.5,1]):   
    divider = make_axes_locatable(axis)
    cax1 = divider.append_axes("top", size="5%", pad=pad)
    orientation = 'vertical' if side in ["right", "left"] else 'horizontal'
    cbar = pyplot.colorbar(colorscheme, cax=cax1, orientation=orientation, label=label, ticks=ticks)
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.ax.set_xticklabels(ticklab)
    return cbar

def make_colorbar2(axis, colorscheme, side="bottom", fontsize=10, pad=3, label='',ticks=[0, 0.5, 1],ticklab=[0,0.5,1]):   
    divider = make_axes_locatable(axis)
    cax1 = divider.append_axes("top", size="5%", pad=3)
    orientation = 'vertical' if side in ["right", "left"] else 'horizontal'
    cbar = pyplot.colorbar(colorscheme, cax=cax1, orientation=orientation, label=label, ticks=ticks)
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.ax.set_xticklabels(ticklab)
    return cbar


def colorbar_on_right(axis, colorscheme, fontsize=10):
    #divider = make_axes_locatable(axis)
    #cax = divider.append_axes("right", size="5%", pad=0.04)
    #cax = space_on_right(axis)
    #cbar = pyplot.colorbar(colorscheme, cax=cax)
    #cbar.ax.tick_params(labelsize=fontsize)
    print('\n\ntoolbox_plotting.colorbar_on_right() is deprecated')
    print('Please use toolbox_plotting.make_colorbar() instead\n\n')
    cbar = make_colorbar(axis, colorscheme, side="right")
    return cbar


def draw_linear_regression(axis, color, x_vals, y_vals):
    slope, intercept, r_value, p_value, std_err = linregress(x_vals, y_vals)
    x1, x2 = min(x_vals), max(x_vals)
    y1, y2 = x1*slope + intercept, x2*slope + intercept
    axis.plot([x1, x2], [y1, y2], color)
    return r_value, p_value, std_err


def distinguishable_colors(number, cmap='hsv'):
    rgba = matplotlib.pyplot.get_cmap(cmap)(numpy.linspace(0.0, 1.0, number))
    html = []
    for old in rgba:
        new = "#"
        for i in range(3):
            temp = hex(int(old[i]*255))[2:]
            if len(temp) == 1:
                temp = "0"+temp
            new += temp
        html.append(new)
    return html

def plot_color_dictionary(color_dict, file_path):
    y_scale = 0.3
    x_scale = 0.15
    height = y_scale*(len(color_dict)+1)
    width = x_scale * (max([len(n) for n in color_dict.keys()]) + 1)
    fig = SlideFigure(height=height, width=width)
    sub = fig.add_subplot("", 111, frameon=False)
    sub.set_xlim([0, width])
    sub.set_ylim([0, height])
    y = height - y_scale
    for name in sorted(color_dict.keys()):
        color = color_dict[name]
        sub.plot([0.0, 1.0], [y]*2, '-', color=color)
        sub.text(1.0 + x_scale, y, name, color='black', va='center', ha='left')
        y -= y_scale
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    fig.save(file_path)
