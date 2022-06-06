#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import math
import numpy
import matplotlib

pi = math.pi


class Arrow:
    def __init__(self):
        self.y2xscale  = 1
        self.linewidth  = 1
        self.zorder     = 0
        self.head_angle = pi/4
        self.head_open  = True
        self.head_size  = 0.01
        self.color      = 'k'
        self.linestyle  = '-'
    def set_points(self, fletch_pos, head_pos):
        self.x1, self.y1 = fletch_pos[0], fletch_pos[1]
        self.x2, self.y2 = head_pos[0], head_pos[1]
        self.angle = math.atan2(self.y2 - self.y1, self.x2 - self.x1)
    def set_defaults(self, color=None, head_angle=None, head_open=None,
            head_size=None, linestyle=None, linewidth=None, y2xscale=None, zorder=None):
        if not (color      == None): self.color      = color
        if not (head_angle == None): self.head_angle = head_angle
        if not (head_open  == None): self.head_open  = head_open
        if not (head_size  == None): self.head_size  = head_size
        if not (linewidth  == None): self.linewidth  = linewidth
        if not (linestyle  == None): self.linestyle  = linestyle
        if not (y2xscale   == None): self.y2xscale   = y2xscale
        if not (zorder     == None): self.zorder     = zorder
    def draw(self, axis, head_size=None, head_angle=None, zorder=None,
                head_open=None, linewidth=None, color=None, linestyle=None):
        linewidth  = self.linewidth  if (linewidth  == None) else linewidth
        zorder     = self.zorder     if (zorder     == None) else zorder
        head_angle = self.head_angle if (head_angle == None) else head_angle
        head_open  = self.head_open  if (head_open  == None) else head_open
        head_size  = self.head_size  if (head_size  == None) else head_size
        color      = self.color      if (color      == None) else color
        linestyle = self.linestyle   if (linestyle  == None) else linestyle
        axis.plot([self.x1, self.x2], [self.y1, self.y2], zorder=zorder,
                      solid_capstyle='round', linestyle=linestyle, lw=linewidth, color=color)
        x_vals = [self.x2 + head_size*math.cos(self.angle + pi - head_angle),
                  self.x2,
                  self.x2 + head_size*math.cos(self.angle - pi + head_angle)]
        pos_diff = self.y2xscale * head_size * math.sin(self.angle + pi - head_angle)
        neg_diff = self.y2xscale * head_size * math.sin(self.angle - pi + head_angle)
        y_vals = [self.y2 + pos_diff, self.y2, self.y2 + neg_diff]
        if head_open:
            axis.plot(x_vals, y_vals, zorder=zorder, solid_capstyle='round',
                                linestyle=linestyle, lw=linewidth, color=color)
        else:
            head = ArrowHead()
            head.set_defaults(edgecolor=self.color, facecolor=self.color,
                        head_angle=self.head_angle, y2xscale=self.y2xscale,
                                    zorder=self.zorder)
            head.set_points((self.x2, self.y2), self.angle, head_size)
            head.draw(axis, facecolor=color, edgecolor=color, zorder=zorder)


class Shape:
    def __init__(self):
        self.y2xscale  = 1
        self.alpha      = 1.0
        self.edgecolor  = 'k'
        self.facecolor  = 'none'
        self.linestyle  = '-'
        self.linewidth  = 1
        self.zorder     = 0
        self.transform  = False
    def set_defaults(self, alpha=None, edgecolor=None, facecolor=None, linestyle=None,
                linewidth=None, transform=None, y2xscale=None, zorder=None):
        if not (alpha      == None): self.alpha      = alpha
        if not (edgecolor  == None): self.edgecolor  = edgecolor
        if not (facecolor  == None): self.facecolor  = facecolor
        if not (linewidth  == None): self.linewidth  = linewidth
        if not (linestyle  == None): self.linestyle  = linestyle
        if not (transform  == None): self.transform  = transform
        if not (y2xscale   == None): self.y2xscale   = y2xscale
        if not (zorder     == None): self.zorder     = zorder
    def set_points(self):
        self.vals = []
    def draw(self, axis, alpha = None, edgecolor=None, facecolor=None, linestyle=None,
                                    transform=None, linewidth=None, zorder=None):
        alpha     = self.alpha     if (alpha     == None) else alpha
        #linestyle = self.linestyle if (linestyle == None) else linestyle
        linewidth = self.linewidth if (linewidth == None) else linewidth
        transform = self.transform if (transform == None) else transform
        facecolor = self.facecolor if (facecolor == None) else facecolor
        edgecolor = self.edgecolor if (edgecolor == None) else edgecolor
        zorder    = self.zorder    if (zorder    == None) else zorder
        if transform:
            axis.add_patch(matplotlib.patches.Polygon(self.vals, alpha=alpha,
                            closed=True, facecolor=facecolor,
                            edgecolor=edgecolor, zorder=zorder, 
                            linewidth=linewidth, transform=axis.transAxes))
        else:
            axis.add_patch(matplotlib.patches.Polygon(self.vals, closed=True,
                            facecolor=facecolor, edgecolor=edgecolor,
                            zorder=zorder, alpha=alpha,
                            linewidth=linewidth))


class CircleSegment(Shape):
    def __init__(self):
        Shape.__init__(self)
        self.npoints = 100
    
    def set_defaults(self, alpha=None, edgecolor=None, facecolor=None, linestyle=None,
                     linewidth=None, npoints=None, transform=None, y2xscale=None, zorder=None):
        Shape.set_defaults(self, alpha=alpha,
                                 edgecolor=edgecolor,
                                 facecolor=facecolor,
                                 linestyle=linestyle,
                                 linewidth=linewidth,
                                 transform=transform,
                                 y2xscale=y2xscale,
                                 zorder=zorder)
        if not (npoints == None):
            self.npoints = npoints
    
    def set_points(self, focus_pos, radius, angles):
        Shape.set_points(self)
        if not len(angles) == 2:
            print('CircleSegment must be given two angles! '+str(angles))
            return
        for angle in numpy.linspace(angles[0], angles[1], self.npoints):
            x = focus_pos[0] + radius*math.cos(angle)
            y = focus_pos[1] + self.y2xscale * radius * math.sin(angle)
            self.vals.append((x, y))


class Circle(CircleSegment):
    def __init__(self):
        CircleSegment.__init__(self)
    def set_points(self, focus_pos, radius):
        CircleSegment.set_points(self, focus_pos, radius, [-pi, pi])


class Rectangle(Shape):
    def __init__(self):
        Shape.__init__(self)
    def set_points(self, sw_corner, ne_corner):
        Shape.set_points(self)
        self.vals.append((sw_corner[0], sw_corner[1]),
                         (ne_corner[0], sw_corner[1]),
                         (ne_corner[0], ne_corner[1]),
                         (sw_corner[0], ne_corner[1]))

'''
class Cell(Shape):
    def __init__(self):
        Shape.__init__(self)
        self.npoints = 100
    def set_defaults(self, edgecolor=None, facecolor=None, linestyle=None,
                     linewidth=None, npoints=None, transform=None, y2xscale=None, zorder=None):
        Shape.set_defaults(self, edgecolor=edgecolor,
                                 facecolor=facecolor,
                                 linestyle=linestyle,
                                 linewidth=linewidth,
                                 transform=transform,
                                 y2xscale=y2xscale,
                                 zorder=zorder)
        if not (npoints == None): self.npoints = npoints
    def set_points(self, focusA, focusB, radius):
        Shape.set_points(self)
        y2xscale = self.y2xscale
        # focusA should be to the left of focusB
        if focusA[0] > focusB[0]:
            temp_focus = focusA
            focusA = focusB
            focusB = temp_focus
        self.focusA, self.focusB, self.radius = focusA, focusB, radius
        # Vertical cell
        if focusA[0] == focusB[0]:
            max_y, min_y = max(focusA[1], focusB[1]), min(focusA[1], focusB[1])
            for angle in numpy.linspace(pi, 0, self.npoints):
                self.x_vals.append(focusA[0] + radius*math.cos(angle))
                y_diff = self.y2xscale * radius * math.sin(angle)
                self.y1_vals.append(max_y + y_diff)
                self.y2_vals.append(min_y - y_diff)
            return
        # Horizontal cell
        if focusA[1] == focusB[1]:
            y = focusA[1]
            for angle in numpy.linspace(pi, pi/2, int(self.npoints/2)):
                self.x_vals.append(focusA[0] + radius*math.cos(angle))
                y_diff = self.y2xscale * radius * math.sin(angle)
                self.y1_vals.append(y + y_diff)
                self.y2_vals.append(y - y_diff)
            for angle in numpy.linspace(pi/2, 0, int(self.npoints/2)):
                self.x_vals.append(focusB[0] + radius*math.cos(angle))
                y_diff = self.y2xscale * radius * math.sin(angle)
                self.y1_vals.append(y + y_diff)
                self.y2_vals.append(y - y_diff)
            return
        # Find the angle that the line from focusA to focusB makes with the
        # horizontal.
        self.spine_angle = math.atan2(focusB[1]-focusA[1], focusB[0]-focusA[0])
        tan_spine = math.tan(self.spine_angle)
        ### Do the part around focusA that is all circular.
        start_angle = - pi
        end_angle = - abs(self.spine_angle) - pi/2
        temp_n = int(self.npoints*abs(start_angle - end_angle)/pi)
        #print temp_n
        for angle in numpy.linspace(start_angle, end_angle, temp_n):
            self.x_vals.append(focusA[0] + radius*math.cos(angle))
            y_diff = self.y2xscale * radius * math.sin(angle)
            self.y1_vals.append(focusA[1] - y_diff)
            self.y2_vals.append(focusA[1] + y_diff)
        ### Do the part around focusA where one side is straight.
        temp_x = focusA[0] + radius*math.cos(end_angle)
        temp_y = math.copysign(radius*math.sin(end_angle), self.spine_angle)
        temp_y = focusA[1] + temp_y
        start_angle = - self.spine_angle - pi/2
        end_angle = self.spine_angle - pi/2
        temp_n = int(self.npoints*abs(start_angle - end_angle)/pi)
        #print temp_n
        for angle in numpy.linspace(start_angle, end_angle, temp_n):
            x = focusA[0] + radius*math.cos(angle)
            self.x_vals.append(x)
            y_diff = self.y2xscale * radius * math.sin(angle)
            if self.spine_angle > 0:
                self.y1_vals.append(temp_y + (x-temp_x)*tan_spine)
                self.y2_vals.append(focusA[1] + y_diff)
            else:
                self.y1_vals.append(focusA[1] - y_diff)
                self.y2_vals.append(temp_y + (x-temp_x)*tan_spine)
        ### Do the part around focusB where one side is straight.
        start_angle = self.spine_angle + pi/2
        end_angle = pi/2 - abs(self.spine_angle)
        temp_x = focusB[0] + radius*math.cos(end_angle)
        temp_y = math.copysign(radius*math.sin(start_angle), self.spine_angle)
        temp_y = focusB[1] - temp_y
        end_angle = pi/2 - self.spine_angle
        for angle in numpy.linspace(start_angle, end_angle, temp_n):
            x = focusB[0] + radius*math.cos(angle)
            self.x_vals.append(x)
            y_diff = self.y2xscale * radius * math.sin(angle)
            if self.spine_angle > 0:
                self.y2_vals.append(temp_y + (x-temp_x)*tan_spine)
                self.y1_vals.append(focusB[1] + y_diff)
            else:
                self.y2_vals.append(focusB[1] - y_diff)
                self.y1_vals.append(temp_y + (x-temp_x)*tan_spine)
        ### Do the part around focusA that is all circular.
        start_angle = pi/2 - abs(self.spine_angle)
        end_angle = 0
        temp_n = int(self.npoints*abs(start_angle - end_angle)/pi)
        #print temp_n
        for angle in numpy.linspace(start_angle, end_angle, temp_n):
            self.x_vals.append(focusB[0] + radius*math.cos(angle))
            y_diff = self.y2xscale * radius * math.sin(angle)
            self.y2_vals.append(focusB[1] - y_diff)
            self.y1_vals.append(focusB[1] + y_diff)


class ArrowHead(Shape):
    def __init__(self):
        Shape.__init__(self)
        self.head_angle = pi/4
    def set_defaults(self, edgecolor=None, facecolor=None, head_angle=None,
                        linestyle=None, linewidth=None,
                        y2xscale=None, zorder=None):
        Shape.set_defaults(self, edgecolor=edgecolor, facecolor=facecolor,
                       linestyle=linestyle, linewidth=linewidth,
                       y2xscale=y2xscale, zorder=zorder)
        if not (head_angle == None): self.head_angle = head_angle
    def set_points(self, corner, spine_angle, side_length):
        Shape.set_points(self)
        self.spine_angle = spine_angle % (2*pi)
        x, y = corner[0], corner[1]
        self.x_vals.append(x)
        pos_x_diff = side_length*math.cos(self.spine_angle + pi - self.head_angle)
        neg_x_diff = side_length*math.cos(self.spine_angle - pi + self.head_angle)
        pos_y_diff = self.y2xscale*side_length*math.sin(self.spine_angle + pi - self.head_angle)
        neg_y_diff = self.y2xscale*side_length*math.sin(self.spine_angle - pi + self.head_angle)
        # Pointing left/right
        if (self.spine_angle == 0) or (self.spine_angle == pi):
            self.x_vals = [x, x + pos_x_diff]
            self.y1_vals = [y, y + pos_y_diff]
            self.y2_vals = [y, y + neg_y_diff]
        # Pointing up/down
        elif (self.spine_angle == 0.5*pi) or (self.spine_angle == 1.5*pi):
            self.x_vals = [x + neg_x_diff, x, x + pos_x_diff]
            self.y1_vals = [y + pos_y_diff, y, y + pos_y_diff]
            self.y2_vals = [y + pos_y_diff, y + pos_y_diff, y + pos_y_diff]
        else:
            self.x_vals = [x + neg_x_diff, x, x + pos_x_diff]
            self.y1_vals = [y + pos_y_diff, y, y + pos_y_diff]
            self.y2_vals = [y + pos_y_diff, y + pos_y_diff, y + pos_y_diff]


class Lemniscate(Shape):
    def __init__(self):
        Shape.__init__(self)
        self.npoints     = 100
    def set_defaults(self, edgecolor=None, facecolor=None, npoints=None,
                                        linestyle=None, linewidth=None,
                                            y2xscale=None, zorder=None):
        Shape.set_defaults(self, edgecolor=edgecolor, facecolor=facecolor,
                       linestyle=linestyle, linewidth=linewidth,
                       y2xscale=y2xscale, zorder=zorder)
        if not (npoints == None): self.npoints = npoints
    def set_points(self, center, width, height):
        Shape.set_points(self)
        for angle in numpy.linspace(0.0, math.pi, self.npoints):
            sin, cos = math.sin(angle), math.cos(angle)
            self.x_vals.append(center[0] + (0.5*width*cos/(1+sin**2)))
            self.y1_vals.append(center[1] + (0.5*height*sin*cos/(1+sin**2)))
            self.y2_vals.append(center[1] - (0.5*height*sin*cos/(1+sin**2)))
'''



class Shape3D:
    def __init__(self):
        self.edgecolor  = 'k'
        self.facecolor  = 'none'
        self.linestyle  = '-'
        self.linewidth  = 1
        self.zorder     = 0
    def set_defaults(self, edgecolor=None, facecolor=None, linestyle=None,
                linewidth=None, transform=None, y2xscale=None, zorder=None):
        if not (edgecolor  == None): self.edgecolor  = edgecolor
        if not (facecolor  == None): self.facecolor  = facecolor
        if not (linewidth  == None): self.linewidth  = linewidth
        if not (linestyle  == None): self.linestyle  = linestyle
        if not (zorder     == None): self.zorder     = zorder
    def set_points(self):
        self.x_vals, self.y_vals, self.z_vals = [], [], []
    def draw(self, axis, edgecolor=None, facecolor=None, linestyle=None,
                                    transform=None, linewidth=None, zorder=None):
        #linestyle = self.linestyle if (linestyle == None) else linestyle
        linewidth = self.linewidth if (linewidth == None) else linewidth
        facecolor = self.facecolor if (facecolor == None) else facecolor
        edgecolor = self.edgecolor if (edgecolor == None) else edgecolor
        zorder    = self.zorder    if (zorder    == None) else zorder
        axis.plot_surface(self.x_vals, self.y_vals, self.z_vals,
                          rstride=1, cstride=1, color=facecolor,
                          edgecolor=edgecolor, zorder=zorder,
                          linewidth=linewidth)


class Sphere(Shape3D):
    def __init__(self):
        Shape3D.__init__(self)
    def set_points(self, center_pos, radius, npoints=20j):
        u, v = numpy.mgrid[0:2*numpy.pi:2*npoints, 0:numpy.pi:npoints]
        self.x_vals = center_pos[0] + radius*numpy.cos(u)*numpy.sin(v)
        self.y_vals = center_pos[1] + radius*numpy.sin(u)*numpy.sin(v)
        self.z_vals = center_pos[2] + radius*numpy.cos(v)
        

class Box(Shape3D):
    def __init__(self):
        Shape3D.__init__(self)
    def set_points(self, x_range, y_range, z_range):
        [x0, x1] = x_range
        [y0, y1] = y_range
        [z0, z1] = z_range
        self.x_vals, self.y_vals, self.z_vals = \
            numpy.mgrid[x0:x1:(x1-x0), y0:y1:(y1-y0), z0:z1:(z1-z0)]