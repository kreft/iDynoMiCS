#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import numpy
import os
import toolbox_basic
import xml.etree.ElementTree as xmlTree


# ResultXMLFile is a general-purpose class describing XML files with data stored
# within as CSV.
class ResultXMLfile:
    def __init__(self, path=None, root_name='result', read_only=False):
        if path == None:
            self.path = path
            self.root = xmlTree.Element(root_name)
            self.tree = xmlTree.ElementTree(self.root)
        else:
            self.path = path
            self.tree = toolbox_basic.get_xml_tree(self.path)
            self.root = self.tree.getroot()
        self.path = path
        self.tree = toolbox_basic.get_xml_tree(self.path)
        if self.tree == None:
            self.root = xmlTree.Element(root_name)
            self.tree = xmlTree.ElementTree(self.root)
        else:
            self.root = self.tree.getroot()
        self.read_only = read_only
        self.header = ''
        try:
            if 'header' in self.root.attrib.keys():
                self.header = self.root.attrib['header']
        except AttributeError: pass
        self.lines = []
        self.subresults = []
    def write(self, output_path=None, prettify=True, update_text=True):
        if self.read_only:
            print('Results file is read-only!')
            return
        if update_text:
            self.update_text()
        self.root.attrib['header'] = self.header
        if output_path == None: output_path = self.path
        else:                   self.path = output_path
        with open(output_path, 'w') as f:
            self.tree.write(f, encoding='utf-8', xml_declaration=True)
        if prettify:
            toolbox_basic.prettify_xml_file(self.path)
    def get_subresult(self, searchpattern):
        sub_root = self.find(searchpattern)
        subresult = ResultXMLfile(path=sub_root)
        self.subresults.append(subresult)
        return subresult
    def find(self, search_pattern):
        return self.tree.find(search_pattern)
    def findall(self, search_pattern):
        return self.tree.findall(search_pattern)
    def display(self):
        xmlTree.dump(self.tree)
    def read_in_text(self, header=None):
        if self.root.text == None: return []
        if not header == None: self.header = header
        self.lines = [SingleCSVline(self.header, line) for line in \
                                                self.root.text.split(';')[:-1]]
        return self.lines
    def get_column_number(self, column_name):
        return self.header.split(',').index(column_name)
    # Checks whether the given column name is already in the header, and adds
    # if it isn't. Returns the column number.
    def add_column_name(self, column_name):
        try:
            return self.get_column_number(column_name)
        except ValueError:
            if self.header == '':
                self.header = column_name
                return 0
            self.header += ','+column_name
            self.root.attrib['header'] = self.header
            return self.get_column_number(column_name)
    # Return all values in a column with the given name.
    def get_column(self, column_name):
        return [line.vars[column_name] for line in self.lines]
    def text(self, column_names=None):
        if column_names == None: column_names = self.header.split(',')
        text = '\n' + ''.join( [ line.text(column_names=column_names) \
                                                    for line in self.lines ] )
        return text
    def update_text(self, recursive=True):
        self.root.text = self.text()
        if recursive:
            for subresult in self.subresults:
                subresult.update_text()
    def append_subresult(self, subresult):
        self.subresults.append(subresult)
        self.root.append(subresult.root)
    def remove_subresult(self, subresult):
        self.subresults.remove(subresult)
        try:
            self.root.remove(subresult.root)
        except ValueError: pass


class SingleCSVline:
    def __init__(self, header, text):
        self.vars = {}
        col_names = header.split(',')
        values = text.replace('\n','').replace(';','').split(',')
        num_cols = len(col_names)
        if not num_cols == len(values):
            print('Header and values have different lengths!')
            print(header+'\n'+text)
            exit()
        for i in range(num_cols):
            self.vars[col_names[i]] = toolbox_basic.typecast_string(values[i])
    def add_column(self, column_name, value):
        self.vars[column_name] = value
    def text(self, column_names=None):
        if column_names == None: column_names = self.vars.keys()
        text = ','.join([str(self.vars[col_name]) for col_name in column_names])
        text += ';\n'
        return text
