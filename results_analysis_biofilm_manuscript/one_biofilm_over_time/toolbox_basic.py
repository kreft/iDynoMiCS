#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import glob
import numpy
import os
import platform
import shutil
import xml.etree.ElementTree as xmlTree
import zipfile


pi = numpy.pi


class LogFile:
    def __init__(self, save_path):
        self.text = ''
        self.path = save_path
    
    def add(self, message):
        print(message)
        self.text += message+'\n'
    
    def save(self):
        with open(self.path, 'w') as f:
            f.write(self.text)


def are_dicts_same(dict1, dict2):
    return is_enclosed_dict(dict1, dict2) and is_enclosed_dict(dict2, dict1)


### TODO
def are_xml_trees_same(tree1, tree2, diffs_allowed=[]):
    treeA = tree1
    rs = treeA.find("./process/param[@name='randomSeed']")
    treeA.getroot().remove(rs)
    exit()
    #xmlTree.dump(tree1)
    def check_parent(parent):
        if parent == None: return
        print(parent.tag, parent.attrib)
        print('#'*10)
        for child in parent:
            check_parent(child)
    check_parent(tree1.getroot())
    return False


def copy_dir(source, destination, overwrite=False, logfile=None):
    if not os.path.isdir(source):
        error_message('Source dir is not valid:', source, logfile=logfile)
    make_dir(destination, logfile=logfile)
    source_file_list = file_list(source)
    for file in source_file_list:
        shutil.copy(file, destination)


def check_path(path):
    path = os.path.expanduser(path)
    path = os.path.abspath(path)
    if os.path.exists(path):
        return os.path.normpath(path)
    elif os.path.lexists(path):
        error_message('Path may be a broken link:', path)
    else:
        error_message('Path may not exist:', path)
    return None


def error_message(message, path, logfile=None):
    out_msg = '#\n# Error!\n# %s\n# %s\n#\n' %(message, path)
    if logfile == None:
        print(out_msg)
    else:
        logfile.add(out_msg)
        logfile.save()


def file_list(path, filetype='*'):
    path = check_path(path)
    if not ('*' in filetype):
        filetype += '*'
    filelist = [filename for filename in glob.glob(path+os.sep+filetype) if
                                 os.path.isfile(os.path.join(path, filename))]
    #if filelist == []:
    #    error_message('Could not find any files in:', path)
    return sorted(filelist)


def find_index_of_row_in_array(array, row):
    counter = 0
    for i_row in array:
        if numpy.array_equal(i_row, row):
            return counter
        counter += 1
    return -1


def find_protocol_file_path(path):
    path = check_path(path)
    if os.path.splitext(path)[1] == '.xml':
        return path
    sim_name = os.path.basename(path)
    potential_files = []
    i = 0
    while i < len(sim_name) and not len(potential_files) == 1:
        potential_files = file_list(path, sim_name[:i]+'*.xml')
        i += 1
    if potential_files == []:
        error_message('Could not find any potential protocol files in:', path)
        return None
    elif len(potential_files) > 1:
        error_message('Found too many potential protocol files in:', path)
        return None
    else:
        return potential_files[0]


def get_key_from_dict_value(dictionary, value):
    for key, val in dictionary.iteritems():
        if val == value:
            return key
    error_message('Could not find '+str(value), 'in '+str(dictionary))
    return None


def get_xml_tree(path):
    if path == None:
        return None
    if isinstance(path, xmlTree.ElementTree):
        return path
    if isinstance(path, xmlTree.Element):
        return xmlTree.ElementTree(path)
    if os.path.isfile(path):
        try:
            tree = xmlTree.parse(path)
        except xmlTree.ParseError:
            error_message('Problem parsing', path)
        return tree
    return None


def is_enclosed_dict(super_dict, sub_dict):
    for key in sub_dict.keys():
        if super_dict.keys().count(key) == 0:
            return False
        if not sub_dict[key] == super_dict[key]:
            return False
    return True


def is_in_range(value, range_values, incl_endpoints=True):
    if not (len(range_values) == 2 and range_values[0] < range_values[1]):
        print('toolbox_basic.is_in_range():')
        print('Range must be defined by two subsequent values!')
        print(range_values)
        return
    if incl_endpoints:
        return (value >= range_values[0]) and (value <= range_values[1])
    else:
        return (value > range_values[0]) and (value < range_values[1])


def load_array(path):
    #path = check_path(path)
    if path[-4:] == '.npy':
        return numpy.load(check_path(path))
    if path[-4:] == '.txt':
        print('toolbox_basic.load_array() TO DO: loading .txt files! '+path)
    if os.path.isfile(path+'.npy'):
        return numpy.load(check_path(path+'.npy'))
    print('toolbox_basic.load_array() TO DO: loading unrecognised files! '+path)


def make_dir(path, logfile=None):
    if not os.path.isdir(path):
        if platform.system == 'Windows':
            try:
                os.mkdir(path)
            except WindowsError:
                error_message('Could not create:', path, logfile=logfile)
        else:
            try:
                os.mkdir(path)
            except OSError:
                error_message('Could not create:', path, logfile=logfile)
    return path


def make_symbolic_link(source, link_name, logfile=None):
    source = check_path(source)
    if platform.system == 'Windows':
        print('toolbox_basic.make_symbolic_link() does not yet support Windows')
    else:
        try:
            os.symlink(source, link_name)
        except OSError:
            error_message('Could not make symbolic link from %s to'%(source),
                                                    link_name, logfile=logfile)


def move_file(source, destination, overwrite=False, logfile=None):
    if not os.path.isfile(source):
        error_message('Source file is not valid:', source, logfile=logfile)
    if not os.path.isdir(os.path.dirname(destination)):
        error_message('Destination may not exist:', destination, logfile=logfile)
    if os.path.isdir(destination):
        destination = os.path.join(destination, os.path.basename(source))
    if os.path.isfile(destination):
        if overwrite:
            os.remove(destination)
            shutil.move(source, destination)
        else:
            msg = 'Cannot move %s to %s without permission to overwrite' \
                    %(source, destination)
            if logfile == None: print(msg)
            else: logfile.add(msg)
    else:
        shutil.move(source, destination)


def are_nearly_equal(a, b, sig_fig=3):
    return (a == b) or (int(a*10**sig_fig) == int(b*10**sig_fig))


def rm_dir(path):
    path = check_path(path)
    if platform.system == 'Windows':
        try:
            shutil.rmtree(path)
            print('Removing '+path)
        except WindowsError:
            error_message('Could not remove:', path)
    else:
        try:
            print('Removing '+path)
            shutil.rmtree(path)
        except OSError:
            error_message('Could not remove:', path)


def pad_filenames(path, filetype='*'):
    path = check_path(path)
    filelist = file_list(path, filetype)
    longest_num = 1
    for filename in filelist:
        open_bracket, close_bracket = filename.rfind('('), filename.rfind(')')
        this_length = close_bracket - open_bracket - 1
        longest_num = max(longest_num, this_length)
    for filename in filelist:
        open_bracket, close_bracket = filename.rfind('('), filename.rfind(')')
        this_length = close_bracket - open_bracket - 1
        padding = '0'*(longest_num - this_length)
        shutil.move(filename,
                  filename[:open_bracket+1]+padding+filename[open_bracket+1:])


def prettify_xml_file(path):
    path = check_path(path)
    with open(path, 'Ur') as f:
        text = f.read()
    text = text.replace('><', '>\n<')
    with open(path, 'w') as f:
        f.write(text)


def save_array(array, save_path):
    numpy.save(save_path+'.npy', array)
    str_array = '\n'.join(' '.join(str(value) for value in row) for row in array)
    with open(save_path+'.txt', 'w') as f:
        f.write(str_array)


def seconds_to_hms(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return (h, m, s)


def subdir_list(path, dirtype='*'):
    path = check_path(path)
    if not ('*' in dirtype):
        dirtype += '*'
    subdirs = [directory for directory in glob.glob(path+os.sep+dirtype) if
                                 os.path.isdir(os.path.join(path, directory))]
    #if subdirs == []:
    #    error_message('Could not find any subdirectories of type '+dirtype+
    #                                                ' in:', path)
    return sorted(subdirs)


def typecast_string(string):
    if not isinstance(string, str):
        return string
    try:
        out = int(string)
    except ValueError:
        try:
            out = float(string)
        except ValueError:
            if string.lower() == 'true': out = True
            elif string.lower() == 'false': out = False
            else: out = string
    return out


def unzip_files(path, overwrite=False):
    path = check_path(path)
    out_dir = path[:-4]
    if zipfile.is_zipfile(path):
        zfile = zipfile.ZipFile(path,'r')
        make_dir(out_dir)
        if (not overwrite) and os.listdir(out_dir):
            print('Skipping extraction of '+path)
        else:
            zfile.extractall(out_dir)
            print('Files from '+path)
            print(' extracted to '+out_dir)
        return out_dir
    else:
        error_message('Path is not a zipfile:', path)


def zip_file(file_path, zip_path=None, rm_file=False):
    file_path = check_path(file_path)
    if os.path.isfile(file_path):
        if zip_path == None:
            zip_path = os.path.splitext(file_path)[0]+'.zip'
    else:
        print('\n# file_path is not a file!\n'+file_path+'\n')
        return
    if zipfile.is_zipfile(zip_path): arg = 'a'
    else: arg = 'w'
    with zipfile.ZipFile(zip_path, arg) as z:
        z.write(file_path, os.path.basename(file_path))
    if rm_file:
        os.remove(file_path)
