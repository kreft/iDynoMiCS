#!/usr/bin/python
# Author: Chinmay Kanchi
import re, os, sys, shutil, glob

def pad_filenames(args):
	if args == []:
		args.append('.')

	file_list = glob.glob(args[0]+os.sep+'*.pov')

	numbers = []
	for filename in file_list:
		open_bracket = filename.rfind('(')
		close_bracket = filename.rfind(')')

		numbers.append(filename[open_bracket+1:close_bracket])

	max_len = max((len(number) for number in numbers))
	
	for filename in file_list:
		open_bracket = filename.rfind('(')
		close_bracket = filename.rfind(')')
		filename_len = len(filename[open_bracket+1:close_bracket])	
		padding = '0'*(max_len-filename_len)
		shutil.move(filename, filename[:open_bracket+1]+padding+filename[open_bracket+1:])

if __name__ == '__main__':
	pad_filenames(sys.argv[1:])
