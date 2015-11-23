#!/usr/bin/env python
import sys
import argparse


def splitHalf(fn):
	with open(fn) as fi:
		# count number of ion blocks in the file
		ion_blocks = sum(1 for line in fi if line.strip().startswith('END IONS'))
		# 'rewind' file
		fi.seek(0)
		ion_blocks_written = 0
		with open(fn + '.1', 'wb') as fo1, open(fn + '.2', 'wb') as fo2:
			fo = fo1
			while 1:
				if ion_blocks_written == ion_blocks / 2:
					# if half of the blocks have been written to file, switch the output file
					fo = fo2
				try:
					line = fi.next()
				except:
					# if end of file is reached, stop
					break
				fo.write(line)
				if line.strip().startswith('END IONS'):
					ion_blocks_written += 1
	pass

def filterIntensities(fn, min_intensity=10.0):
	with open(fn) as fi, open(fn + '.filtered', 'wb') as fo:
		for line in fi:
			if line[0].isdigit() and len(line.strip().split()) == 2:
				if float(line.strip().split()[1]) < min_intensity:
					continue
			fo.write(line)

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument('--split-half', action='store_true')
	ap.add_argument('--min-intensity', type=float, default=10.0)
	ap.add_argument('mgffile')
	args = ap.parse_args()


	if args.split_half:
		splitHalf(args.mgffile)
	else:
		filterIntensities(args.mgffile, min_intensity=args.min_intensity)


	pass

if __name__ == '__main__': main()
