#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
# CCP4 map Output Module
# Gael Goret for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr
###################################################################################
# -----------------------------------------------------------------------------
# Write an MRC 2000/ccp4 format file.
#
# Header contains four byte integer or float values:
#
# 1	NX	number of columns (fastest changing in map)	
# 2	NY	number of rows					
# 3	NZ	number of sections (slowest changing in map)	
# 4	MODE	data type :					
# 		0	image : signed 8-bit bytes range -128 to 127
# 		1	image : 16-bit halfwords		
# 		2	image : 32-bit reals			
# 		3	transform : complex 16-bit integers	
# 		4	transform : complex 32-bit reals	
# 5	NXSTART	number of first column in map			
# 6	NYSTART	number of first row in map			
# 7	NZSTART	number of first section in map			
# 8	MX	number of intervals along X			
# 9	MY	number of intervals along Y			
# 10	MZ	number of intervals along Z			
# 11-13	CELLA	cell dimensions in angstroms			
# 14-16	CELLB	cell angles in degrees				
# 17	MAP# axis corresp to cols (1,2,3 for X,Y,Z)		
# 18	MAPR	axis corresp to rows (1,2,3 for X,Y,Z)		
# 19	MAPS	axis corresp to sections (1,2,3 for X,Y,Z)	
# 20	DMIN	minimum density value				
# 21	DMAX	maximum density value				
# 22	DMEAN	mean density value				
# 23	ISPG	space group number 0 or 1 (default=0)		
# 24	NSYMBT	number of bytes used for symmetry data (0 or 80)
# 25-49   EXTRA	extra space used for anything			
# 50-52	ORIGIN  origin in X,Y,Z used for transforms		
# 53	MAP	character string 'MAP ' to identify file type	
# 54	MACHST	machine stamp					
# 55	RMS	rms deviation of map from mean density		
# 56	NLABL	number of labels being used			
# 57-256 LABEL(20,10) 10 80-character text labels		


import sys, time, numpy as np


#--------------------------------------------------------------------------------------------------------
# CCP4 Writer
#--------------------------------------------------------------------------------------------------------

def write_ccp4_grid_data(volume,path):
	dtype = volume.dtype
	mtype = closest_ccp4_type(dtype)
	f = open(path, 'wb')

	header = ccp4_header(volume, mtype)
	f.write(header)

	I, J, K = volume.shape
	for k in range(K):
		matrix = volume[:,:,k]
		if dtype != mtype:
			matrix = matrix.astype(dtype)
		f.write(matrix.tostring())
	# Put matrix statistics in header
	header = ccp4_header(volume, dtype)
	f.seek(0)
	f.write(header)
	f.close()
	
#--------------------------------------------------------------------------------------------------------

def ccp4_header(volume, value_type):
	size = volume.shape
	from numpy import float32, int16, int8, int32
	if value_type == np.int8:
		mode = 0
	elif value_type == np.int16:
		mode = 1
	else:
		mode = 2

	cell_size = map(lambda a,b: a*b, (1,1,1), size)

	if np.little_endian:
		machst = 0x00004144
	else:
		machst = 0x11110000

	ver_stamp = '3DSRS %s' % time.asctime()
	labels = [ver_stamp[:80]]

	nlabl = len(labels)
	# Make ten 80 character labels.
	labels.extend(['']*(10-len(labels)))
	labels = [l + (80-len(l))*'\0' for l in labels]
	labelstr = ''.join(labels)

	dmin = volume.min()
	dmax = volume.max()
	dmean = volume.mean()

	strings = [
			binary_string(size, np.int32),  # nx, ny, nz
			binary_string(mode, np.int32),  # mode
			binary_string((0,0,0), np.int32), # nxstart, nystart, nzstart
			binary_string(size, np.int32),  # mx, my, mz
			binary_string(cell_size, np.float32), # cella
			binary_string((90,90,90), np.float32), # cellb
			binary_string((1,2,3), np.int32), # mapc, mapr, maps
			binary_string((dmin, dmax, dmean), np.float32), # dmin, dmax, dmean
			binary_string(0, np.int32), # ispg
			binary_string(0, np.int32), # nsymbt
			binary_string([0]*25, np.int32), # extra
			binary_string((0,0,0), np.float32), # origin
			'MAP ', # map
			binary_string(machst, np.int32), # machst
			binary_string(0, np.float32), # rms
			binary_string(nlabl, np.int32), # nlabl
			labelstr,
			]
	header = ''.join(strings)
	return header

#--------------------------------------------------------------------------------------------------------    

def binary_string(values, dtype):
	return np.array(values, dtype).tostring()

#--------------------------------------------------------------------------------------------------------

def closest_ccp4_type(dtype):
	if dtype in (np.float32, np.float64, np.float, np.int32, np.int, np.uint32, np.uint, np.uint16):
		ctype = np.float32
	elif dtype in (np.int16, np.uint8):
		ctype = np.int16
	elif dtype in (np.int8, np.int0, np.character):
		ctype = np.int8
	else:
		raise TypeError, ('Volume data has unrecognized type %s' % dtype)
	return ctype
    
#------------------------------------------------------------------------------   
    
def main():
	volume = np.load(sys.argv[1])
	fname = sys.argv[2]
	write_ccp4_grid_data(volume, fname)
	sys.exit()


if __name__ == "__main__":
    main()

