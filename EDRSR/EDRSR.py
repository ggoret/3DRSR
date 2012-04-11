#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
# 3D Reciprocal Space Reconstruction
# Gael Goret and Alessandro Mirone for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr, mirone@esrf.fr
###################################################################################
try :
	import fabio
except :
	print 'Warning : fabio module could not be initialised'
try:
	raise Exception
	from mayavi import mlab
except:
	print 'Warning : mayavi module could not be initialised'
try:
	import scipy as np
except:
	import numpy as np
import sys, time
__version__=0.1

import fillvolume

#--------------------------------------------------------------------------------------------------------
# Logo
#--------------------------------------------------------------------------------------------------------

def display_logo():
	print "                                                            " 
	print "                             ___          ___          ___   "
	print "                            /\  \        /\  \        /\  \   "
	print "                           /  \  \      /  \  \      /  \  \   "
	print "__________________________/ /\ \  \____/ /\ \  \____/ /\ \  \__ " 
	print "\       __     _____     /  \~\ \  \  _\ \~\ \  \  /  \~\ \  \ \ "
	print " \    /'__`\  /\  _ `\  / /\ \ \ \__\/\ \ \ \ \__\/ /\ \ \ \__\ \ " 
	print "  \  /\_\L\ \ \ \ \/\ \ \/_|  \/ /  /\ \ \ \ \/__/\/_|  \/ /  /  \ "
	print "   \ \/_/_\_<_ \ \ \ \ \   | |  /  /  \ \ \ \__\     | |  /  /    \ "
	print "    \  /\ \L\ \ \ \ \_\ \  | |\/__/    \ \/ /  /     | |\/__/      \ "
	print "     \ \ \____/  \ \____/  | |  |       \  /  /      | |  |         \ "
	print "      \ \/___/    \/___/    \|__|        \/__/        \|__|   v %3.1f  \ "%__version__
	print "       \______________________________________________________________\ "
	print "                                                                         "
	print "      3DRSR : a software for volumetric Reconstruction of Reciprocal Space "
	print ""

#--------------------------------------------------------------------------------------------------------
# Rotation, Projection and Orientation Matrix
#--------------------------------------------------------------------------------------------------------


def Rotation(angle, rotation_axis=0 ):
	if type(rotation_axis)==type("") :
		rotation_axis={"x":0,"y":1,"z":2}
	assert((rotation_axis>=0 and rotation_axis<=3))
	angle = np.radians(angle)
	ret_val = np.zeros([3,3],"d")
	i1=rotation_axis
	i2=(rotation_axis+1)%3
	i3=(rotation_axis+2)%3
	ret_val[i1,i1  ] =  1
	ret_val[i2,i2  ] =  np.cos(angle)
	ret_val[i3,i3  ] =  np.cos(angle)
	ret_val[i2,i3  ] =  np.sin(angle)
	ret_val[i3,i2  ] = -np.sin(angle)
	return ret_val

#--------------------------------------------------------------------------------------------------------	
	
def Prim_Rot_of_RS(omega,phi,kappa,alpha,beta,omega_offset):
	""" 
	Primary rotation of reciprocal space. 
	Omega, kappa, phi are the nominal values of the angle
	"""
	tmp = Rotation(omega + omega_offset,2)
	tmp = np.dot(tmp, Rotation(alpha,1))
	tmp = np.dot(tmp, Rotation(kappa,2))
	tmp = np.dot(tmp, Rotation(-alpha,1))
	tmp = np.dot(tmp, Rotation(beta,1))
	tmp = np.dot(tmp, Rotation(phi,2))
	tmp = np.dot(tmp, Rotation(-beta,1))
	return tmp

	r = R3(omega + omega_offset).dot(R2(alpha)).dot(R3(kappa)).dot(R2(-alpha)).dot(R2(beta)).dot(R3(phi)).dot(R2(-beta))
	return r

#--------------------------------------------------------------------------------------------------------
	
def DET(theta, theta_offset, d1,d2):
	""" 
	Rotation matrix for the detector 
	theta is the nominal theta value and d1,D2 are the tilts of detector
	"""
	det = Rotation(theta + theta_offset,2).dot(Rotation(d2,1)).dot(Rotation(d1,0))
	return det
#--------------------------------------------------------------------------------------------------------	

def Snd_Rot_of_RS(r1,r2,r3):
	""" Secondary rotation of reciprocal space (to orient the crystallographic axis in a special way) """
	u = Rotation(r3,2).dot(Rotation(r2,1)).dot(Rotation(r1,0))
	return u

#--------------------------------------------------------------------------------------------------------

def P0(dist,b2):
	""" Primary projection of pixel coordinates (X,Y) to the reciprocal space. """
	B = Rotation(b2,1) # Beam tilt matrix
	p0 = np.dot(B,   [ dist, 0,0 ] ) 
	return p0
	
#--------------------------------------------------------------------------------------------------------
# Corrections
#--------------------------------------------------------------------------------------------------------

MD0 = np.array([[1,0],[0,1]  ],dtype = np.int32)
MD1 = np.array([[-1,0],[0,1] ],dtype = np.int32)
MD2 = np.array([[1,0],[0,-1] ],dtype = np.int32)
MD3 = np.array([[-1,0],[0,-1]],dtype = np.int32)
MD4 = np.array([[0,1],[1,0]  ],dtype = np.int32)
MD5 = np.array([[0,-1],[1,0] ],dtype = np.int32)
MD6 = np.array([[0,1],[-1,0] ],dtype = np.int32)
MD7 = np.array([[0,-1],[-1,0]],dtype = np.int32)

#--------------------------------------------------------------------------------------------------------
# Parameter Class
#--------------------------------------------------------------------------------------------------------

class Parameters(object):
	def __init__(self):
		self.pixel_size = 0.172 # mm
		detector_distance = 174.42 # mm	
		self.dist = detector_distance * 1.72  # mm
		self.beam_tilt_angle = 0.99075 # deg

		self.det_origin_X =  1733.19127
		self.det_origin_Y =  1712.50841
		
		self.lmbda = 0.67018 # Wavelength users specified (ang)
		self.kappa = -134.
		self.alpha = 50.
		self.beta = 0.
		self.omega = 57.
		self.theta = 0.
		self.phi = 0. #(n-1)*0.1 where n is the number of image
	
		self.omega_offset = -0.19777
		self.theta_offset = 0.39804
	
		self.d1 = -0.41144
		self.d2 = 1.17097
	
		self.r1 = -88.788
		self.r2 = 2.257
		self.r3 = 69.629
		
		self.pol_degree = 1. # polarisation degree
		self.normal_to_pol = np.array([0,0,1]) # normal to polarisation plane
		
		self.cube_dim = 256 + 1
		
		#DEBUG OPTIONS :
		self.cpt_corr = True

#--------------------------------------------------------------------------------------------------------
# Miscelaneous
#--------------------------------------------------------------------------------------------------------	          

def mag_max(l):
	"""
	l is a list of coordinates
	"""
	return max(np.sqrt(np.sum(l*l,axis=-1)))

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

	return 

#--------------------------------------------------------------------------------------------------------
# Projection
#--------------------------------------------------------------------------------------------------------
"""
def project_image(data,                            float32 (dim1,dim2) 
		  p0,                              float(dim1,dim2,3)
		  Q0,                              float(dim1,dim2,3)
		  XY_array_tmp,                    float(dim1,dim2,2) 
		  P_total_tmp,                     float(dim1,dim2,3)
		  P_total_tmp_modulus,             float(dim1,dim2)
		  Qmax,                            float
		  params,                          object
		  Volume,                          float(nz,ny,nz   )
		  Mask,                            float(nz,ny,nz   ) 
		  C3,                              float32 (dim1,dim2) 
		  POL_tmp                          float32 (dim1,dim2)  
		  ):
"""

def project_image(data,p0,Q0,XY_array_tmp,P_total_tmp,P_total_tmp_modulus,Qmax,params, Volume, Mask,C3,  POL_tmp):
	
	dim1,dim2 = data.shape


	R = Prim_Rot_of_RS(params.omega,params.phi,params.kappa,params.alpha,params.beta,params.omega_offset)
	U = Snd_Rot_of_RS(params.r1,params.r2,params.r3)
	
	time4 = time.time()
	print 'Primary Rotation : Q0 -> Q'
	Q = np.tensordot (  Q0 , R.T , axes=([2],[1]))
	time5 = time.time()
	print '-> t = %.2f s'%(time5-time4)
	print 'Secondary Rotation : Q -> Qfin'
	Qfin = np.tensordot (Q , U.T , axes=([2],[1]))
	time6 = time.time()
	print '-> t = %.2f s'%(time6-time5)
	print Qfin.shape
	print '--------------------------------------'
	
	cube_dim= params.cube_dim

	print 'CUBE DIMs =', cube_dim

	dqx = dqy = dqz = Qmax
		
	q0x = q0y = q0z = 0
	
	comment="""
	# Qfin float (dim1,dim2,3)
	#  Volume,                          float(nz,ny,nz   )
	#  Mask,                            float(nz,ny,nz   ) 
	#  C3,                              float32 (dim1,dim2) 
	#  POL_tmp                          float32 (dim1,dim2)  
	# (data,                            float32 (dim1,dim2) 
"""
	print " COMMENT " ,  comment
	print  "Qfin    ",  Qfin   .shape            
	print  "Volume  ",  Volume .shape          
	print  "Mask    ",  Mask   .shape         
	print  "C3      ",  C3     .shape         
	print  "POL_tmp ",  POL_tmp.shape            
	print  "data    ",  data   .shape   

	print "  ENTREE DE FUNC_SOMME >>>>>>>>>>>>>>>>>>>>>>>>>>" 
	Qfin=Qfin.astype(np.float32)
	print Volume.dtype
	print Mask.dtype
	print Qfin.dtype
	print " ----------- " 
	print data.dtype
	print POL_tmp.dtype
	print C3.dtype
	#func_somme(cube_dim  , q0x , q0y , q0z,  dqx , dqy , dqz , Qfin,Volume,Mask, data,POL_tmp,C3  )
	fillvolume.func_somme(q0x, q0y, q0z, dqx, dqy, dqz, Volume, Mask,  Qfin, data, POL_tmp, C3)
	print "  SORTIE DE FUNC_SOMME <<<<<<<<<<<<<<<<<<<<<<<<<<" 
	
def func_somme(cube_dim  , q0x , q0y , q0z,  dqx , dqy , dqz , Qfin,Volume,Mask, data, POL_tmp, C3):
	dim1,dim2 = data.shape

	print 'RECIPROCAL SPACE CENTER  =', q0x, q0y, q0z
	print '-----------------------------'
	
	print 'Computation of 3D Volume Indices ...'
	
	time9 = time.time()
	I_array_tmp = (np.floor((cube_dim-1)//2 * (1 + (Qfin[:,:,0]-q0x)/dqx))).astype(np.int32)
	J_array_tmp = (np.floor((cube_dim-1)//2 * (1 + (Qfin[:,:,1]-q0y)/dqy))).astype(np.int32)
	K_array_tmp = (np.floor((cube_dim-1)//2 * (1 + (Qfin[:,:,2]-q0z)/dqz))).astype(np.int32)
	
	I_array = I_array_tmp[np.where(I_array_tmp<= cube_dim)].reshape(dim1,dim2)
	J_array = J_array_tmp[np.where(J_array_tmp<= cube_dim)].reshape(dim1,dim2)
	K_array = K_array_tmp[np.where(K_array_tmp<= cube_dim)].reshape(dim1,dim2)
	time10 = time.time()
	print '-> t = %.2f s'%(time10-time9)
	

	print 'Summing Corrected Intensity into the Volume ...'
	Volume[I_array,J_array,K_array] += data/(POL_tmp*C3) # Data Correction
	Mask[I_array,J_array,K_array] += 1

	return I_array,J_array,K_array
	
#--------------------------------------------------------------------------------------------------------
# Main
#--------------------------------------------------------------------------------------------------------

def main():
	display_logo()
	flist = sys.argv[1:]
	time0 = time.time()
	print 'Reading File ...'
	
	#DEBUG OPTIONS
	rendering = False
	file_saving = True
	
	img = fabio.open(flist[0])
	time1 = time.time()
	print '-> t = %.2f s'%(time1-time0)
	data = img.data
	dim1,dim2 = data.shape
	
	"""
	dim1,dim2 = 2527,2463 
	data = np.ones((dim1,dim2),dtype = np.int32)*10000
	"""
	print 'Image Dimension :',dim1,dim2

	print 'Setting Parameters ...'
	params = Parameters()

	p0 = P0(params.dist,params.beam_tilt_angle)
	MD = MD0
	
	time2 = time.time()
	print 'Computation of Initial Projection Coordinates Q0'
	Q0 = np.zeros((dim1,dim2,3),dtype = np.float32) 
	MD_pix_tmp = params.pixel_size*MD

	X_array_tmp = np.zeros( (dim1,2  ))
	X_array_tmp [:,0] = np.arange(dim1) - (params.det_origin_X - 1725 + dim1/2.0 )

	Y_array_tmp = np.zeros( (dim2,2  ))
	Y_array_tmp [:,1] = np.arange(dim2) - (params.det_origin_Y - 1725 + dim2/2.0 )

	X_array_tmp=np.tensordot(   X_array_tmp , MD_pix_tmp ,axes=([1],[1])   )
	Y_array_tmp=np.tensordot(   Y_array_tmp ,MD_pix_tmp , axes=([1],[1])   )
	
	XY_array_tmp =  X_array_tmp[:,None,:] + Y_array_tmp[None,:,:]  
	
	P_total_tmp = np.zeros((dim1, dim2 ,3),dtype = np.float32)
	P_total_tmp[:,:, 0 ] = -params.dist
	P_total_tmp[:,:,1:3] = XY_array_tmp 
	
	P_total_tmp=np.tensordot(P_total_tmp, DET(params.theta, params.theta_offset, params.d1,params.d2),   axes=([2],[1]))
	
	P_total_tmp_modulus = np.sqrt(np.sum(P_total_tmp*P_total_tmp,axis=-1))
	
	Q0_tmp = P_total_tmp.T/P_total_tmp_modulus.T
	Q0 = (Q0_tmp.T  + p0/params.dist)/params.lmbda
	
	time3 = time.time()
	print '-> t = %.2f s'%(time3-time2)	
	
	time6 = time.time()
	if params.cpt_corr :
		print 'Computation of Polarisation Correction ...'
		P0xn = np.cross(p0,params.normal_to_pol,axis=0)
		P0xn_modulus =  np.sqrt(np.sum(P0xn*P0xn,axis = -1))
		# np.tensordot(P0xn,P_total_tmp,axes=([0],[2]) 
		POL_tmp = params.pol_degree*(1-( (P0xn*P_total_tmp).sum(axis=-1))/(P0xn_modulus*P_total_tmp_modulus))**2	
	
		POL_tmp += 	(1-params.pol_degree)*(1-(np.tensordot(params.normal_to_pol,P_total_tmp,axes=([0],[2]))/P_total_tmp_modulus)**2)
	else :
		print 'Computation of Polarisation Correction : Canceled'

	
	POL_tmp=POL_tmp.astype(np.float32)


	time7 = time.time()
	print '-> t = %.2f s'%(time7-time6)	
	
	if params.cpt_corr:
		print 'Computation of Flux Density and Parallax Correction ...'
		C3 = ( params.dist**3/(params.dist**2+np.sum (XY_array_tmp*XY_array_tmp, axis = -1 ) )**(3/2)).astype(np.float32)
	else:
		print 'Computation of Flux Density and Parallax Correction : Canceled'
		
	time8 = time.time()
	print '-> t = %.2f s'%(time8-time7)
	
	print 'Estimation of Qmax ...'
	# Estimation of Qmax
	corners  = np.array([Q0[0,0,:],Q0[0,dim2-1,:],Q0[dim1-1,0,:],Q0[dim1-1,dim2-1,:]])
	Qmax = mag_max(corners) # maximal magnitude for reciprocal vector corresponding to the corners pixels
	print '-----------------------------'
	print 'Qmax = ', Qmax
	
	cube_dim= params.cube_dim
	Volume = np.zeros((cube_dim,cube_dim,cube_dim),dtype = np.float32)
	Mask = np.zeros((cube_dim,cube_dim,cube_dim),dtype = np.float32)
	total = len(flist)
	nbfile = 0.
	for fname in flist:
		
		timeI0 = time.time()
		
		img = fabio.open(fname)
		print 'Working on image %s'%fname
		data = img.data.astype(np.float32)

		"""
		dim1,dim2 = 2527,2463 
		data = np.ones((dim1,dim2),dtype = np.int32)*10000
		"""
		params.phi = nbfile*0.1
		# I_array,J_array,K_array = 
		project_image(data,  p0, Q0, XY_array_tmp, P_total_tmp, P_total_tmp_modulus, Qmax, params,
							Volume, Mask, C3, POL_tmp)
		# print 'Summing Corrected Intensity into the Volume ...'
		# Volume[I_array,J_array,K_array] += data/(POL_tmp*C3) # Data Correction
		# Mask[I_array,J_array,K_array] += 1
		nbfile += 1
		print '##################################'
		print 'Progression : %6.2f %% '%((nbfile/total)*100.)
		timeI1 = time.time()	
		print '-> time for this image = %.2f s'%(timeI1-timeI0)
		print '##################################'
		
	time11 = time.time()
	print '3D Intensity Distribution : Done'
	threeDid = Volume[np.where(Volume>0)]
	print threeDid
	vmin,vmax = min(threeDid),max(threeDid)
	print '-> Total time = %.2f s'%(time11-time0)
	if file_saving:
		np.save('intensity_distribution',Volume)
		np.save('mask',Mask)
	if rendering:
		print 'Rendering...'
		s = mlab.pipeline.volume(mlab.pipeline.scalar_field(Volume),vmin=-1, vmax=vmax+1)
		s.scene.background = (1,1,1)
		mlab.show()
	print 'Normal END'
	


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print('Usage : python 3DRSR.py image-file(s)')
	else:
		main()
	sys.exit()
    
    
