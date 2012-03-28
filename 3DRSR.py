#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
# 3D Reciprocal Space Reconstruction
# Gael Goret and Alessandro Mirone for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr, mirone@esrf.fr
###################################################################################

import fabio
import sys, time
import numpy as np
from numpy.linalg import norm

np.matrix=None

__version__=0.0

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
# Rotation Matrix
#--------------------------------------------------------------------------------------------------------


def Rotation(angle, rotation_axis=0 ):
	if type(rotation_axis)==type("") :
		rotation_axis={"x":0,"y":1,"z":2}
	assert((rotation_axis>=0 and rotation_axis<=3))
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
	
def R(omega,phi,kappa,alpha,beta,omega_offset):
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

def U(r1,r2,r3):
	""" Secondary rotation of reciprocal space (to orient the crystallographic axis in a special way) """
	u = Rotation(r3,2).dot(Rotation(r2,1)).dot(Rotation(r1,0))
	return u

#--------------------------------------------------------------------------------------------------------
# Projection and Orientation
#--------------------------------------------------------------------------------------------------------

def P0(dist,b2):
	""" Beam tilt matrix """
	B = Rotation(b2,1)
	""" Primary projection of pixel coordinates (X,Y) to the reciprocal space. """
	# p0 = BM.dot(np.matrix([[dist],[0],[0]]))
	p0 = np.dot(B,   [ dist, 0,0 ] ) 
	return p0
	
#--------------------------------------------------------------------------------------------------------

def P(theta, theta_offset,d1,d2,dist,x,y):
	res = np.dot(DET(theta, theta_offset, d1,d2) , [-dist , x , y])
	return res

#--------------------------------------------------------------------------------------------------------
# Corrections
#--------------------------------------------------------------------------------------------------------

MD0 = np.array([[1,0],[0,1]  ],"d")
MD1 = np.array([[-1,0],[0,1] ],"d")
MD2 = np.array([[1,0],[0,-1] ],"d") 
MD3 = np.array([[-1,0],[0,-1]],"d")
MD4 = np.array([[0,1],[1,0]  ],"d")
MD5 = np.array([[0,-1],[1,0] ],"d")
MD6 = np.array([[0,1],[-1,0] ],"d")
MD7 = np.array([[0,-1],[-1,0]],"d")

#--------------------------------------------------------------------------------------------------------
# Xcalibur parameter Class
#--------------------------------------------------------------------------------------------------------

class XCalibur_parameters(object):
	def __init__(self):
		self.pixel_size = 0.172 # mm
		detector_distance = 174.42 # mm	
		self.dist = detector_distance * 1.72  # mm
		self.beam_tilt_angle = 0.99075 # deg
		self.det_origin_X =  1733.19127
		self.det_origin_Y =  1712.50841
		self.lmbda = 0.67018 # Wavelength users specified (ang)
		
#--------------------------------------------------------------------------------------------------------
# Miscelaneous
#--------------------------------------------------------------------------------------------------------	          

def sup_pow_2(x):
	p=1
	while(p<x):
		p=p*2
	return p
	
#--------------------------------------------------------------------------------------------------------

def mag_max(l):
	"""
	l is a list of coordinated 
	"""
	maxi = 0
	for c in l:
		mag = norm(c)
		if mag > maxi:
			maxi = mag
	return maxi
		
#--------------------------------------------------------------------------------------------------------
# Main
#--------------------------------------------------------------------------------------------------------

def main():
	display_logo()
	time0 = time.time()
	print 'Reading File ...'
	"""
	img = fabio.open('feo1_1_00001.cbf')
	time1 = time.time()
	print '-> t = %.2f s'%(time1-time0)
	data = img.data
	dim1,dim2 = img.dim1,img.dim2
	"""
	dim1,dim2 = 1200,1024
	print 'Image Dimension :',dim2,dim1
	data = np.ones((dim2,dim1),dtype = np.int32)*10000
	

	print 'Setting Parameters ...'
	params = XCalibur_parameters()

	#PREPARATION STEP 1
	p0 = P0(params.dist,params.beam_tilt_angle)
	MD = MD0

	#PREPARATION STEP 2-3-4
	kappa = -134.
	alpha = 50.
	beta = 50.
	omega = 57.
	theta = 0.
	phi = 0. #(n-1)*0.1 where n is the number of image
	
	omega_offset = -0.19777
	theta_offset = 0.39804
	
	d1 = -0.41144
	d2 = 1.17097
	
	r1 = -88.788
	r2 = 2.257
	r3 = 69.629
	time2 = time.time()
	print 'Computation of Initial Projection Coordinates Q0'
	Q0 = np.zeros((dim2,dim1,3),dtype = np.float32) 
	MD_pix_tmp = params.pixel_size*MD

	X_array_tmp = np.zeros( (dim1,2  ))
	X_array_tmp [:,0] = np.arange(dim1) - (params.det_origin_X +dim1/2.0 )

	Y_array_tmp = np.zeros( (dim2,2  ))
	Y_array_tmp [:,1] = np.arange(dim2) - (params.det_origin_Y +dim2/2.0 )

	X_array_tmp=np.tensordot(   X_array_tmp , MD_pix_tmp ,axes=([1],[1])   )
	Y_array_tmp=np.tensordot(   Y_array_tmp ,MD_pix_tmp , axes=([1],[1])   )
	
	YX_array_tmp =  Y_array_tmp[:,None,:] + X_array_tmp[None,:,:] 
	
	P_total_tmp = np.zeros((dim2, dim1 ,3),dtype = np.float32)
	P_total_tmp[:,:,1:3] = YX_array_tmp 
	P_total_tmp[:,:, 0 ] = -params.dist
	
	P_total_tmp=np.tensordot(P_total_tmp, DET(theta, theta_offset, d1,d2),   axes=([2],[1]))

	P_total_tmp_modulus = norm( P_total_tmp)
	
	
	
	Q0 = ((P_total_tmp/P_total_tmp_modulus) - p0/params.dist)/params.lmbda
	
	time3 = time.time()
	print '-> t = %.2f s'%(time3-time2)	
	
	
	Prim_Rot_of_RS = R(omega,phi,kappa,alpha,beta,omega_offset)
	Snd_Rot_of_RS = U(r1,r2,r3)
	
	time4 = time.time()
	print 'Primary Rotation : Q0 -> Q'
	Q = np.tensordot (  Q0 , Prim_Rot_of_RS.T , axes=([2],[1]))
	time5 = time.time()
	print '-> t = %.2f s'%(time5-time4)
	print 'Secondary Rotation : Q -> Qfin'
	Qfin = np.tensordot (Q , Snd_Rot_of_RS.T , axes=([2],[1]))
	time6 = time.time()
	print '-> t = %.2f s'%(time6-time5)
	print Qfin.shape
	print '--------------------------------------'
	
	pol_degree = 1. # polarisation degree
	normal_to_pol = np.array([0,0,1]) # normal to polarisation plane
	
	#CORRECTION OPTIONS
	cpt_pol = False 
	cpt_c3 = False
	
	if cpt_pol :
		print 'Computation of Polarisation Correction ...'
		P0xn = np.cross(p0,normal_to_pol,axis=0)
		NormP0xn =  norm(P0xn)
		# np.tensordot(P0xn,P_total_tmp,axes=([0],[2]) 
		POL_tmp = pol_degree*(1-( (P0xn*P_total_tmp).sum(axis=-1))/(NormP0xn*P_total_tmp_modulus))**2	
	
		POL_tmp += 	(1-pol_degree)*(1-(np.tensordot(normal_to_pol,P_total_tmp,axes=([0],[2]))/P_total_tmp_modulus)**2)
		print 'POL_tmp.shape',POL_tmp.shape
	else :
		print 'Computation of Polarisation Correction : Canceled'
		#POL_tmp = np.ones((dim2,dim1),dtype = np.float32)
	
	time7 = time.time()
	print '-> t = %.2f s'%(time7-time6)	
	
	if cpt_c3:
		print 'Computation of Flux Density and Parallax Correction ...'
		C3 =  params.dist**3/(params.dist**2+np.sum (YX_array_tmp*YX_array_tmp, axis = -1 ) )**(3/2)
		print 'C3.shape',C3.shape
	else:
		print 'Computation of Flux Density and Parallax Correction : Canceled'
		#C3 = np.ones((dim2,dim1),dtype = np.float32)
		
	time8 = time.time()
	print '-> t = %.2f s'%(time8-time7)
	print 'Estimation of Qmax ...'
	# Estimation of Qmax
	corners  = np.array([Q0[0,0,:],Q0[0,dim1-1,:],Q0[dim2-1,0,:],Q0[dim2-1,dim1-1,:]])
	Qmax = mag_max(corners) # maximal magnitude for reciprocal vector corresponding to the corners pixels
	print '-----------------------------'
	print 'Qmax = ', Qmax
	
	#cube_dim = sup_pow_2(2*Qmax) + 1 # closest power of 2 > 2Qmax + 1 (to be sure that we have a symetric center)

	cube_dim=128 + 1

	print 'CUBE DIMs =', cube_dim
	
	dqx = dqy = dqz = Qmax
		
	q0x = q0y = q0z = 0
	
	print 'RECIPROCAL SPACE CENTER  =', q0x, q0y, q0z
	print '-----------------------------'
	
	print 'Computation of 3D Volume Indices ...'
	QxFilter = np.where((Qfin[:,:,0]-q0x) <= dqx)
	QyFilter = np.where((Qfin[:,:,1]-q0y) <= dqy)
	QzFilter = np.where((Qfin[:,:,2]-q0z) <= dqz)

	Qxfin = np.zeros((dim2,dim1),dtype = np.float32)
	Qyfin = np.zeros((dim2,dim1),dtype = np.float32)
	Qzfin = np.zeros((dim2,dim1),dtype = np.float32)
	
	Qxfin[QxFilter] = Qfin[QxFilter[0],QxFilter[1],0]
	Qyfin[QyFilter] = Qfin[QyFilter[0],QyFilter[1],1]
	Qzfin[QzFilter] = Qfin[QzFilter[0],QzFilter[1],2]
	
	time9 = time.time()
	I_array = ( np.floor ( np.sqrt(cube_dim-1) * (1 + (Qxfin - q0x)/dqx) ) ).astype(np.int32) 
	J_array = ( np.floor ( np.sqrt(cube_dim-1) * (1 + (Qyfin - q0y)/dqy) ) ).astype(np.int32)
	K_array = ( np.floor ( np.sqrt(cube_dim-1) * (1 + (Qzfin - q0z)/dqz) ) ).astype(np.int32)
	
	time10 = time.time()
	print '-> t = %.2f s'%(time10-time9)
	
	print '-------------------'
	print I_array
	print '-------------------'
	print J_array
	print '-------------------'
	print K_array
	print '-------------------'
	
	Volume = np.zeros((cube_dim,cube_dim,cube_dim),dtype = np.float32)
	
	print 'Filling up the Volume with Corrected Intensity ...'
	#Intensity = (data/np.tensordot(POL_tmp,C3,axes=([1],[1])))#Wrong
	Intensity = data
	print 'Intensity.shape',Intensity.shape
	print 'Volume.shape',Volume.shape
	Volume[I_array,J_array,K_array] = Intensity 
	
	time11 = time.time()
	print '-> t = %.2f s'%(time11-time10)
	
	print '3D Intensity Distribution :'
	print np.where(Volume>0)
		
	sys.exit()


if __name__ == "__main__":
    main()
    
    
    
