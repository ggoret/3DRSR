#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
# 3D Reciprocal Space Reconstruction
# Gael Goret and Alessandro Mirone for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr, mirone@esrf.fr
###################################################################################

import sys
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

class Projection(object):
	def __init__(self,theta, theta_offset,d1,d2,dist,lmbda,beam_tilt_angle,pixel_size,MD,det_origin_X,det_origin_Y,dim1,dim2,omega,phi,kappa,alpha,beta,omega_offset,r1,r2,r3):
		self.theta = theta
		self.theta_offset = theta_offset
		self.d1 = d1
		self.d2 = d2
		self.dist = dist
		
		self.beam_tilt_angle = beam_tilt_angle
		self.pixel_size = pixel_size
		self.MD = MD
		self.det_origin_X = det_origin_X
		self.det_origin_Y = det_origin_Y
		self.dim1 = dim1
		self.dim2 = dim2
		
		self.omega = omega
		self.omega_offset = omega_offset
		self.phi = phi
		self.kappa = kappa
		self.alpha = alpha
		self.beta = beta
		self.lmbda = lmbda 
		
		self.r1 = r1
		self.r2 = r2
		self.r3 = r3
		
	def set_theta(self,theta):
		self.theta = theta

	
	def Q0(self,X,Y):
		x,y = np.array((self.pixel_size*self.MD*np.array([[X-X0(self.det_origin_X,self.dim1)],[Y-Y0(self.det_origin_Y,self.dim2)]])).flatten())[0]
		p = P(self.theta, self.theta_offset,self.d1,self.d2,self.dist,x,y)
		q0 = (p/norm(p)+P0(self.dist,self.beam_tilt_angle)/self.dist)/self.lmbda
		return q0

		
	def Q(self,X,Y):
		q = R(self.omega,self.phi,self.kappa,self.alpha,self.beta,self.omega_offset).transpose().dot(self.Q0(X,Y))
		return q

		
	def Qfin(self,X,Y):
		qfin = U(self.r1,self.r2,self.r3).transpose().dot(self.Q(X,Y))
		return qfin
		
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

class Correction(object):
	def __init__(self,theta, theta_offset,d1,d2,dist,lmbda,beam_tilt_angle,pixel_size,MD,det_origin_X,det_origin_Y,dim1,dim2,pol_degree,n):
		self.theta = theta
		self.theta_offset = theta_offset
		self.d1 = d1
		self.d2 = d2
		self.dist = dist
		self.lmbda = lmbda
		self.beam_tilt_angle = beam_tilt_angle
		self.pixel_size = pixel_size
		self.MD = MD
		self.det_origin_X = det_origin_X
		self.det_origin_Y = det_origin_Y
		self.dim1 = dim1
		self.dim2 = dim2
		self.pol_degree = pol_degree # degree of polarisation
		self.n = n # normal to the polarisation plane
	
	def set_theta(self,theta):
		self.theta = theta

	
	def POL(self,X,Y,n):
		"""
		Polarisation correction
		n is the normal to the polarization plane (vector)
		"""
		x,y = np.array((self.pixel_size*self.MD*np.matrix([[X-X0(self.det_origin_X,self.dim1)],[Y-Y0(self.det_origin_Y,self.dim2)]])).flatten())[0]
		p = P(self.theta, self.theta_offset,self.d1,self.d2,self.dist,x,y) # p == p(X,Y)
		p0 = P0(self.dist,self.beam_tilt_angle)
		
		pol = self.pol_degree *((1-float(((np.cross(p0,self.n,axis=0).transpose()*p))/norm(np.cross(p0,n,axis=0))*norm(p))**2)+(1-self.pol_degree)*(1-((float(n.transpose()*p)/norm(p))**2)))
		return pol


	def C3(self,X,Y):
		x,y = np.array((self.pixel_size*self.MD*np.matrix([[X-X0(self.det_origin_X,self.dim1)],[Y-Y0(self.det_origin_Y,self.dim2)]])).flatten())[0]
		c3 = self.dist**3/(self.dist**2+x**2+y**2)**(3/2)
		return c3
		
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
		mag = norm(np.matrix(c))
		if mag > maxi:
			maxi = mag
	return maxi
		
#--------------------------------------------------------------------------------------------------------
# Main
#--------------------------------------------------------------------------------------------------------

def main():
	display_logo()
	
	dim1 = 32
	dim2 = 32

	data = np.ones((dim1,dim2),dtype = np.int32)*100000

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
		
	Q0 = np.zeros((dim2,dim1,3),dtype = np.float32) 
	MD_pix_tmp = params.pixel_size*MD

	Q_for_x = np.zeros((dim1 ,3),dtype = np.float32) 
	X_array_tmp = np.zeros( (dim1,2  ))
	X_array_tmp [:,0] = np.arange(dim1) - (params.det_origin_X +dim1/2.0 )

	Q_for_y = np.zeros((dim2 ,3),dtype = np.float32) 
	Y_array_tmp = np.zeros( (dim2,2  ))
	Y_array_tmp [:,1] = np.arange(dim2) - (params.det_origin_Y +dim2/2.0 )

	X_array_tmp=np.tensordot( MD_pix_tmp ,  X_array_tmp , axes=([1],[1])   )
	Y_array_tmp=np.tensordot( MD_pix_tmp ,  Y_array_tmp , axes=([1],[1])   )
	
	YX_array_tmp =  Y_array_tmp[:,None,:] + X_array_tmp[None,:,:] 
	
	print YX_array_tmp.shape
	P_total_tmp = np.zeros((dim2, dim1 ,3),dtype = np.float32)
	P_total_tmp[:,:, 1:3] =  YX_array_tmp 
	P_total_tmp[:,:,  0 ] =  -DIST
	P_total_tmp=np.tensordot( DET  , P_total_tmp  , [1],[1]   )
	P_total_tmp_modulus = np.sqrt( np.sum(  P_total_tmp   , axis=-1))
	
	
	Q0 = ((P_total_tmp.T/P_total_tmp_modulus.T).T - p0/params.dist)/params.lmbda
	
	print Q0

	"""

	X-X0(self.det_origin_X,self.dim1)


	x,y = np.array((*np.array([[X-X0(self.det_origin_X,self.dim1)],[Y-Y0(self.det_origin_Y,self.dim2)]])).flatten())[0]
	p = P(self.theta, self.theta_offset,self.d1,self.d2,self.dist,x,y)
	q0 = (p/norm(p)+P0(self.dist,self.beam_tilt_angle)/self.dist)/self.lmbda
	


		for Y in range(dim2):
			Q0[X,Y,:] = Proj.Q0(X,Y).flatten()


	for X in range(dim1):
		for Y in range(dim2):
			Q0[X,Y,:] = Proj.Q0(X,Y).flatten()
	print 'Q0(X,Y) :\n',Q0
	
	# Construction 3D intensity distribution
	# --------------------------------------
	
	pol_degree = 1. # polarisation degree
	n = np.matrix([[0],[0],[1]]) # normal to polarisation plane
	
	Corr = Correction(theta, theta_offset,d1,d2,params.dist,params.lmbda,params.beam_tilt_angle,params.pixel_size,MD,params.det_origin_X,params.det_origin_Y,dim1,dim2,pol_degree,n)
	
	# Estimation of the cube size :
	corners  = [Q0[0,0,:],Q0[0,dim2-1,:],Q0[dim1-1,0,:],Q0[dim1-1,dim2-1,:]]
	
	Qmax = mag_max(corners) # maximal magnitude for reciprocal vector corresponding to the corners pixels
	
	print 'Qmax : ', Qmax
	
	cube_dim = sup_pow_2(2*Qmax) + 1 # closest power of 2 > 2Qmax + 1 (to be sure that we have a symetric center)
	
	print 'CUBE DIMs :', cube_dim
	
	A = np.zeros((cube_dim,cube_dim,cube_dim),dtype = np.float32)
	B = np.zeros((cube_dim,cube_dim,cube_dim), dtype = np.int32)
	
	dqx = dqy = dqz = cube_dim//2
		
	q0x, q0y, q0z = 0,0,0	
	
	print 'RECIPROCAL SPACE CENTER :', q0x, q0y, q0z
	print '-----------------------------'
	
	ent = 0
	for X in range(dim1):
		for Y in range(dim2):
			qx,qy,qz = np.array(Proj.Qfin(X,Y).flatten())[0]
			#print qx,qy,qz
			print qx,qy,qz
			i ,j ,k = np.floor(np.array([qx+dqx, qy+dqy, qz+dqz]))
			#print i,j,k
			if i <= dqx and j <= dqy and k <= dqz:
				ent +=1
				A[i,j,k] = data[X,Y]/(Corr.POL(X,Y,n)*Corr.C3(X,Y))
	print '3D Intensity Distribution :'	
	print '%d values'%ent	
	print A
	
	"""
		
	sys.exit()


if __name__ == "__main__":
    main()
    
    
    
