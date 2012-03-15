#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
# 3D Reciprocal Space Reconstruction
# Gael Goret and Alessandro Mirone for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr, mirone@esrf.fr
###################################################################################

import sys
import numpy as np 

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

def R1(x):
	r1 = np.matrix([[1, 0,0], [0,np.cos(x),np.sin(x)],[0,-np.sin(x),np.cos(x)]])
	return r1

def R2(x):
	r2 = np.matrix([[np.cos(x), 0,-np.sin(x)], [0,1,0],[np.cos(x),0,np.cos(x)]])
	return r2
	
def R3(x):
	r3 = np.matrix([[np.cos(x),np.sin(x) ,0], [-np.sin(x),np.cos(x),0],[0,0,1]])
	return r3
	
def R(omega,phi,kappa,alpha,beta,omega_offset):
	""" 
	Primary rotation of reciprocal space. 
	Omega, kappa, phi are the nominal values of the angle
	"""
	r = R3(omega + omega_offset).dot(R2(alpha)).dot(R3(kappa)).dot(R2(-alpha)).dot(R2(beta)).dot(R3(phi)).dot(R2(beta))
	return r
	
def DET(theta, theta_offset, d1,d2):
	""" 
	Rotation matrix for the detector 
	theta is the nominal theta value and d1,D2 are the tilts of detector
	"""
	det = R3(theta + theta_offset).dot(R2(d2)).dot(d1)
	return det
	
def U(r1,r2,r3):
	""" Secondary rotation of reciprocal space (to orient the crystallographic axis in a special way) """
	u = R3(r3).dot(R2(r2)).dot(R1(r1))
	return u

#--------------------------------------------------------------------------------------------------------
# Projection and Orientation
#--------------------------------------------------------------------------------------------------------

def B(b2):
	""" Beam tilt matrix """
	b = R2(b2)
	return b

def P0(dist,b2):
	""" Primary projection of pixel coordinates (X,Y) to the reciprocal space. """
	p0 = B(b2).dot(np.matrix([[dist],[0],[0]]))
	return p0
	
def X0(X,vertical_size):
	x0 = X - 1725 + vertical_size
	return x0
	
def Y0(Y,horizontal_size):
	y0 =  Y - 1725 + horizontal_size
	return y0
	
MD0 = np.matrix([[1,0],[0,1]])
MD1 = np.matrix([[-1,0],[0,1]])
MD2 = np.matrix([[1,0],[0,-1]]) 
MD3 = np.matrix([[-1,0],[0,-1]])
MD4 = np.matrix([[0,1],[1,0]])
MD5 = np.matrix([[0,-1],[1,0]])
MD6 = np.matrix([[0,1],[-1,0]])
MD7 = np.matrix([[0,-1],[-1,0]])

#--------------------------------------------------------------------------------------------------------
# Corrections
#--------------------------------------------------------------------------------------------------------

def POL(X,Y,n):
	"""
	Polarisation correction
	n is the normal to the polarization plane (vector)
	"""
	#TODO
	return
#--------------------------------------------------------------------------------------------------------
# Xcalibur parameter Class
#--------------------------------------------------------------------------------------------------------

class XCalibur_parameters(object):
	def __init__(self):
		self.pixel_size = 0.172 # mm
		detector_distance = 174.42 # mm	
		self.dist = detector_distance * 1.72  # mm
		self.beam_tilt_angle = 0.99075 # deg
#--------------------------------------------------------------------------------------------------------
# Main
#--------------------------------------------------------------------------------------------------------	          

def main():
	display_logo()
	
	data = np.ones((dim1,dim2),dtype = np.int32) 
	
	dim1 = 8
	dim2 = 8

	params = XCalibur_parameters()

	#STEP 1
	p0 = P0(params.dist,params.beam_tilt_angle)
	MD = MD0

	#STEP 2
	coords = np.zeros((dim1,dim2,2),dtype = np.float32)
	for X in range(dim1):
		for Y in range(dim2):
			coords[X,Y,:] = (params.pixel_size*MD0*np.matrix([[X0(X,dim1)],[Y0(Y,dim2)]])).flatten()
	
	print coords
	kappa = -134
	omega = 57
	theta = 0
	phi = 0
	
	omega_offset = omega
	theta_offset = theta
	
	d1 = -0.41144
	d2 = 1.17097
	
	lmbda = 0.67018 # ang
	
	#STEP 3
	p = lambda x,y : DET(theta, theta_offset, d1,d2).dot(np.matrix([[-params.dist],[x],[y]]))
	
	#STEP 4
	Q0 = lambda X,Y : ((p(coords[X,Y,0],coords[X,Y,1]))/np.abs(p(coords[X,Y,0],coords[X,Y,1]))+p0/params.dist)/lmbda
	 
	
	Q0XY = np.zeros((dim1,dim2,3),dtype = np.float32)
	for X in range(dim1):
		for Y in range(dim2):
			Q0XY[X,Y,:] = Q0(X,Y).flatten()
	print Q0XY
		
		
	sys.exit()


if __name__ == "__main__":
    main()
    
    
    
