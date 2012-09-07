#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
# parameters extractor from 3D Reciprocal Space Reconstruction programm
# Gael Goret for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr
###################################################################################
import sys, time
import numpy as np
__version__ = 0.1

class Parameters(object):
    def __init__(self):
        #example :
        self.pixel_size = None # mm
        self.dist =  None # np.float32(detector_distance * 1.72)  # mm
        self.beam_tilt_angle = None # deg

        self.det_origin_X = None
        self.det_origin_Y = None

        self.lmbda = None # Wavelength users specified (ang)
        self.kappa = None
        self.alpha = None
        self.beta = None
        self.omega = None
        self.theta = None
        self.phi = None #(n-1)*0.1 where n is the number of image

        self.omega_offset = None #-
        self.theta_offset = None

        self.d1 = None
        self.d2 = None

        self.r1 = None
        self.r2 = None
        self.r3 = None

        self.pol_degree = np.float32(1.) # polarisation degree
        self.normal_to_pol = [0, 0, 1] # normal to polarisation plane
        
        self.cube_dim = None
        
        self.angular_step = None
        #DEBUG OPTIONS :
        self.cpt_corr = True

def main():
    parfname = sys.argv[1]
    cfgfname = sys.argv[2]
    parf = open(parfname,'r')
    cfgf = open(cfgfname,'r')
    param = Parameters()
    for l in parf:
        if ('ALPHA (DEG)' in l) and ('BETA (DEG)' in l):
            param.alpha, param.beta = float(l.split()[4]),np.float32(l.split()[7])
        if 'WAVELENGTH USERSPECIFIED (ANG)' in l:
            param.lmbda = np.float32(l.split()[6])
        if 'DETECTOR ZERO (PIX, 1X1 BINNING)' in l:
            param.det_origin_X, param.det_origin_Y = np.float32(l.split()[8]),np.float32(l.split()[10])
        if 'X-RAY BEAM ORIENTATION (DEG)' in l:
            param.beam_tilt_angle = np.float32(l.split()[7])
        if 'DETECTOR ROTATION (DEG)' in l:
            param.d1,param.d2 = np.float32(l.split()[6]),np.float32(l.split()[8])
        if 'DETECTOR DISTANCE (MM)' in l:
            param.dist = np.float32(l.split()[5])*1.72 # constant ?
        if 'SOFTWARE ZEROCORRECTION (DEG)' in l:
            param.omega_offset, param.theta_offset = np.float32(l.split()[6]),np.float32(l.split()[8])
    for l in cfgf:
         if '1st tilt of U matrix' in l:
            param.r1 = np.float32(l[1:l.find('-')])
         if '2nd tilt of U matrix' in l:
            param.r2 = np.float32(l[1:l.find('-')])
         if '3rd tilt of U matrix' in l:
            param.r3 = np.float32(l[1:l.find('-')])
         if 'oscillation size in degrees' in l :
            param.angular_step = np.float32(l[1:l.find('-')])

    parf.close()
    cfgf.close()
    print('\n')
    print('------------------------------------------------------------')
    print('                  3DRSR Start-up Procedure v%3.1f '%__version__)
    print('------------------------------------------------------------\n')
    print('Initial 3DRSR Parameter Set (auto-assignment) :')
    for attr in dir(param):
        if '__' not in attr:
            exec("tmp = param.%s"%attr)
            print "param.%s = "%attr, tmp
    print('------------------------------------------------------------')  
    for attr in dir(param):
        if '__' not in attr:
            exec("tmp = param.%s"%attr)
            if tmp == None:
                counterror = 0
                while 1:
                    print('\n%s has not been defined by auto-assignment.\nSet a value for the parameter : '%attr)
                    try:
                        tmp = np.float32(raw_input('%s = '%attr))
                    except:
                        if(counterror>4):
                            raise  Exception('Error : No more try ! bye bye !')
                        else:
                            print('\nEntry for %s is not correct. Try again you have still %d try '%(attr,4-counterror))
                            counterror+=1	        
                    if tmp != None:
                        break
                print('------------------------------------------------------------')
                exec("param.%s = tmp"%attr)
    print('Final 3DRSR Parameter Set :')
    for attr in dir(param):
        if '__' not in attr:
            exec("tmp = param.%s"%attr)
            print "param.%s = "%attr, tmp
    print('------------------------------------------------------------')
    
    confout = open('3DRSR.conf','w')
    
    confout.write('#\n#3DRSR Configuration file\n#\n\n')
    for attr in dir(param):
        if '__' not in attr:
            exec("tmp = param.%s"%attr)
            confout.write("%s = %s\n"%(attr, tmp))
    confout.close()
    print('Configuration file 3DRSR.conf written and ready to be used !')
    

     
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print('Usage : python extract_param.py parameter_file.par Crysalis_config_file.cfg')
    else:
        main()
    sys.exit()

