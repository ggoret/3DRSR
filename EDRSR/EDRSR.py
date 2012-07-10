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


import numpy as np
import sys, time
__version__ = 0.1

import fillvolume_cy as fillvolume

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
    print "      \ \/___/    \/___/    \|__|        \/__/        \|__|   v %3.1f  \ " % __version__
    print "       \______________________________________________________________\ "
    print "                                                                         "
    print "      3DRSR : a software for volumetric Reconstruction of Reciprocal Space "
    print ""

#--------------------------------------------------------------------------------------------------------
# Rotation, Projection and Orientation Matrix
#--------------------------------------------------------------------------------------------------------


def Rotation(angle, rotation_axis=0):
    if type(rotation_axis) == type("") :
        rotation_axis = {"x":0, "y":1, "z":2}[rotation_axis]
    assert((rotation_axis >= 0 and rotation_axis <= 3))
    angle = np.radians(angle)
    ret_val = np.zeros([3, 3], np.float32)
    i1 = rotation_axis
    i2 = (rotation_axis + 1) % 3
    i3 = (rotation_axis + 2) % 3
    ret_val[i1, i1  ] = 1
    ret_val[i2, i2  ] = np.cos(angle)
    ret_val[i3, i3  ] = np.cos(angle)
    ret_val[i2, i3  ] = np.sin(angle)
    ret_val[i3, i2  ] = -np.sin(angle)
    return ret_val

#--------------------------------------------------------------------------------------------------------

def Arb_Rot(angle, rot_axis):
    angle = np.radians(angle)
    assert(len(rot_axis)==3)
    x,y,z = rot_axis
    rot_mat = np.array([[ 1 + (1-np.cos(angle))*(x*x-1) ,
                         -z*np.sin(angle)+(1-np.cos(angle))*x*y, 
                          y*np.sin(angle)+(1-np.cos(angle))*x*z ],
                          
                        [ z*np.sin(angle)+(1-np.cos(angle))*x*y ,
                          1 + (1-np.cos(angle))*(y*y-1),
                         -x*np.sin(angle)+(1-np.cos(angle))*y*z ],
                        
                        [-y*np.sin(angle)+(1-np.cos(angle))*x*z,
                          x*np.sin(angle)+(1-np.cos(angle))*y*z,
                          1 + (1-np.cos(angle))*(z*z-1) ]])
    return rot_mat

#--------------------------------------------------------------------------------------------------------

def Mirror(mirror_axis):
    x,y,z = {"x":[-1,1,1], "y":[1,-1,1], "z":[1,1,-1],0:[-1,1,1], 1:[1,-1,1], 2:[1,1,-1]}[mirror_axis]
    mat = np.array([[x,0,0],[0,y,0],[0,0,z]])
    return mat 

#--------------------------------------------------------------------------------------------------------

def ascii_gen2operator(generator):
    I = np.identity(3)
    gendict = {  '1':I,
            '1~':-I,
            '2(x)':Arb_Rot(180,[1,0,0]),
            '2(y)':Arb_Rot(180,[0,1,0]),
            '2(z)':Arb_Rot(180,[0,0,1]),
            '2(110)':np.array([[0,1,0],[1,0,0],[0,0,-1]]),
            '3(z)':Arb_Rot(120,[0,0,1]),
            '3(111)':np.array([[0,1,0],[0,0,1],[1,0,0]]),
            '4(z)':Arb_Rot(90,[0,0,1]),
            '4(z)~':-I*Arb_Rot(90,[0,0,1]),
            'm(x)':Mirror(0),
            'm(y)':Mirror(1),
            'm(z)':Mirror(2)
          }
          
    return gendict[generator]

#--------------------------------------------------------------------------------------------------------

def existing_combination(candidate_op,op_list,sigma):
    for op in op_list:
        absdiff = abs(candidate_op - op)
        if absdiff.sum()<sigma:
            return True
    return False
    
#--------------------------------------------------------------------------------------------------------

def genlist2oplist(genlist):
    base_op_list = [ascii_gen2operator(gen) for gen in genlist]
    I = np.identity(3)
    combined_op_list = [I]
    while 1:
        new_add_in_cycle = False
        for op in combined_op_list:
            for base_op in base_op_list:
                newop = np.dot(op,base_op)
                if not existing_combination(newop,combined_op_list,0.01):
                    combined_op_list.append(newop)
                    new_add_in_cycle = True
        if not new_add_in_cycle:
            break
    return np.array(combined_op_list).astype(np.float32)

#--------------------------------------------------------------------------------------------------------

def Prim_Rot_of_RS(omega, phi, kappa, alpha, beta, omega_offset):
    """ 
    Primary rotation of reciprocal space. 
    Omega, kappa, phi are the nominal values of the angle
    """
    tmp = Rotation(omega + omega_offset, 2)
    tmp = np.dot(tmp, Rotation(alpha, 1))
    tmp = np.dot(tmp, Rotation(kappa, 2))
    tmp = np.dot(tmp, Rotation(-alpha, 1))
    tmp = np.dot(tmp, Rotation(beta, 1))
    tmp = np.dot(tmp, Rotation(phi, 2))
    tmp = np.dot(tmp, Rotation(-beta, 1))
    return tmp

#--------------------------------------------------------------------------------------------------------

def DET(theta, theta_offset, d1, d2):
    """ 
    Rotation matrix for the detector 
    theta is the nominal theta value and d1,D2 are the tilts of detector
    """
    tmp = Rotation(theta + theta_offset, 2)
    tmp = np.dot(tmp, Rotation(d2, 1))
    tmp = np.dot(tmp, Rotation(d1, 0))
    return tmp

#--------------------------------------------------------------------------------------------------------    

def Snd_Rot_of_RS(r1, r2, r3):
    """ Secondary rotation of reciprocal space (to orient the crystallographic axis in a special way) """
    tmp = Rotation(r3, 2)
    tmp = np.dot(tmp, Rotation(r2, 1))
    tmp = np.dot(tmp, Rotation(r1, 0))
    return tmp

#--------------------------------------------------------------------------------------------------------

def P0(dist, b2):
    """ Primary projection of pixel coordinates (X,Y) to the reciprocal space. """
    B = Rotation(b2, 1) # Beam tilt matrix
    p0 = np.dot(B, [ dist, 0, 0 ])
    return p0

#--------------------------------------------------------------------------------------------------------
# Corrections
#--------------------------------------------------------------------------------------------------------

MD0 = np.array([[1, 0], [0, 1]  ], dtype=np.int32)
MD1 = np.array([[-1, 0], [0, 1] ], dtype=np.int32)
MD2 = np.array([[1, 0], [0, -1] ], dtype=np.int32)
MD3 = np.array([[-1, 0], [0, -1]], dtype=np.int32)
MD4 = np.array([[0, 1], [1, 0]  ], dtype=np.int32)
MD5 = np.array([[0, -1], [1, 0] ], dtype=np.int32)
MD6 = np.array([[0, 1], [-1, 0] ], dtype=np.int32)
MD7 = np.array([[0, -1], [-1, 0]], dtype=np.int32)

#--------------------------------------------------------------------------------------------------------
# Parameter Class
#--------------------------------------------------------------------------------------------------------

class Parameters(object):
    def __init__(self):
        #Feo : 
        """
        self.pixel_size = np.float32(0.172) # mm
        detector_distance = np.float32(174.42) # mm    
        self.dist = np.float32(detector_distance * 1.72)  # mm
        self.beam_tilt_angle = np.float32(0.99075) # deg

        self.det_origin_X = np.float32(1712.50841)
        self.det_origin_Y = np.float32(1733.19127)

        self.lmbda = np.float32(0.67018) # Wavelength users specified (ang)
        self.kappa = np.float32(-134.)
        self.alpha = np.float32(50.)
        self.beta = np.float32(0.)
        self.omega = np.float32(57.)
        self.theta = np.float32(0.)
        self.phi = np.float32(0.) #(n-1)*0.1 where n is the number of image

        self.omega_offset = np.float32(-0.19777)
        self.theta_offset = np.float32(0.39804)

        self.d1 = np.float32(-0.41144)
        self.d2 = np.float32(1.17097)

        self.r1 = np.float32(-88.788)
        self.r2 = np.float32(2.257)
        self.r3 = np.float32(69.629)

        self.pol_degree = np.float32(1.) # polarisation degree
        self.normal_to_pol = np.array([0, 0, 1], dtype=np.float32) # normal to polarisation plane
        """
        #GdFeSb :
        self.pixel_size = np.float32(0.172) # mm
        detector_distance = np.float32(174.42) # mm    
        self.dist = np.float32(detector_distance * 1.72)  # mm
        self.beam_tilt_angle = np.float32(1.06345) # deg

        self.det_origin_X = np.float32(1710.12344)
        self.det_origin_Y = np.float32(1734.16990)

        self.lmbda = np.float32(0.68870) # Wavelength users specified (ang)
        self.kappa = np.float32(-134.)
        self.alpha = np.float32(50.)
        self.beta = np.float32(0.)
        self.omega = np.float32(57.)
        self.theta = np.float32(0.)
        self.phi = np.float32(0.) #(n-1)*0.1 where n is the number of image

        self.omega_offset = np.float32(-0.31320)# -
        self.theta_offset = np.float32(0.33914)

        self.d1 = np.float32(-0.17552)
        self.d2 = np.float32(1.27807)

        self.r1 = np.float32(53.898)
        self.r2 = np.float32(29.653)
        self.r3 = np.float32(54.696)

        self.pol_degree = np.float32(1.) # polarisation degree
        self.normal_to_pol = np.array([0, 0, 1], dtype=np.float32) # normal to polarisation plane
        
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
    return max(np.sqrt(np.sum(l * l, axis= -1)))
    
#--------------------------------------------------------------------------------------------------------

class ProgressBar:
    def __init__(self, duration):
        self.duration = duration
        self.prog_bar = '[]'
        self.fill_char = '|'
        self.width = 60
        self.__update_amount(0)

    def update_time(self, elapsed_secs):
        self.__update_amount((elapsed_secs / float(self.duration)) * 100.0)
        self.prog_bar += ' %d/%d images' % (elapsed_secs, self.duration)
        
    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) / 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])
        
    def __str__(self):
        return str(self.prog_bar)
        
#--------------------------------------------------------------------------------------------------------
# CCP4 Writer
#--------------------------------------------------------------------------------------------------------

def write_ccp4_grid_data(volume, path):
    dtype = volume.dtype
    mtype = closest_ccp4_type(dtype)
    f = open(path, 'wb')

    header = ccp4_header(volume, mtype)
    f.write(header)

    I, J, K = volume.shape
    for k in range(K):
        matrix = volume[:, :, k]
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

    cell_size = map(lambda a, b: a * b, (1, 1, 1), size)

    if np.little_endian:
        machst = 0x00004144
    else:
        machst = 0x11110000

    ver_stamp = '3DSRS %s' % time.asctime()
    labels = [ver_stamp[:80]]

    nlabl = len(labels)
    # Make ten 80 character labels.
    labels.extend([''] * (10 - len(labels)))
    labels = [l + (80 - len(l)) * '\0' for l in labels]
    labelstr = ''.join(labels)

    dmin = volume.min()
    dmax = volume.max()
    dmean = volume.mean()

    strings = [
            binary_string(size, np.int32), # nx, ny, nz
            binary_string(mode, np.int32), # mode
            binary_string((0, 0, 0), np.int32), # nxstart, nystart, nzstart
            binary_string(size, np.int32), # mx, my, mz
            binary_string(cell_size, np.float32), # cella
            binary_string((90, 90, 90), np.float32), # cellb
            binary_string((1, 2, 3), np.int32), # mapc, mapr, maps
            binary_string((dmin, dmax, dmean), np.float32), # dmin, dmax, dmean
            binary_string(0, np.int32), # ispg
            binary_string(0, np.int32), # nsymbt
            binary_string([0] * 25, np.int32), # extra
            binary_string((0, 0, 0), np.float32), # origin
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
    if dtype in (np.float32, np.float64, np.int32, np.uint32, np.uint8, np.uint16):
        ctype = np.float32
    elif dtype in (np.int16, np.uint8):
        ctype = np.int16
    elif dtype in (np.int8, np.int0, np.character):
        ctype = np.int8
    else:
        raise TypeError, ('Volume data has unrecognized type %s' % dtype)

    return ctype
    
#--------------------------------------------------------------------------------------------------------

def padd_mar(data):
    dim1,dim2 = data.shape
    
    if dim1 <= 2300 and  dim2<= 2300:
        size = 2300
    else:
        size = 3450
            
    left = (size - dim1)//2
    right = size - (dim1 + left)
    up = (size - dim2)//2
    down = size - (dim2 + up) 
    out = np.zeros((size,size))
    
    if left > 0 : #pad
        outlm = left 
        inlm =  0
    else: #crop
        outlm = 0
        inlm = - left
    if right > 0: #pad
        outrm = -right
        inrm = dim1
    else: #crop
        outrm = size
        inrm = right
    if up > 0 : #pad
        outum = up 
        inum =  0
    else: #crop
        outum = 0
        inum = - up
    if down > 0: #pad
        outdm = - down
        indm = dim2
    else: #crop
        outdm = size
        indm = down
    
    out[outlm:outrm,outum:outdm] = data[inlm:inrm,inum:indm]
    return out,size

#--------------------------------------------------------------------------------------------------------
# Projection (whole volume)
#--------------------------------------------------------------------------------------------------------

def project_image(data, p0, Q0, XY_array_tmp, P_total_tmp, P_total_tmp_modulus, Qmax, params, Volume, Mask, C3, POL_tmp,Filter,sym_3_fold_axis,sym_3_fold_ops):

    dim1, dim2 = data.shape

    R = Prim_Rot_of_RS(params.omega, params.phi, params.kappa, params.alpha, params.beta, params.omega_offset)
    U = Snd_Rot_of_RS(params.r1, params.r2, params.r3)


    Q = np.tensordot (Q0 , R.T , axes=([2], [1]))
    Qfin = np.tensordot (Q , U.T , axes=([2], [1]))

    cube_dim = params.cube_dim

    dqx = dqy = dqz = Qmax

    q0x = q0y = q0z = 0
    fillvolume.func_somme(q0x, q0y, q0z, dqx, dqy, dqz, Volume, Mask, Qfin, data, POL_tmp, C3,Filter,sym_3_fold_axis,sym_3_fold_ops)

#--------------------------------------------------------------------------------------------------------
# Projection (slice)
#--------------------------------------------------------------------------------------------------------

def project_image_and_slice(data, p0, Q0, XY_array_tmp, P_total_tmp, P_total_tmp_modulus, Qmax, G, dQ0, dQ1, dQ2, Qoff, params, Image, Mask, C3, POL_tmp, Filter, apply_sym,ops_list):

    dim1, dim2 = data.shape

    R = Prim_Rot_of_RS(params.omega, params.phi, params.kappa, params.alpha, params.beta, params.omega_offset)
    U = Snd_Rot_of_RS(params.r1, params.r2, params.r3)

    Q = np.tensordot(Q0 , R.T , axes=([2],[1]))
    Qfin = np.tensordot(Q , U.T , axes=([2],[1]))
    
    #Qn = np.tensordot(Qfin, G, axes=([2],[1]))
    Qoff = np.dot(Qoff,G)
    
    fillvolume.slice_somme(dQ0,dQ1,dQ2,Qoff,Image,Mask,Qfin,data,POL_tmp,C3,Filter, apply_sym,ops_list,G)


#--------------------------------------------------------------------------------------------------------
# 3D Main
#--------------------------------------------------------------------------------------------------------

generators = ['1','1~','2(x)','2(y)','2(z)','2(110)','3(z)','3(111)','4(z)','4(z)~','m(x)','m(y)','m(z)']

def main3d(Filter_fname,flist):
    print(np.array(generators))
    print('---------------------------------------------------------------')
    counterror = 0
    while 1:
        if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
        print('Please enter a set of generators (a list of generator separated by comas) : ')
        try:
            ascii_gen_list = raw_input('generator_list = ').replace(',',' ').split()
        except:
            if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
            else:
                print('Entry for generator_list is not correct. Try Again you have still %d try '%(4-counterror))
                counterror+=1
        all_gen_in_list = True	        
        for gen in  ascii_gen_list:
            if gen not in generators:
                print ('Generator %s is not in the list. Try Again you have still %d try '%(gen,4-counterror))
                counterror+=1
                all_gen_in_list = False
        if all_gen_in_list:
            break
    print('------------------------------------------------------------')
    
    sym_3_fold_axis = 0
    sym_3_fold_ops = np.zeros((2,3,3)).astype(np.float32)
    ascii_3op = '3(z)'
    if ascii_3op in ascii_gen_list:
        id3 = ascii_gen_list.index(ascii_3op)
        del ascii_gen_list[id3]
        sym_3_fold_axis = 1
        sym_3_fold_ops = np.array([Arb_Rot(120,[0,0,1]),Arb_Rot(240,[0,0,1])]).astype(np.float32)
        print '------------------------------------------------------------'
        print '3 fold axis symmetry averaging will be computed direclty during the reconstruction to avoid the loss of resolution due to quadratic grid.'
        print '------------------------------------------------------------'
        
    ops_list = genlist2oplist(ascii_gen_list)

    print('------------------------------------------------------------')
    counterror = 0
    interp_factor = None
    while 1:
        try :
            print('\nPlease enter the interpolation factor. [default = 1 (no interpolation)] : ')
            itpfct = raw_input('interp_factor = ')
            if itpfct.strip()=='':
                interp_factor = 1
            else:
                interp_factor = int(itpfct)
        except:
            if(counterror>4):
                raise  Exception('Error', 'No more try ! bye bye !')
            else:
                print('\nEntry for interpolation factor is not correct. Try again you have still %d try '%(4-counterror))
                counterror+=1
        if interp_factor != None:
            if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
            break
    print('------------------------------------------------------------')
    
    time0 = time.time()
    print 'Reading File ...'

    #DEBUG OPTIONS
    rendering = False
    file_saving = True

    img = fabio.open(flist[0])
    time1 = time.time()
    print '-> t = %.2f s' % (time1 - time0)
    data = img.data
    dim1, dim2 = data.shape

    """
    dim1,dim2 = 2527,2463 
    data = np.ones((dim1,dim2),dtype = np.int32)*10000
    """
    print 'Image Dimension :', dim1, dim2

    print 'Setting Parameters ...'
    params = Parameters()

    p0 = P0(params.dist, params.beam_tilt_angle)
    MD = MD4

    time2 = time.time()
    print 'Computation of Initial Projection Coordinates Q0'
    Q0 = np.zeros((dim1, dim2, 3), dtype=np.float32)
    MD_pix_tmp = params.pixel_size * MD

    X_array_tmp = np.zeros((dim1, 2))
    X_array_tmp [:, 0] = np.arange(dim1) - (params.det_origin_X - 1725 + dim1 / 2.0)

    Y_array_tmp = np.zeros((dim2, 2))
    Y_array_tmp [:, 1] = np.arange(dim2) - (params.det_origin_Y - 1725 + dim2 / 2.0)

    X_array_tmp = np.tensordot(X_array_tmp , MD_pix_tmp , axes=([1], [1]))
    Y_array_tmp = np.tensordot(Y_array_tmp , MD_pix_tmp , axes=([1], [1]))

    XY_array_tmp = X_array_tmp[:, None, :] + Y_array_tmp[None, :, :]

    P_total_tmp = np.zeros((dim1, dim2 , 3), dtype=np.float32)
    P_total_tmp[:, :, 0 ] = -params.dist
    P_total_tmp[:, :, 1:3] = XY_array_tmp

    P_total_tmp = np.tensordot(P_total_tmp, DET(params.theta, params.theta_offset, params.d1, params.d2), axes=([2], [1]))

    P_total_tmp_modulus = np.sqrt(np.sum(P_total_tmp * P_total_tmp, axis= -1))

    Q0_tmp = P_total_tmp.T / P_total_tmp_modulus.T
    Q0 = ((Q0_tmp.T + p0 / params.dist) / params.lmbda).astype(np.float32)

    time3 = time.time()
    print '-> t = %.2f s' % (time3 - time2)

    time6 = time.time()
    if params.cpt_corr :
        print 'Computation of Polarisation Correction ...'
        
        P0xn = np.cross(p0, params.normal_to_pol, axis=0)
        P0xn_modulus = np.sqrt(np.sum(P0xn * P0xn, axis= -1))
        POL_tmp = params.pol_degree * (1 - ((P0xn * P_total_tmp).sum(axis= -1) / (P0xn_modulus * P_total_tmp_modulus)) ** 2)
        POL_tmp += (1 - params.pol_degree) * (1 - ((params.normal_to_pol * P_total_tmp).sum(axis=-1) / P_total_tmp_modulus) ** 2)
         
        #POL_tmp = ((params.dist ** 2 + (XY_array_tmp[:,:,1]*XY_array_tmp[:,:,1])) / (params.dist ** 2 + np.sum (XY_array_tmp * XY_array_tmp, axis= -1)))
        
    else :
        print 'Computation of Polarisation Correction : Canceled'


    POL_tmp = POL_tmp.astype(np.float32)


    time7 = time.time()
    print '-> t = %.2f s' % (time7 - time6)

    if params.cpt_corr:
        print 'Computation of Flux Density and Parallax Correction ...'
        C3 = (params.dist ** 3 / (params.dist ** 2 + np.sum (XY_array_tmp * XY_array_tmp, axis= -1)) ** (3 / 2)).astype(np.float32)
    else:
        print 'Computation of Flux Density and Parallax Correction : Canceled'

    time8 = time.time()
    print '-> t = %.2f s' % (time8 - time7)

    print 'Estimation of Qmax ...'
    # Estimation of Qmax
    corners = np.array([Q0[0, 0, :], Q0[0, dim2 - 1, :], Q0[dim1 - 1, 0, :], Q0[dim1 - 1, dim2 - 1, :]])
    Qmax = mag_max(corners) # maximal magnitude for reciprocal vector corresponding to the corners pixels
    print '-----------------------------'
    print 'Qmax = ', Qmax

    cube_dim = params.cube_dim
    Volume = np.zeros((cube_dim, cube_dim, cube_dim), dtype=np.float32)
    Mask = np.zeros((cube_dim, cube_dim, cube_dim), dtype=np.float32)
    total = len(flist)
    p = ProgressBar(total)
    nbfile = 0.
    print 'Loading images filter ...'
    f = open(Filter_fname,'rb')
    raw_data = f.read()
    Filter = np.fromstring(raw_data, np.uint8).reshape((dim1,dim2)).astype(np.float32)
    angular_step = 0.1
    #interp_factor = 1
    for fname in flist:
        timeI0 = time.time()
        img = fabio.open(fname)
        print 'Working on image %s' % fname
        data = img.data.astype(np.float32)
        timeI1 = time.time()
        print '-> time for opening image = %.2f s' % (timeI1 - timeI0)
        for j in range(interp_factor):
            params.phi = nbfile *  angular_step + (j/float(interp_factor))* angular_step
            project_image(data, p0, Q0, XY_array_tmp, P_total_tmp, P_total_tmp_modulus, Qmax, params, Volume, Mask, C3, POL_tmp,Filter,sym_3_fold_axis,sym_3_fold_ops)
            print('interpolation #%d on %d'%(j+1,interp_factor))
        nbfile += 1
        timeI2 = time.time()
        print '-> time for projecting image = %.2f s' % (timeI2 - timeI1)
        print '-> total time for this image = %.2f s' % (timeI2 - timeI0)
        print 'Global Progression : '
        
        p.update_time(nbfile)
        print '------------------------------------------------------------'
        print p
        print '------------------------------------------------------------'
        print '\n'
        

    time11 = time.time()
    print '3D Intensity Distribution : Done'
    threeDid = Volume[np.where(Volume > 0)]
    print threeDid
    vmin, vmax = min(threeDid), max(threeDid)
    print '-> Total time = %.2f s' % (time11 - time0)
##################################
    skip_avg = False
    if len(ops_list)==1:
        absdiff = abs(ops_list[0] - np.identity(3))
        if absdiff.sum()<0.001:#if op is only identity then no need to avg
            skip_avg = True
    if not skip_avg:
        print '------------------------------------------------------------'
        print 'NCS 3D intensity distribution averaging procedure START'
        print '------------------------------------------------------------'
        Cpt = Mask[:]
        Mask[np.where(Mask > 0)] = 1
        NCS_AVG = np.zeros_like(Volume)
        symopnb = 0
        for op in ops_list:
            print 'Application of symmetry operators #%s' % symopnb
            symopnb += 1
            print op
            tmp_intens,tmp_cpt = fillvolume.orient_volume(Volume, Mask, op)
            NCS_AVG += tmp_intens
            Cpt += tmp_cpt
        mapout = np.zeros_like(NCS_AVG)
        filter_ids = np.where(NCS_AVG!=0)
        mapout[filter_ids] = NCS_AVG[filter_ids]/Cpt[filter_ids]
    write_ccp4_grid_data( mapout, 'ncs_averaged_map.ccp4')
###################################
    print 'Normal END'

#--------------------------------------------------------------------------------------------------------
# 2D Main
#--------------------------------------------------------------------------------------------------------

def main2d(Filter_fname,flist):
    print('------------------------------------------------------------')
    counterror = 0
    vect_a = np.array([])
    vect_bp = np.array([])
    while 1:
        print('\nIn order to define the orientation of cut please enter vectors a and b (3 scalars separated by space for each entry) : ')
        try:
            vect_a = np.array(raw_input('a = ').replace(',',' ').split()).astype(np.float32)
            vect_bp = np.array(raw_input('b = ').replace(',',' ').split()).astype(np.float32)
        except:
            if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
            else:
                print('\nEntry for a and/or b is not correct. Try again you have still %d try '%(4-counterror))
                counterror+=1	        
        if vect_a.size != 3 or vect_bp.size != 3:
            if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
            print ('\nDimensions of a and/or b are not correct. Try again you have still %d try '%(4-counterror))
            counterror+=1
        else:
            break
    print('------------------------------------------------------------')
    counterror = 0
    dQ0 = None
    dQ1 = None
    dQ2 = None
    while 1:
        try :
            print('\nPlease enter the thickness of cut : ')
            dQ0 = float(raw_input('dQ0 = '))
            print('\nPlease enter the extent of sampling : ')
            dQ1 = float(raw_input('dQ1 = '))
            dQ2 = float(raw_input('dQ2 = '))
        except:
            if(counterror>4):
                raise  Exception('Error', 'No more try ! bye bye !')
            else:
                print('\nAt least one of the Entrys for dQ0, dQ1, dQ2 is not correct. Try again you have still %d try '%(4-counterror))
                counterror+=1
        if dQ0 != None and dQ1 != None  and dQ2 != None:
            if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
            break
    print('------------------------------------------------------------')
    counterror = 0
    Qoff = np.array([])
    while 1:
        print('\nPlease enter the offset vector : ')
        try:
            Qoff = np.array(raw_input('Qoff = ').replace(',',' ').split()).astype(np.float32)
        except:
            if(counterror>4):
                raise  Exception('Error', 'No more try ! bye bye !')
                print('\nEntry for Qoff is not correct. Try again you have still %d try '%(4-counterror))
                counterror+=1	        
        if Qoff.size != 3:
            if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
            print ('\nDimensions of a and/or b are not correct. Try again you have still %d try '%(4-counterror))
            counterror+=1
        else:
            break
    print('------------------------------------------------------------')
    counterror = 0
    interp_factor = None
    while 1:
        try :
            print('\nPlease enter the interpolation factor. [default = 1 (no interpolation)] : ')
            itpfct = raw_input('interp_factor = ')
            if itpfct.strip()=='':
                interp_factor = 1
            else:
                interp_factor = int(itpfct)
        except:
            if(counterror>4):
                raise  Exception('Error', 'No more try ! bye bye !')
            else:
                print('\nEntry for interpolation factor is not correct. Try again you have still %d try '%(4-counterror))
                counterror+=1
        if interp_factor != None:
            if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
            break
    print('------------------------------------------------------------')
    print vect_a
    print vect_bp
    vect_a = vect_a/np.linalg.norm(vect_a)
    vect_c = np.cross(vect_a,vect_bp)
    vect_c = vect_c/np.linalg.norm(vect_c)
    vect_b = np.cross(vect_c,vect_a)
    
    G = np.array([vect_a,vect_b,vect_c])
    
    print('------------------------------------------------------------')
    
    print G
    
    print('------------------------------------------------------------')
    
    print(np.array(generators))
    print('---------------------------------------------------------------')
    counterror = 0
    while 1:
        if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
        print('Please enter a set of generators (a list of generator separated by comas) : ')
        try:
            ascii_gen_list = raw_input('generator_list = ').replace(',',' ').split()
        except:
            if(counterror>4):
                raise  Exception('Error : No more try ! bye bye !')
            else:
                print('Entry for generator_list is not correct. Try Again you have still %d try '%(4-counterror))
                counterror+=1
        all_gen_in_list = True	        
        for gen in  ascii_gen_list:
            if gen not in generators:
                print ('Generator %s is not in the list. Try Again you have still %d try '%(gen,4-counterror))
                counterror+=1
                all_gen_in_list = False
        if all_gen_in_list:
            break
    print('------------------------------------------------------------')
    
    ops_list = genlist2oplist(ascii_gen_list)
    apply_sym=0
    if len(ops_list)>1:
        apply_sym=1
    
        
    
    time0 = time.time()
    print 'Reading File ...'

    #DEBUG OPTIONS
    rendering = False
    file_saving = True

    img = fabio.open(flist[0])
    time1 = time.time()
    print '-> t = %.2f s' % (time1 - time0)
    data = img.data
    dim1, dim2 = data.shape

    print 'Image Dimension :', dim1, dim2

    print 'Setting Parameters ...'
    params = Parameters()

    p0 = P0(params.dist, params.beam_tilt_angle)
    MD = MD4

    time2 = time.time()
    print 'Computation of Initial Projection Coordinates Q0'
    Q0 = np.zeros((dim1, dim2, 3), dtype=np.float32)
    MD_pix_tmp = params.pixel_size * MD

    X_array_tmp = np.zeros((dim1, 2))
    X_array_tmp [:, 0] = np.arange(dim1) - (params.det_origin_X - 1725 + dim1 / 2.0)

    Y_array_tmp = np.zeros((dim2, 2))
    Y_array_tmp [:, 1] = np.arange(dim2) - (params.det_origin_Y - 1725 + dim2 / 2.0)

    X_array_tmp = np.tensordot(X_array_tmp , MD_pix_tmp , axes=([1], [1]))
    Y_array_tmp = np.tensordot(Y_array_tmp , MD_pix_tmp , axes=([1], [1]))

    XY_array_tmp = X_array_tmp[:, None, :] + Y_array_tmp[None, :, :]

    P_total_tmp = np.zeros((dim1, dim2 , 3), dtype=np.float32)
    P_total_tmp[:, :, 0 ] = -params.dist
    P_total_tmp[:, :, 1:3] = XY_array_tmp

    P_total_tmp = np.tensordot(P_total_tmp, DET(params.theta, params.theta_offset, params.d1, params.d2), axes=([2], [1]))

    P_total_tmp_modulus = np.sqrt(np.sum(P_total_tmp * P_total_tmp, axis= -1))

    Q0_tmp = P_total_tmp.T / P_total_tmp_modulus.T
    Q0 = ((Q0_tmp.T + p0 / params.dist) / params.lmbda).astype(np.float32)

    time3 = time.time()
    print '-> t = %.2f s' % (time3 - time2)

    time6 = time.time()
    if params.cpt_corr :
        print 'Computation of Polarisation Correction ...'
        
        P0xn = np.cross(p0, params.normal_to_pol, axis=0)
        P0xn_modulus = np.sqrt(np.sum(P0xn * P0xn, axis= -1))
        POL_tmp = params.pol_degree * (1 - ((P0xn * P_total_tmp).sum(axis= -1) / (P0xn_modulus * P_total_tmp_modulus)) ** 2)
        POL_tmp += (1 - params.pol_degree) * (1 - ((params.normal_to_pol * P_total_tmp).sum(axis=-1) / P_total_tmp_modulus) ** 2)        
    else :
        print 'Computation of Polarisation Correction : Canceled'

    POL_tmp = POL_tmp.astype(np.float32)

    time7 = time.time()
    print '-> t = %.2f s' % (time7 - time6)

    if params.cpt_corr:
        print 'Computation of Flux Density and Parallax Correction ...'
        C3 = (params.dist ** 3 / (params.dist ** 2 + np.sum (XY_array_tmp * XY_array_tmp, axis= -1)) ** (3 / 2)).astype(np.float32)
    else:
        print 'Computation of Flux Density and Parallax Correction : Canceled'

    time8 = time.time()
    print '-> t = %.2f s' % (time8 - time7)

    print 'Estimation of Qmax ...'
    # Estimation of Qmax
    corners = np.array([Q0[0, 0, :], Q0[0, dim2 - 1, :], Q0[dim1 - 1, 0, :], Q0[dim1 - 1, dim2 - 1, :]])
    Qmax = mag_max(corners) # maximal magnitude for reciprocal vector corresponding to the corners pixels
    print '-----------------------------'
    print 'Qmax = ', Qmax
    size = 2300
    Image = np.zeros((size,size), dtype=np.float32)#dimension can change !
    Mask = np.zeros((size,size), dtype=np.float32)
    total = len(flist)
    p = ProgressBar(total)
    nbfile = 0.
    print 'Loading images filter ...'
    f = open(Filter_fname,'rb')
    raw_data = f.read()
    Filter = np.fromstring(raw_data, np.uint8).reshape((dim1,dim2)).astype(np.float32)
    angular_step = 0.1
    #interp_factor = 1
    for fname in flist:

        timeI0 = time.time()

        img = fabio.open(fname)
        print 'Working on image %s' % fname
        data = img.data.astype(np.float32)
        
        timeI1 = time.time()
        print '-> time for opening image = %.2f s' % (timeI1 - timeI0)
        
        for j in range(interp_factor):
            params.phi = nbfile *  angular_step + (j/float(interp_factor))* angular_step
            project_image_and_slice(data, p0, Q0, XY_array_tmp, P_total_tmp, P_total_tmp_modulus, Qmax, G, dQ0, dQ1, dQ2, Qoff, params, Image, Mask, C3, POL_tmp, Filter, apply_sym,ops_list)
            print('interpolation #%d on %d'%(j+1,interp_factor))
        nbfile += 1
        timeI2 = time.time()
        print '-> time for projecting image = %.2f s' % (timeI2 - timeI1)
        print '-> total time for this image = %.2f s' % (timeI2 - timeI0)
        print 'Global Progression : '
        
        p.update_time(nbfile)
        print '------------------------------------------------------------'
        print p
        print '------------------------------------------------------------'
        print '\n'
        

    time11 = time.time()
    print '2D Intensity Distribution : Done'
    print '-> Total time = %.2f s' % (time11 - time0)
    mapout = np.zeros_like(Image)
    mapout[np.where(Image!=0)] = Image[np.where(Image!=0)]/Mask[np.where(Image!=0)]
    mardata,size = padd_mar(mapout)
    out = fabio.mar345image.mar345image(data=mardata, header=img.header)
    out.write('map2D.mar%d'%size)
    np.save('map2D.npy', mapout)
    np.save('mask2D.npy', Mask)
    print 'Normal END'

#--------------------------------------------------------------------------------------------------------
# Main
#--------------------------------------------------------------------------------------------------------

def main():
    display_logo()
    Filter_fname = sys.argv[1]
    flist = sys.argv[2:]
    counterror = 0
    while 1:
	    choice = raw_input('Would you like to reconstruct the whole volume [v] or a slice [s] : ')
	    if choice.strip()=='':
	        choice = 'v'
	        print('Default choice : volume reconstruction')
	        break
	    elif choice not in ['V','v','s','S']:
		    if(counterror<4):
			    print ('key %s are not a correct. Try Again you have still %d try '%(choice,4-counterror))
		    else:
			    raise  Exception('Error','key %s are not correct. No more try'%choice)
		    counterror+=1
	    else:
		    break
    if choice in ['V','v']:
        main3d(Filter_fname,flist)
    elif choice in ['S','s']:
        main2d(Filter_fname,flist)



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Usage : python 3DRSR.py image-file(s)')
    else:
        main()
    sys.exit()


