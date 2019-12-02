#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
# MDSP : Multi-Dimentional Symetrization Programm 
# Gael Goret for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr
###################################################################################
try :
    import fabio
except :
    print 'Warning : fabio module could not be initialised'
import numpy as np
import sys, time
import rotsum
__version__ = 0.1


#--------------------------------------------------------------------------------------------------------
# Logo
#--------------------------------------------------------------------------------------------------------

def display_logo():
    print "\n-----------------------------------------------------------------------"
    print " _______  ______   ______  ______  |  ______  ______   ______  _______ "
    print "(_______)(______) / _____)(_____ \ | / _____)(_____ \ (______)(_______)"
    print " _  _  _  _     _( (____   _____) )|( (_____   ____) )_     _  _  _  _ "
    print "| ||_|| || |   | |\____ \ |  ____/ | \____  | / ____/| |   | || ||_|| |"
    print "| |   | || |__/ / _____) )| |      |      | |( (_____ \ \__| || |   | |"
    print "|_|   |_||_____/ (______/ |_|      |      |_| \______) \_____||_|   |_|"
    print ""
    print " MDSP : Multi-Dimensional Symmetrization Program : v %3.1f " % __version__
    print "-----------------------------------------------------------------------\n"


#--------------------------------------------------------------------------------------------------------
# Rotation, Projection and Orientation Matrix
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

def Arb_Rot_2d(angle):
    angle = np.radians(angle)
    rot_mat = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    return rot_mat

#--------------------------------------------------------------------------------------------------------

def Mirror(mirror_axis):
    x,y,z = {"x":[-1,1,1], "y":[1,-1,1], "z":[1,1,-1],0:[-1,1,1], 1:[1,-1,1], 2:[1,1,-1]}[mirror_axis]
    mat = np.array([[x,0,0],[0,y,0],[0,0,z]])
    return mat 

#--------------------------------------------------------------------------------------------------------

generators3d = ['1~','2(x)','2(y)','2(z)','2(110)','3(z)','3(111)','4(z)','4(z)~','m(x)','m(y)','m(z)']

generators2d = ['-I','2','3','4','6','m(h)','m(v)']

def ascii_gen2operator(generator):
    gendict = {  '1':np.identity(3),
            '1~':-np.identity(3),
            '-I':-np.identity(2),
            '2':Arb_Rot_2d(180),
            '2(x)':Arb_Rot(180,[1,0,0]),
            '2(y)':Arb_Rot(180,[0,1,0]),
            '2(z)':Arb_Rot(180,[0,0,1]),
            '2(110)':np.array([[0,1,0],[1,0,0],[0,0,-1]]),
            '3':Arb_Rot_2d(120),
            '3(z)':Arb_Rot(120,[0,0,1]),
            '3(111)':np.array([[0,1,0],[0,0,1],[1,0,0]]),
            '4':Arb_Rot_2d(90),
            '4(z)':Arb_Rot(90,[0,0,1]),
            '4(z)~':-np.identity(3)*Arb_Rot(90,[0,0,1]),
            '6':Arb_Rot_2d(60),
            'm(h)':np.array([[-1,0],[0,1]]),
            'm(v)':np.array([[1,0],[0,-1]]),
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

def genlist2oplist(genlist,dim):
    base_op_list = [ascii_gen2operator(gen) for gen in genlist]
    I = np.identity(dim)
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
# CCP4 Reader
#--------------------------------------------------------------------------------------------------------

def read_ccp4(fname):
    f = open(fname,'rb')
    f.seek(0,2)                              # go to end of file
    file_size = f.tell()
    f.seek(0,0)                              # beginning
    esize = np.array((), np.int32).itemsize
    count = 3
    string = f.read(esize * count)
    dim1, dim2, dim3 = np.fromstring(string, np.int32)
    
    headersize = 256
    bytecode = 4
    f.seek(headersize*bytecode)
    raw = f.read()
    v = np.fromstring(raw,np.float32)
    v = v.reshape(dim1,dim2,dim3)
    v = v.transpose((1, 2, 0))
    return v
    
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
# Main 2D
#--------------------------------------------------------------------------------------------------------

def main2d(fname,outfname,debug):

    print "-----------------------------------------------------------------------"
    print(np.array(generators2d))
    print "-----------------------------------------------------------------------"
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
            if gen not in generators2d:
                print ('Generator %s is not in the list. Try Again you have still %d try '%(gen,4-counterror))
                counterror+=1
                all_gen_in_list = False
        if all_gen_in_list:
            break
    print "-----------------------------------------------------------------------\n"
    
    f = fabio.open(fname)
    img = f.data.astype(np.float32)
    ops_list = genlist2oplist(ascii_gen_list,2)
    print ops_list
    
    skip_avg = False
    if len(ops_list)==1:
        absdiff = abs(ops_list[0] - np.identity(2))
        if absdiff.sum()<0.001:#if op is only identity then no need to avg
            skip_avg = True
    mapout = img
    if not skip_avg:
        print '------------------------------------------------------------'
        print 'NCS 3D intensity distribution averaging procedure START'
        print '------------------------------------------------------------'
        symopnb = 0
        Mask  = np.ones_like(img) 
        ncs_avg_intens = np.zeros_like(img)
        cpt = np.zeros_like(img)
        for op in ops_list:
            print 'Application of symmetry operators #%s' % symopnb
            symopnb += 1
            print op
            tmp_intens,tmp_cpt = rotsum.orient_image(img, Mask, op)
            ncs_avg_intens += tmp_intens
            cpt += tmp_cpt
        mapout = np.zeros_like(img)
        filter_ids = np.where(ncs_avg_intens!=0)
        mapout[filter_ids] = ncs_avg_intens[filter_ids]/cpt[filter_ids]
    out = fabio.mar345image.mar345image(data=mapout.astype(np.int32), header=f.header)
    out.write(outfname + '.mar')
    np.save(outfname, mapout)
    
    
#--------------------------------------------------------------------------------------------------------
# Main 3D
#--------------------------------------------------------------------------------------------------------


def main3d(fname,outfname,debug):
    print "-----------------------------------------------------------------------"
    print(np.array(generators3d))
    print "-----------------------------------------------------------------------"
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
            if gen not in generators3d:
                print ('Generator %s is not in the list. Try Again you have still %d try '%(gen,4-counterror))
                counterror+=1
                all_gen_in_list = False
        if all_gen_in_list:
            break
    print "-----------------------------------------------------------------------\n"
    
    Volume = read_ccp4(fname)   
    ops_list = genlist2oplist(ascii_gen_list,3)
    skip_avg = False
    if len(ops_list)==1:
        absdiff = abs(ops_list[0] - np.identity(3))
        if absdiff.sum()<0.001:#if op is only identity then no need to avg
            skip_avg = True
    mapout = Volume
    if not skip_avg:
        print '------------------------------------------------------------'
        print 'NCS 3D intensity distribution averaging procedure START'
        print '------------------------------------------------------------'
        symopnb = 0
        Mask  = np.ones_like(Volume) 
        ncs_avg_intens = np.zeros_like(Volume)
        cpt = np.zeros_like(Volume)
        T = 0
        for op in ops_list:
            print 'Application of symmetry operators #%s' % symopnb
            symopnb += 1
            print op
            t = time.time()
            if debug:
                tmp_intens,tmp_cpt = rotsum.orient_volume_opt(Volume, Mask, op)
            else:
                tmp_intens,tmp_cpt = rotsum.orient_volume(Volume, Mask, op)
            T += time.time() - t
            ncs_avg_intens += tmp_intens
            cpt += tmp_cpt
        print 'Total Time = ',T
        mapout = np.zeros_like(Volume)
        filter_ids = np.where(ncs_avg_intens!=0)
        mapout[filter_ids] = ncs_avg_intens[filter_ids]/cpt[filter_ids]
    write_ccp4_grid_data( mapout, outfname + '.ccp4')
    
#--------------------------------------------------------------------------------------------------------
# Main
#--------------------------------------------------------------------------------------------------------

def main():
    display_logo()
    fname = sys.argv[1]
    outfname = sys.argv[2]
    try:
        debug = bool(int(sys.argv[3]))
    except:
        debug = False
    
    counterror = 0
    while 1:
	    choice = raw_input('What kind of data you would like to symmetrize : \n3D volume ? -> [v] \nor \n2D image ? -> [i] \n')
	    if choice.strip()=='':
	        choice = 'v'
	        print('Default choice : 3D volume')
	        break
	    elif choice not in ['V','v','i','I']:
		    if(counterror<4):
			    print ('key %s are not a correct. Try Again you have still %d try '%(choice,4-counterror))
		    else:
			    raise  Exception('Error','key %s are not correct. No more try'%choice)
		    counterror+=1
	    else:
		    break
    if choice in ['V','v']:
        main3d(fname,outfname,debug)
    elif choice in ['I','i']:
        main2d(fname,outfname,debug)

#--------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Usage : python MDSP.py data_to_symmetrize output_file_name')
    else:
        main()
    sys.exit()

