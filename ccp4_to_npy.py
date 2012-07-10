#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
# ccp4 to npy
# Gael Goret for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr
###################################################################################

import numpy as np
import sys

headersize = 256
bytecode = 4 #float32


def main():
    fname = sys.argv[1]
    try:
        dim1 = int(sys.argv[2])
        dim2 = int(sys.argv[3])
        dim3 = int(sys.argv[4])
    except:
        raise Exception('Problem casting dims into integer')
    
    f = open(fname,'rb')
    f.seek(headersize*bytecode)
    raw = f.read()
    v = np.fromstring(raw,np.float32)
    v = v.reshape(dim1,dim2,dim3)
    v = v.transpose((1, 2, 0))
    np.save(fname[:-5],v)
    sys.exit()
    
#--------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print('Usage : python ccp4_to_npy.py 3D_volume.ccp4 dim1 dim2 dim3')
    else:
        main()
    
