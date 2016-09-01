# 3DRSR
3DRSR is a High-throughput low-resolution 2D/3D data reconstruction package. Called 3DRSR for 3D Reciprocal Space Reconstruction, the software is Capable of computing 3D low resolution reconstruction of reciprocal space, but also 2D high-resolution slices of the volume (images) and a lot more ...
The package is made of two main python modules and some tools/scripts written also in python :
An intelligent configuration file parser, which extracts from one XCALIBUR parameter file and from one in-house parameter file (format : see wiki). This program called 3DRSR_config generate a 3DRSR.conf ascii file required as input of 3DRSR.
The main program, 3DRSR, compute 2D/3D reconstuction, but also includes :

a multi-format reading (cf FabIO library), some corrections for polarization/flux density/parallax, 2D/3D symmetrization procedure, a naive interpolation, CCP4 output for intensity maps / mar2300-3450 for images.
ccp4_to_npy.py allow to convert the 3DRSR 3D map output format (CCP4) to NumPy binary file format.
npy_to_ccp4.py allow to convert .npy file (NumPy binary file) into the CCP4 3D map format (reading by such software as Chimera, Pymol, …)
Generator.py is tool allowing to generate a set of point group symmetry operations, by combining, a set of generator.
MDSP.py is a multidimensional (2D, 3D) symmetrization program, which compute a symmetrizing volume/image given as input a volume/image and a set of generators.

Licence: GPL
    
Wiki there -> https://forge.epn-campus.eu/projects/3drsr/wiki
