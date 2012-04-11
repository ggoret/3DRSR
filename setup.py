
import sys,os
import glob
import platform
import string
from distutils.core import Extension, setup
try:
    import numpy
except ImportError:
    text  = "You must have numpy installed.\n"
    text += "See http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=175103\n"
    raise ImportError, text


import distutils.sysconfig
global EDRSR_INSTALL_DIR
global EDRSR_SCRIPTS_DIR


jn = os.sep.join

# Append cvs tag if working from cvs tree
if os.path.isdir('.git') and os.path.isfile(os.sep.join(['.git', 'logs', 'HEAD'])):
    logs  = open( os.sep.join(['.git', 'logs', 'HEAD']) ).read()
    logs=string.replace(logs, '"""' , '"-"-"')
    open( jn( ["EDRSR", "logs.py"]) , "w").write('logs="""%s""" '%logs)

    os.system("git diff > diffoutput" )
    diffoutput= open("diffoutput").read()

    diffoutput=string.replace(diffoutput, '"""' , '"-"-"')
    open( jn(["EDRSR", "diffoutput.py"]), "w").write('diffoutput="""%s""" '%diffoutput)
    

packages = ['EDRSR', 'EDRSR.EDRSR_c' ]

if sys.platform == "win32":
    define_macros = [('WIN32',None)]
    script_files = [glob.glob( 'EDRSR/scripts/*')]
else:
    define_macros = []
    script_files = glob.glob( 'EDRSR/scripts/*')


files_to_check = glob.glob(jn(['EDRSR','EDRSR_c','*.i']))

for file in files_to_check:
    stat_i = os.stat( file ).st_mtime
    ok=0
    wrapfile = file[:-2]+"_wrap.c"
    if os.path.isfile( wrapfile):
        stat_cxx = os.stat( wrapfile  ).st_mtime

        print stat_cxx
        print stat_i
        
        if stat_cxx >   stat_i:
            ok=1
    if not ok :
        print " GENERATING SWIG WRAPPER FOR ", file
        
        if sys.platform == "win32":
            raise "  da farsi "
        else:
            os.system("swig  -python "+ file  )
                    

def touch_older(dependencies):
    for key in dependencies.keys():
        stat_target = os.stat( key ).st_mtime

        deps = dependencies[key]

        for dep in deps:
            stat_dep = os.stat(dep  ).st_mtime
            if stat_dep >   stat_target:
                if sys.platform == "win32":
                    os.system("copy "+ key+ " " + key+"__tmp"  )
                    os.system("copy "+ key+ "__tmp " + key  )
                    os.system("rm  "+ key+ "__tmp "  )
                else:
                    os.system("touch "+ key  )
            



def build_EDRSR_c_CPU(ext_modules,internal_clibraries):

    if sys.platform == "win32" :
        raise "fissa libs "
    else:  
        sdk_libs = ['gcc','m']

    c_sorgenti_tmp = ["fillvolume.c","fillvolume_wrap.c" ]

    c_sorgenti=[jn(['EDRSR','EDRSR_c',tok]) for tok in c_sorgenti_tmp]


    module  = Extension(name = 'EDRSR.EDRSR_c._fillvolume',
                        sources = c_sorgenti ,
                        depends=["EDRSR/EDRSR_c/fillvolume.h"], 
                        define_macros = define_macros,
                        include_dirs = [ "EDRSR/EDRSR_c/", numpy.get_include()])
    ext_modules.append(module)




ext_modules = []
internal_clibraries = []

build_EDRSR_c_CPU(ext_modules,internal_clibraries )


# data_files fix from http://wiki.python.org/moin/DistutilsInstallDataScattered
from distutils.command.install_data import install_data

from distutils.command.install_scripts import install_scripts
class smart_install_scripts(install_scripts):
    def run (self):
        global EDRSR_SCRIPTS_DIR, cudaversion, internal_clibraries
        
        #I prefer not to translate the python used during the build
        #process for the case of having an installation on a disk shared
        #by different machines and starting python from a shell script
        #that positions the environment
        from distutils import log
        from stat import ST_MODE
        
        install_cmd = self.get_finalized_command('install')
        #This is to ignore the --install-scripts keyword
        #I do not know if to leave it optional ...
        if False:
            self.install_dir = os.path.join(getattr(install_cmd, 'install_lib'), 'EDRSR')
            self.install_dir = os.path.join(self.install_dir, 'bin')        
        else:
            self.install_dir = getattr(install_cmd, 'install_scripts')
            self.install_lib_dir = getattr(install_cmd, 'install_lib')
            
        EDRSR_SCRIPTS_DIR = self.install_dir        
        if sys.platform != "win32":
            print "EDRSR scripts to be installed in %s" %  self.install_dir
        self.outfiles = self.copy_tree(self.build_dir, self.install_dir)
        self.outfiles = []
        for filein in glob.glob('EDRSR/scripts/*'):
            
            filedest = os.path.join(self.install_dir, os.path.basename(filein))
            filedest_py =  os.path.basename(filein)+".py"

            
            if os.path.exists(filedest):
                os.remove(filedest)
                
            moddir = os.path.join(getattr(install_cmd,'install_lib'), "EDRSR")

            f = open(filein, 'r')
            modfile = f.readline().replace("\n","")
            f.close()

            text  = "#!/bin/bash\n"
            text += "export PYTHONPATH=%s:%s/EDRSR_c/:${PYTHONPATH}\n" %( moddir,moddir)
            # text +="cudaversion=%s\n"%cudaversion
            # text +="nofbits=%s\n"%nofbits
            # text +="cudadir=`ls -td /usr/local/cuda-${cudaversion}* | head -1`\n"
            # text +="n=`ls -td ${cudadir}/lib64 | wc -l`\n"

#             text +="""

# if [  "${cudadir}" != "" ]; then
# if [  "$nofbits" == "64"   -a    "$n" == "1"   ] ; then
#     export LD_LIBRARY_PATH=${cudadir}/lib64:${LD_LIBRARY_PATH}
# else
#     export LD_LIBRARY_PATH=${cudadir}/lib:${LD_LIBRARY_PATH}
# fi
# fi
#             """

            text += "export LD_LIBRARY_PATH=%s:${LD_LIBRARY_PATH}\n" % moddir
            text += "export PATH=%s:${PATH}\n" % self.install_dir
            text += "exec python %s $*\n" %    os.path.join(moddir,filedest_py )
            f=open(filedest, 'w')
            f.write(text)
            f.close()
            #self.copy_file(filein, filedest)
            self.outfiles.append(filedest)



            

        for libfile in internal_clibraries:
            if os.path.isfile( libfile  ):
                os.system("rm %s/%s"% (moddir, libfile))
                assert(not os.path.isfile( "%s/%s"%(moddir,libfile)  ))
                os.system("cp %s %s" % ( libfile,    moddir  ) )
                assert(os.path.isfile( "%s/%s"%(moddir, libfile )  ))
                
                
        if os.name == 'posix':
            # Set the executable bits (owner, group, and world) on
            # all the scripts we just installed.
            for file in self.get_outputs():
                if self.dry_run:
                    log.info("changing mode of %s", file)
                else:
                    mode = ((os.stat(file)[ST_MODE]) | 0555) & 07777
                    log.info("changing mode of %s to %o", file, mode)
                    os.chmod(file, mode)

                    

class smart_install_mylibs(install_scripts):
    def run (self):
        print "CI SIAMO !!!!! "

   
description = ""
long_description = """
"""

distrib = setup(name="EDRSR",
                license = "GPL - Please read LICENSE.GPL for details",
                # version= logs,
                description = description,
                author = "Gael Goret, Alessandro Mirone",
                author_email="mirone@esrf.fr",
                url = "http://forge.epn-campus.fr/projects/EDRSR",
                long_description = long_description,
                packages = packages,
                platforms='any',
                ext_modules = ext_modules,
##                data_files = data_files,
##                package_data = package_data,
##                package_dir = {'': 'lib'},
    #            lib_files=["a"], 
                cmdclass = {
    'install_scripts':smart_install_scripts},
                scripts=script_files,
    

                )




















