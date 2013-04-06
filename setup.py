"""
    Setup for SasView
    #TODO: Add checks to see that all the dependencies are on the system
"""
import sys
import os
import platform
import shutil
from setuptools import setup, Extension, find_packages
from distutils.command.build_ext import build_ext

try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
    NUMPY_INC = get_numpy_include_dirs()[0]
except:
    try:
        import numpy
        NUMPY_INC = os.path.join(os.path.split(numpy.__file__)[0], 
                                 "core","include")
    except:
        msg = "\nNumpy is needed to build SasView. "
        print msg, "Try easy_install numpy.\n  %s" % str(sys.exc_value)
        sys.exit(0)

# Manage version number ######################################
import sansview
VERSION = sansview.__version__
##############################################################

package_dir = {}
package_data = {}
packages = []
ext_modules = []

# Remove all files that should be updated by this setup
# We do this here because application updates these files from .sansview
# except when there is no such file
# Todo : make this list generic
plugin_model_list = ['polynominal5.py', 'sph_bessel_jn.py', 
                     'sum_Ap1_1_Ap2.py', 'sum_p1_p2.py', 
                     'testmodel_2.py', 'testmodel.py',
                     'polynominal5.pyc', 'sph_bessel_jn.pyc', 
                     'sum_Ap1_1_Ap2.pyc', 'sum_p1_p2.pyc', 
                     'testmodel_2.pyc', 'testmodel.pyc', 'plugins.log']
sans_dir = os.path.join(os.path.expanduser("~"),'.sasview')
if os.path.isdir(sans_dir):
    f_path = os.path.join(sans_dir, "sasview.log")
    if os.path.isfile(f_path):
        os.remove(f_path)
    f_path = os.path.join(sans_dir, "serialized_cat.p")
    if os.path.isfile(f_path):
        os.remove(f_path)
    f_path = os.path.join(sans_dir, 'config', "custom_config.py")
    if os.path.isfile(f_path):
        os.remove(f_path)
    f_path = os.path.join(sans_dir, 'plugin_models')
    if os.path.isdir(f_path):
        for file in os.listdir(f_path): 
            if file in plugin_model_list:
                file_path =  os.path.join(f_path, file)
                os.remove(file_path)
                    
# 'sys.maxsize' and 64bit: Not supported for python2.5
is_64bits = False
if sys.version_info >= (2, 6):
    is_64bits = sys.maxsize > 2**32
    
    
enable_openmp = True                    

if sys.platform =='darwin' and not is_64bits:
    # Disable OpenMP
    enable_openmp = False

# Options to enable OpenMP
copt =  {'msvc': ['/openmp'],
         'mingw32' : ['-fopenmp'],
         'unix' : ['-fopenmp']}
lopt =  {'msvc': ['/MANIFEST'],
         'mingw32' : ['-fopenmp'],
         'unix' : ['-lgomp']}

# Platform-specific link options
platform_lopt = {'msvc' : ['/MANIFEST']}

class build_ext_subclass( build_ext ):
    def build_extensions(self):
        # Get 64-bitness
        c = self.compiler.compiler_type
        print "Compiling with %s (64bit=%s)" % (c, str(is_64bits))
        
        # OpenMP build options
        if enable_openmp:
            if copt.has_key(c):
               for e in self.extensions:
                   e.extra_compile_args = copt[ c ]
            if lopt.has_key(c):
                for e in self.extensions:
                    e.extra_link_args = lopt[ c ]
                    
        # Platform-specific build options
        if platform_lopt.has_key(c):
            for e in self.extensions:
                e.extra_link_args = platform_lopt[ c ]

        build_ext.build_extensions(self)


# sans.invariant
package_dir["sans.invariant"] = "sansinvariant/src/sans/invariant"
packages.extend(["sans.invariant"])

# sans.guiframe
guiframe_path = os.path.join("sansguiframe", "src", "sans", "guiframe")
package_dir["sans.guiframe"] = guiframe_path
package_dir["sans.guiframe.local_perspectives"] = os.path.join(guiframe_path, 
                                                        "local_perspectives")
package_data["sans.guiframe"] = ['images/*', 'media/*']
packages.extend(["sans.guiframe", "sans.guiframe.local_perspectives"])
# build local plugin
for dir in os.listdir(os.path.join(guiframe_path, "local_perspectives")):
    if dir not in ['.svn','__init__.py', '__init__.pyc']:
        package_name = "sans.guiframe.local_perspectives." + dir
        packages.append(package_name)
        package_dir[package_name] = os.path.join(guiframe_path, 
                                                 "local_perspectives", dir)

# sans.dataloader
package_dir["sans.dataloader"] = os.path.join("sansdataloader", 
                                              "src", "sans", "dataloader")
package_data["sans.dataloader.readers"] = ['defaults.xml']
packages.extend(["sans.dataloader","sans.dataloader.readers"])

# sans.calculator
package_dir["sans.calculator"] = "sanscalculator/src/sans/calculator"
packages.extend(["sans.calculator"])

# sans.corfunc
package_dir["sans.corfunc"] = "corfunc/src/sans/corfunc"
packages.extend(["sans.corfunc"])
    
# sans.pr
numpy_incl_path = os.path.join(NUMPY_INC, "numpy")
srcdir  = os.path.join("pr_inversion", "src", "sans", "pr", "c_extensions")


    
package_dir["sans.pr.core"] = srcdir
package_dir["sans.pr"] = os.path.join("pr_inversion", "src","sans", "pr")
packages.extend(["sans.pr","sans.pr.core"])
ext_modules.append( Extension("sans.pr.core.pr_inversion",
                              sources = [ os.path.join(srcdir, "Cinvertor.c"),
                                         os.path.join(srcdir, "invertor.c"),
                                         ],
                              include_dirs=[numpy_incl_path],
                              ) )
        
# sans.fit (park integration)
package_dir["sans.fit"] = "park_integration/src/sans/fit"
packages.append("sans.fit")

# inversion view
package_dir["sans.perspectives"] = "inversionview/src/sans/perspectives"
package_dir["sans.perspectives.pr"] = "inversionview/src/sans/perspectives/pr"
packages.extend(["sans.perspectives","sans.perspectives.pr"])
package_data["sans.perspectives.pr"] = ['images/*']

# Invariant view
package_dir["sans.perspectives"] = os.path.join("invariantview", "src", 
                                                "sans", "perspectives")
package_dir["sans.perspectives.invariant"] = os.path.join("invariantview", 
                                    "src", "sans", "perspectives", "invariant")
                
package_data['sans.perspectives.invariant'] = [os.path.join("media",'*')]
packages.extend(["sans.perspectives","sans.perspectives.invariant"]) 

# Fitting view
fitting_path = os.path.join("fittingview", "src", "sans", 
                            "perspectives", "fitting")
package_dir["sans.perspectives"] = os.path.join("fittingview", 
                                            "src", "sans", "perspectives"),
package_dir["sans.perspectives.fitting"] = fitting_path
package_dir["sans.perspectives.fitting.plugin_models"] = \
                                os.path.join(fitting_path, "plugin_models")
package_data['sans.perspectives.fitting'] = ['media/*','plugin_models/*']
packages.extend(["sans.perspectives", "sans.perspectives.fitting", 
                 "sans.perspectives.fitting.plugin_models"])

# Calculator view
package_dir["sans.perspectives"] = "calculatorview/src/sans/perspectives"
package_dir["sans.perspectives.calculator"] = os.path.join("calculatorview", 
                                "src", "sans", "perspectives", "calculator")
package_data['sans.perspectives.calculator'] = ['images/*', 'media/*']
packages.extend(["sans.perspectives", "sans.perspectives.calculator"])
     
# Calculator view
package_dir["sans.perspectives"] = "corfuncview/src/sans/perspectives"
package_dir["sans.perspectives.corfunc"] = os.path.join("corfuncview", 
                                "src", "sans", "perspectives", "corfunc")
package_data['sans.perspectives.corfunc'] = ['images/*', 'media/*']
packages.extend(["sans.perspectives", "sans.perspectives.corfunc"])     
     
# Data util
package_dir["data_util"] = "sansutil"
packages.extend(["data_util"])

# Plottools
package_dir["danse"] = os.path.join("plottools", "src", "danse")
package_dir["danse.common"] = os.path.join("plottools", "src", 
                                           "danse", "common")
package_dir["danse.common.plottools"] = os.path.join("plottools", 
                                    "src", "danse", "common", "plottools")
packages.extend(["danse", "danse.common", "danse.common.plottools"])

# Park 1.2.1
package_dir["park"]="park-1.2.1/park"
packages.extend(["park"])
package_data["park"] = ['park-1.2.1/*.txt', 'park-1.2.1/park.epydoc']
ext_modules.append( Extension("park._modeling",
                              sources = [ os.path.join("park-1.2.1", 
                                                "park", "lib", "modeling.cc"),
                                         os.path.join("park-1.2.1", 
                                                "park", "lib", "resolution.c"),
                                         ],
                              ) )

# Sans models
includedir  = os.path.join("sansmodels", "include")
igordir = os.path.join("sansmodels", "src", "libigor")
c_model_dir = os.path.join("sansmodels", "src", "c_models")
smear_dir  = os.path.join("sansmodels", "src", "c_smearer")
gen_dir  = os.path.join("sansmodels", "src", "c_gen")
wrapper_dir  = os.path.join("sansmodels", "src", "python_wrapper", "generated")
model_dir = os.path.join("sansmodels", "src", "sans","models")

if os.path.isdir(wrapper_dir):
    for file in os.listdir(wrapper_dir): 
        file_path =  os.path.join(wrapper_dir, file)
        os.remove(file_path)
else:
    os.makedirs(wrapper_dir)
sys.path.append(os.path.join("sansmodels", "src", "python_wrapper"))
from wrapping import generate_wrappers
generate_wrappers(header_dir = includedir, 
                  output_dir = model_dir,
                  c_wrapper_dir = wrapper_dir)

IGNORED_FILES = [".svn"]
if not os.name=='nt':
    IGNORED_FILES.extend(["gamma_win.c","winFuncs.c"])


EXTENSIONS = [".c", ".cpp"]

def append_file(file_list, dir_path):
    """
    Add sources file to sources
    """
    for f in os.listdir(dir_path):
        if os.path.isfile(os.path.join(dir_path, f)):
            _, ext = os.path.splitext(f)
            if ext.lower() in EXTENSIONS and f not in IGNORED_FILES:
                file_list.append(os.path.join(dir_path, f)) 
        elif os.path.isdir(os.path.join(dir_path, f)) and \
                not f.startswith("."):
            sub_dir = os.path.join(dir_path, f)
            for new_f in os.listdir(sub_dir):
                if os.path.isfile(os.path.join(sub_dir, new_f)):
                    _, ext = os.path.splitext(new_f)
                    if ext.lower() in EXTENSIONS and\
                         new_f not in IGNORED_FILES:
                        file_list.append(os.path.join(sub_dir, new_f)) 
        
model_sources = []
append_file(file_list=model_sources, dir_path=igordir)
append_file(file_list=model_sources, dir_path=c_model_dir)
append_file(file_list=model_sources, dir_path=wrapper_dir)

smear_sources = []
append_file(file_list=smear_sources, dir_path=smear_dir)


package_dir["sans"] = os.path.join("sansmodels", "src", "sans")
package_dir["sans.models"] = model_dir

package_dir["sans.models.sans_extension"] = os.path.join("sansmodels", "src", 
                                            "sans", "models", "sans_extension")
            
package_data['sans.models'] = [os.path.join('media', "*.*")]
package_data['sans.models'] += [os.path.join('media','img', "*.*")]

packages.extend(["sans","sans.models","sans.models.sans_extension"])
    
smearer_sources = [os.path.join(smear_dir, "smearer.cpp"),
                  os.path.join(smear_dir, "smearer_module.cpp")]
geni_sources = [os.path.join(gen_dir, "sld2i_module.cpp")]
if os.name=='nt':
    smearer_sources.append(os.path.join(igordir, "winFuncs.c"))
    geni_sources.append(os.path.join(igordir, "winFuncs.c"))
ext_modules.extend( [ Extension("sans.models.sans_extension.c_models",
                                sources=model_sources,                 
                                include_dirs=[igordir, includedir, 
                                              c_model_dir, numpy_incl_path],
                                ),       
                    # Smearer extension
                    Extension("sans.models.sans_extension.smearer",
                              sources = smearer_sources,
                              include_dirs=[igordir, 
                                            smear_dir, numpy_incl_path],
                              ),
                    
                    Extension("sans.models.sans_extension.smearer2d_helper",
                              sources = [os.path.join(smear_dir, 
                                                "smearer2d_helper_module.cpp"),
                                         os.path.join(smear_dir, 
                                                "smearer2d_helper.cpp"),],
                              include_dirs=[smear_dir, numpy_incl_path],
                              ),
                    
                    Extension("sans.models.sans_extension.sld2i",
                              sources = [os.path.join(gen_dir, 
                                                "sld2i_module.cpp"),
                                         os.path.join(gen_dir, 
                                                "sld2i.cpp"),
                                         os.path.join(c_model_dir, 
                                                "libfunc.c"),
                                         os.path.join(c_model_dir, 
                                                "librefl.c"),],
                              include_dirs=[gen_dir, includedir, 
                                            c_model_dir, numpy_incl_path],
                              )
                    ] )
        
# SasView

package_dir["sans.sansview"] = "sansview"
package_data['sans.sansview'] = ['images/*', 'media/*', 'test/*', 
                                 'default_categories.p']
packages.append("sans.sansview")

#required = ['lxml>=2.2.2', 'numpy>=1.4.1', 'matplotlib>=0.99.1.1', 
#            'wxPython>=2.8.11', 'pil',
#            'periodictable>=1.3.0', 'scipy>=0.7.2']
required = ['lxml','periodictable>=1.3.1','pyparsing<2.0.0']

if os.name=='nt':
    required.extend(['html5lib', 'reportlab'])
else:
    required.extend(['pil'])
   
 # Set up SasView    
setup(
    name="sasview",
    version = VERSION,
    description = "SasView application",
    author = "University of Tennessee",
    author_email = "sansdanse@gmail.com",
    url = "http://danse.chem.utk.edu",
    license = "PSF",
    keywords = "small-angle x-ray and neutron scattering analysis",
    download_url = "https://sourceforge.net/projects/sansviewproject/files/",
    package_dir = package_dir,
    packages = packages,
    package_data = package_data,
    ext_modules = ext_modules,
    install_requires = required,
    zip_safe = False,
    entry_points = {
                    'console_scripts':[
                                       "sasview = sans.sansview.sansview:run",
                                       ]
                    },
    cmdclass = {'build_ext': build_ext_subclass }
    )   
