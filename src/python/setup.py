#Install via pip with 'pip install -install -e .'
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
import os

sourcefiles = ["cython_wrapper.pyx"]

cpp_dir = os.getenv('cpp_dir')

ext_modules = [Extension("max_mutual_info", sourcefiles,
                library_dirs = [cpp_dir],
                define_macros = [('NOSHORTS',None)],
                include_dirs = [numpy.get_include(), cpp_dir],
                libraries = ['cpp_version'],
                language="c++",
                extra_compile_args=["-O3"])]

for e in ext_modules:
    e.cython_directives = {'language_level': "3"} #all are Python-3

setup(
    name = 'Sensor Placement via Mutual Information SDP',
    version = '0.0.1',
    author = ['Robert Bassett', 'Erik Vargas', 'Jefferson Huang'],
    author_email = ['robert.bassett@nps.edu', 'erik.vargas@nps.edu', 'jefferson.huang@nps.edu'],
    url = 'TBD',
    license = 'CC0',
    description = 'Fast rounder for Mutual information SDP.',
    install_requires=['numpy', 'cython', 'python>=3.0']
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    zip_safe=False
    )
            
