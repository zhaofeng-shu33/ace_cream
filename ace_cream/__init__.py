import os
from distutils.sysconfig import get_python_lib
# For gfortran+msvc combination, extra shared libraries may exist (stored by numpy.distutils)
if os.name == 'nt':
    extra_dll_dir = os.path.join(get_python_lib(), 'ace_cream', '.libs')
    if os.path.isdir(extra_dll_dir):
        os.environ['PATH'] += os.pathsep + extra_dll_dir
    else:
        raise ImportError('ace_cream/.libs not found')
from .ace_cream import ace_cream
from .ace_cream import f_mapping