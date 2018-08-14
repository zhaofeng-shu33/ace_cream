import os
from distutils.sysconfig import get_python_lib
# For gfortran+msvc combination, extra shared libraries may exist (stored by numpy.distutils)
if os.name == 'nt':
    HAS_DLL = False
    extra_dll_dir = os.path.join(get_python_lib(), 'ace_cream', '.libs')
    if os.path.isdir(extra_dll_dir):
        HAS_DLL = True
    else:
        extra_dll_dir = os.path.join(get_python_lib(), 'ace_internal', '.libs')
        if os.path.isdir(extra_dll_dir):
            HAS_DLL = True
    if(HAS_DLL):
        os.environ['PATH'] += os.pathsep + extra_dll_dir
    else:
        raise ImportError('directory %s/ace_cream/.libs or %s/ace_internal/.libs not found'%(get_python_lib(),get_python_lib()))
from .ace_cream import ace_cream
from .ace_cream import f_mapping