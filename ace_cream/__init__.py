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
        find_ace_internal_pyd = False
        for i in os.listdir(get_python_lib()):
            if(i.find('ace_internal')>=0 and i.find('pyd')>=0): 
                os.environ['PATH'] += get_python_lib()
                find_ace_internal_pyd = True
                break
        if not(find_ace_internal_pyd):
            raise ImportError('directory %s/ace_cream/.libs or %s/ace_internal/.libs not found'%(get_python_lib(),get_python_lib()))
from .ace_cream import ace_cream
from .ace_cream import f_mapping