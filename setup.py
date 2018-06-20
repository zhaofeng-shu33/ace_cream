#!/usr/bin/python

from numpy.distutils.core import Extension, setup

ext_internal = Extension(name = 'ace_internal', sources = ['ace.pyf',
                'ace.f', 'avas.f', 'rlsmo.f'
                ])
ext = Extension(name = 'ace_cream.ace_cream', sources = 'ace_cream.py')                
if __name__ == '__main__':
    setup(name = 'ace_cream',
          version = '0.2',
          description = 'alternating conditional\
 expectation algorithm',
          author = 'zhaofeng-shu33',
          author_email = '616545598@qq.com',
          maintainer = 'zhaofeng-shu33',
          maintainer_email = '616545598@qq.com',
          long_description = 'wrapper of Fortran implementation of\
 alternating conditional expectation algorithm',          
          ext_modules = [ext_internal],
          license = 'Apache License Version 2.0',
          packages = ['ace_cream']
          )
