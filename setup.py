#!/usr/bin/python

from numpy.distutils.core import Extension, setup
with open("README.rst") as fh:
    long_description = fh.read()
ext_internal = Extension(name = 'ace_internal', sources = ['ace.pyf',
                'ace.f', 'avas.f', 'rlsmo.f'
                ])
ext = Extension(name = 'ace_cream.ace_cream', sources = 'ace_cream.py')                
if __name__ == '__main__':
    setup(name = 'ace_cream',
          version = '0.4.post2',
          description = 'Alternating Conditional\
 Expectation Algorithm',
          author = 'zhaofeng-shu33',
          author_email = '616545598@qq.com',
          url = 'https://github.com/zhaofeng-shu33/ace_cream',
          maintainer = 'zhaofeng-shu33',
          maintainer_email = '616545598@qq.com',
          long_description = long_description,  
          install_requires = ['numpy'],
          ext_modules = [ext_internal],
          license = 'Apache License Version 2.0',
          packages = ['ace_cream'],
          classifiers = (
              "Development Status :: 4 - Beta",
              "Programming Language :: Fortran",
              "Programming Language :: Python :: 3.6",
              "Programming Language :: Python :: 3.5",
              "Programming Language :: Python :: 2.7",
              "Operating System :: OS Independent",
          ),
          )
