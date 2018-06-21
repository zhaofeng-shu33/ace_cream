Alternating Conditional Expectation Algorithm
=============================================

|Build Status|

This project provides a wrapper program of Python for ACE algorithm
implementation of Fortran.

How to build
------------

You need ``numpy`` and fortran compiler to build from source.

Windows
~~~~~~~

-  Install `Visual
   C++ <https://blogs.msdn.microsoft.com/vcblog/2017/03/07/msvc-the-best-choice-for-windows/>`__
   toolchain.

-  Download MinGW-w64 from
   `sourceforge <https://sourceforge.net/projects/mingw-w64/files/latest/download?source=typ_redirect>`__,
   which provides the necessary fortran compiler

-  Install MinGW-w64 and add ``{install_dir}\mingw64\bin`` path to
   environment variable (make ``gfortran`` accessible from command
   line).

-  (for conda environment) Add ``{install_dir}\Anaconda3\Scripts`` to
   environment variable (make ``f2py`` accessible from command line).

Mac
~~~

You can use package manager to install ``gfortran`` (included within gnu
compiler collection). For example, with ``Homebrew`` you can use

.. code:: shell

    brew install gcc

Ubuntu
~~~~~~

To install ``gfortran``, use the default package manager:

.. code:: shell

    sudo apt-get install gfortran

Run ``python setup.py install`` from command line at the project root
directory.

change log
----------

v0.1 initial commit v0.2 modify to relative import in ``__init__.py``
v0.3 add support for multiple columns of x and other directions of
transformation v0.4 add ``f_mapping`` function and unittests for this
function ## License

Apache License Version 2.0

.. |Build Status| image:: https://travis-ci.org/zhaofeng-shu33/ace_cream.svg?branch=master
   :target: https://travis-ci.org/zhaofeng-shu33/ace_cream
