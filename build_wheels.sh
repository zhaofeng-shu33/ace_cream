#!/bin/bash
set -e -x

# Install a system package required by our library
# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install numpy
    "${PYBIN}/pip" wheel /io/ -w /io/build
done

# Bundle external shared libraries into the wheels
for whl in /io/build/*.whl; do
    auditwheel repair "$whl" -w /io/build/
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install ace_cream --no-index -f /io/build
    "${PYBIN}/python" /io/tests/ace_cream_test.py
done
