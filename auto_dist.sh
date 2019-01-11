#! /bin/bash
# automatic uploads to pypi.org after version changes

rm -r dist
python setup.py sdist bdist
# make .rst from .md
pandoc --from=markdown --to=rst --output=README.rst README.md
twine upload -u zhaofeng-shu33 dist
