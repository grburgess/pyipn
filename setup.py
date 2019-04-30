#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""


from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
from setuptools import setup, find_packages, Command, Extension
import os
import io
import sys
from shutil import rmtree



class UploadCommand(Command):
    """Support setup.py upload."""

    description = "Build and publish the package."
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print("\033[1m{0}\033[0m".format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status("Removing previous builds...")
            rmtree(os.path.join(here, "dist"))
        except OSError:
            pass

        self.status("Building Source and Wheel (universal) distribution...")
        os.system("{0} setup.py sdist bdist_wheel --universal".format(sys.executable))

        self.status("Uploading the package to PyPI via Twine...")
        os.system("twine upload dist/*")

        self.status("Pushing git tags...")
        os.system("git tag v{0}".format(about["__version__"]))
        os.system("git push --tags")

        sys.exit()



# Package meta-data.
NAME = "pyipn"
DESCRIPTION = "A generic IPN simulator"
URL = 'https://github.com/grburgess/pyipn'
EMAIL = 'jburgess@mpe.mpg.de'
AUTHOR = "J. Michael Burgess"
REQUIRES_PYTHON = ">=2.7.0"
VERSION = None

REQUIRED = [
    "numpy",
    "scipy",
    "ipython",
    "matplotlib",
    "h5py",
    "pandas",
    "Cython",
    'astropy',
    'numba',
    'pyyaml'

]




#extra_files = find_data_files("pyipn/data")



with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = REQUIRED

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author=AUTHOR,
    author_email=EMAIL,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description=DESCRIPTION,
    install_requires=requirements,
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pyipn',
    name='pyipn',
    packages=find_packages(include=['pyipn'],exclude=('tests')),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/grburgess/pyipn',
    version='0.1.0',
    zip_safe=False,
    cmdclass={"upload": UploadCommand},
#        package_data={"": extra_files},
)









