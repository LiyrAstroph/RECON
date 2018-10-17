import os
from setuptools import setup
from setuptools.extension import Extension

os.environ["CC"] = "mpicc"

include_dirs = ["/home/liyropt/Projects/GIT/DNest","/home/liyropt/Projects/LCRecon/Code"]
library_dirs = ["/home/liyropt/Projects/GIT/DNest","/home/liyropt/Projects/LCRecon/Code"]

if os.name == 'nt':  # Windows, assumming MSVC compiler
  libraries = ['dnest']
  compiler_args = ['/Ox', '/fp:fast']
elif os.name == 'posix':  # UNIX, assumming GCC compiler
  libraries = ['m', 'dnest', 'recon']
  compiler_args = ['-O3', '-ffast-math']


try:
  from Cython.Build import cythonize
except ImportError:
  raise RuntimeError('Cython not found.')

extensions = cythonize([
  Extension("cyrecon", 
	  sources=["cyrecon.pyx",],
	  extra_compile_args=compiler_args,
    include_dirs=include_dirs,
    libraries=libraries,
    library_dirs=library_dirs
    ),
  ], annotate=False)

setup(
	name="cyrecon",
	packages="cyrecon",
	ext_modules = extensions,
  description = 'A package for measuring spectral power and reconstructing time series in AGN',
  author = 'Yan-Rong Li',
  author_email = 'liyanrong@mail.ihep.ac.cn',
	)