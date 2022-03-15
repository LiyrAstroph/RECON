import os
from setuptools import setup
from setuptools.extension import Extension
import pkgconfig
from glob import glob
import numpy

# if CC is not set, use the default value
if not os.environ.get("CC"):
  os.environ["CC"] = "mpicc"

def configure_mpi():
  """
  get configurations of mpi
  """
  if pkgconfig.exists('mpich'):
    mpiconf = pkgconfig.parse('mpich')
  elif pkgconfig.exists('ompi'):
    mpiconf = pkgconfig.parse('ompi')
  else:
    raise SystemError("Not found MPICH or OpenMPI installed.")

  return mpiconf

mpiconf = configure_mpi()

basedir = os.path.dirname(os.path.abspath(__file__))
homedir = os.environ['HOME']
cdnestdir = "/home/liyropt/Projects/GIT/CDNest"
include_dirs = [basedir, os.path.join(basedir, "src"), cdnestdir, numpy.get_include(),] + mpiconf['include_dirs']
library_dirs = [basedir, cdnestdir] + mpiconf['library_dirs']

src = [os.path.join(basedir, "python", "cyrecon", "cyrecon.pyx")] + glob(os.path.join(basedir, "src", "*.c"))

headerfiles = [os.path.join(basedir, "python", "cyrecon", "cyrecon.pyd")] + glob(os.path.join(basedir, "src", "*.h"))

if os.name == 'nt':  # Windows, assumming MSVC compiler
  libraries = ['dnest']
  compiler_args = ['/Ox', '/fp:fast']
  link_args = []
elif os.name == 'posix':  # UNIX, assumming GCC compiler
  libraries = ['m', 'c', 'gsl', 'gslcblas', 'fftw3', 'dnest'] + mpiconf['libraries']
  compiler_args = ['-O3', '-ffast-math']
  link_args = []


try:
  from Cython.Build import cythonize
except ImportError:
  raise RuntimeError('Cython not found.')

extensions = cythonize([
  Extension("cyrecon.cyrecon", 
	  sources=src,
    depends=headerfiles,
	  extra_compile_args=compiler_args,
    extra_link_args=link_args,
    include_dirs=include_dirs,
    libraries=libraries,
    library_dirs=library_dirs
    ),
  ], annotate=False)

setup(
	name="cyrecon",
  version="0.1.0",
	packages=["cyrecon"],
  package_dir={"":"python"},
	ext_modules = extensions,
  description = 'A package for measuring spectral power and reconstructing time series in AGN',
  author = 'Yan-Rong Li',
  author_email = 'liyanrong@mail.ihep.ac.cn',
  setup_requires=['numpy', 'mpi4py', 'pkgconfig'],
  install_requires=['numpy', 'mpi4py', 'pkgconfig'],
  license="GSL",
	)
