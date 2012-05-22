from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(ext_modules=[Extension(
                   name="talesf",
                   sources=["talesf.pyx"],
                   language="c",
                   #include_dirs = ["../src"],
                   #library_dirs = ["../"],
                   include_dirs = ["/usr/include"],
                   library_dirs = ["/usr/lib"],
                   libraries = ["talesf"],
                   )],
      cmdclass={'build_ext': build_ext})