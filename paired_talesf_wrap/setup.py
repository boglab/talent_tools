from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(ext_modules=[Extension(
                   name="paired_talesf",
                   sources=["paired_talesf.pyx"],
                   language="c",
                   include_dirs = ["/usr/include/bcutils", "/usr/include/pairedtalesf"],
                   library_dirs = ["/usr/lib"],
                   libraries = ["bcutils", "pairedtalesf"],
                   )],
      cmdclass={'build_ext': build_ext})