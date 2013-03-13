from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("gtf", ["src/gtf.pyx", "src/gtf.pxd"])]

setup(
  name='pygtf',
  cmdclass={'build_ext': build_ext},
  ext_modules=ext_modules,
  setup_requires=["cython==0.18"],
  install_requires=["argparse"],
  entry_points={
    'console_scripts': [
        'flux-utils = commands:main'
    ]
  }
)
