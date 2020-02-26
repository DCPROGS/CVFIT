from setuptools import setup

setup(name='cvfit',
      version='1.0.0',
      description="DC's lab curve fitting tools",
      url='https://github.com/DCPROGS/CVFIT',
      keywords='Hill fit',
      author='Remis Lape',
      author_email='',
      license='GPL-2.0',
      packages=['cvfit'],
      install_requires=[
          'numpy',
          'pandas',
          'markdown',
          'matplotlib',
          'PyQt5',
          'xlrd',
      ],
      zip_safe=False)