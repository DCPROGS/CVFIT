import pkg_resources
import numpy as np

from cvfit import data

def test_loading_from_excel():
    filename = "./Example/Example.xlsx"
    set0 = data.read_sets_from_Excel(filename, 2, 0, 2)[0]
    assert set0.X.all() == np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5]).all()

def test_imports():
    dependencies = [
      'scipy',
      'numpy',
      'xlrd',
      'markdown',
      'matplotlib',
    ]
    # throw a DistributionNotFound or VersionConflict
    # exception if a dependency is not met
    pkg_resources.require(dependencies)