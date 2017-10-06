
import os
def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

from distutils.core import setup 
extra_files = package_files('pygear')
setup(name='pygear',
      version='0.25',
      package_data={'': extra_files},
      packages=['pygear'])

