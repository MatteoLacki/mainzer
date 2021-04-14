# This Python file uses the following encoding: utf-8
from setuptools import setup, find_packages

setup(  name='mainzer',
        packages=find_packages(),
        version='0.0.1',
        description='Description.',
        long_description='Long description.',
        author='MatteoLacki',
        author_email='matteo.lacki@gmail.com',
        url='https://github.com/MatteoLacki/mainzer.git',
        keywords=['top down mass spectrometry', 'intensity annotation'],
        classifiers=['Development Status :: 1 - Planning',
                     'License :: OSI Approved :: BSD License',
                     'Intended Audience :: Science/Research',
                     'Topic :: Scientific/Engineering :: Chemistry',
                     'Programming Language :: Python :: 3.6',
                     'Programming Language :: Python :: 3.7'],
        install_requires=['numpy',
                          'pandas >= 1.2',
                          'IsoSpecPy',
                          'networkx',
                          'pyteomics',
                          'lxml',
                          'ncls',
                          'tqdm',
                          'scipy',
                          'pyteomics',
                          'matplotlib',
                          'aa2atom',
                          'nogui'],
        scripts = [
            'bin/lipido.py'
        ]
)
