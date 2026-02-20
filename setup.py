"""Lace package setup."""
from setuptools import setup

setup(name='Lace',
      version='2.0.0',
      packages=['lace'],
      license='GPL3',
      python_requires='>=3.12',
      install_requires=['pandas>=2.2',
                        'networkx>=3.4',
                        'numpy>=1.26',
                        'tqdm>=4.66'],
      entry_points={
          'console_scripts': [
              'BuildSuperTranscript=lace.build_supertranscript:main',
              'Lace_Checker=lace.checker:main',
              'Lace=lace.lace_run:main',
              'Mobius-as=lace.mobius_as:main',
              'Mobius=lace.mobius:main',
              'STViewer=lace.st_viewer:main'
          ]
      },
)
