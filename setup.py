from setuptools import setup

setup(name='Lace',
      version='2.0.0',
      packages=['Lace'],
      license='GPL3',
      python_requires='>=3.12',
      install_requires=['pandas>=2.2',
                        'networkx>=3.4',
                        'numpy>=1.26',
                        'tqdm>=4.66'],
      entry_points={
          'console_scripts': [
              'BuildSuperTranscript=Lace.BuildSuperTranscript:main',
              'Lace_Checker=Lace.Checker:main',
              'Lace=Lace.Lace_run:main',
              'Mobius-as=Lace.Mobius_as:main',
              'Mobius=Lace.Mobius:main',
              'STViewer=Lace.STViewer:main'
          ]
      },
)
