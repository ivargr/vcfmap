
from setuptools import setup

setup(name='vcfmap',
      version='0.0.1',
      description='Connecting vcf entries and haplotype information to graph edges',
      url='http://github.com/ivargr/vcfmap',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=['vcfmap'],
      zip_safe=False,
      install_requires=['numpy', 'tqdm'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points = {
            'console_scripts': ['vcfmap=vcfmap.command_line_interface:main'],
      })