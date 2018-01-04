from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='cpac_outputs_to_bids',
      version='0.1',
      description='Converts CPAC outputs into BIDS derivative format.',
      long_description=readme(),
      classifiers=['Development Status :: 3 = Alpha',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Medical Science Apps.'],
      keywords='functional MRI fMRI informatics preprocessing cpac',
      url='https://github.com/FCP-INDI/cpac_outputs_to_bids',
      author='R. Cameron Craddock',
      author_email='cameron.craddock@gmail.com',
      license='MIT',
      packages=['cpac_output_to_bids'],
      install_requires=[
          'numpy',
          'scipy',
          'pandas',
          'yaml'
      ],
      include_package_data=True,
      zip_safe=False)
