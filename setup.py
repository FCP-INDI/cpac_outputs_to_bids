from subprocess import call

from setuptools import Command, find_packages, setup


def readme():
    with open('README.md') as f:
        return f.read()


class RunTests(Command):
    """Run all tests."""
    description = 'run tests'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        """Run all tests!"""
        errno = call(['py.test', '--cov=cpac_output_to_bids', '--cov-report=term-missing'])
        raise SystemExit(errno)


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
          'pyyaml'
      ],
      extras_require={
          'test': ['coverage', 'pytest', 'pytest-cov'],
      },
      entry_points={
          'console_scripts': [
              'cpb=cpac_output_to_bids.cli:main',
          ],
      },
      cmdclass={'test': RunTests},
      include_package_data=True,
      zip_safe=False)
