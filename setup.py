from setuptools import setup, find_packages


setup(name='survival',
      version='0.0.0',
      description='Survival analysis.',
      url='https://github.com/ihmeuw-msca/Survival',
      author='ihmeuw-msca',
      author_email='ihme.math.sciences@gmail.com',
      license='MIT',
      package_dir={'': 'src'},
      packages=find_packages(where='src'),
      install_requires=['numpy',
                        'scipy',
                        'pandas',
                        'pytest'],
      zip_safe=False)
