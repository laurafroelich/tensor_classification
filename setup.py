from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
      name='tensor_classification',
      version='0.0.1',
      description='Code for tensor classification',
      license="MIT",
      long_description=long_description,
      author='Laura Froelich and Emil Froelich',
      author_email='laura.frolich@gmail.com',
      url="https://github.com/laurafroelich/tensor_classification",
      packages=['tensor_classification']  #same as name
      )
