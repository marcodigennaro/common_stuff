from setuptools import setup
from setuptools import find_packages

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='mylib',
    version='0.1.0',
    description='a demonstration of a python package',
    long_description=long_description,
    author='Marco Di Gennaro',
    author_email='m.di.gennaro@outlook.com',
    readme='README.md',
    url='https://github.com/marcodigennaro/common_stuff',
    packages=find_packages(exclude=('test*'))
)


if __name__ == '__main__':
    exit(setup())
