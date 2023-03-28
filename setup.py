from setuptools import find_packages, setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='common_stuff',
    version="0.0.1",
    description="common functions to be used for data analysis",
    package_dir={"": "src"},
    author="Marco Di Gennaro",
    author_email="m.di.gennaro@outlook.com",
    url="https://mdigennaro.space",
    python_requires=">=3.7"
)