from setuptools import setup

with open("README.md", "r") as f:
    long_description = f.read()

if __name__ == "__main__":
    setup(
        name='mylib',
        version="0.0.1",
        description="common functions to be used for data analysis",
        package_dir={"": "src"},
        author="Marco Di Gennaro",
        author_email="m.di.gennaro@outlook.com",
        url="https://mdigennaro.space",
        python_requires=">=3.7"
    )

