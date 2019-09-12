import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scfurl",
    version="0.1",
    author="Scott Furlan",
    author_email="scottfurlan@gmail.com",
    description="A python wrapper for fast processing of sequencing data using CutnRun or CutnTag",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/scfurl/henipipe.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

## run this to make package: python3 setup.py sdist bdist_wheel

## run this to upload to pypi: python3 -m twine upload dist/*
