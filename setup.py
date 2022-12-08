import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

    setuptools.setup(
    name="henipipe",
    version="2.2.1",
    author="Scott Furlan",
    author_email="scottfurlan@gmail.com",
    description="A python wrapper for fast and parallel processing of sequencing data using CutnRun or CutnTag",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/scfurl/henipipe.git",
    packages=setuptools.find_packages(),
    package_data={'henipipe': ['henipipe/data/genomes.json']},
    install_requires=['six >= 1'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.5',
    entry_points={'console_scripts': [
        'henipipe = henipipe.__main__:run_henipipe',
        'samTobed = henipipe.samTobed:run_sam2bed',
        'pyWriter = henipipe.pyWriter:run_pyWriter',
        'auc = henipipe.auc:run_auc',
    ]},
    )


"""
#TEST PyPI
## run this to make package: python3 setup.py sdist bdist_wheel
## run this to upload to test pypi: twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#TEST INSTALL
# pip install --index-url https://test.pypi.org/henipipe/ henipipe


#FOR REAL
cd ~/computation/develop/henipipe/
rm dist/*
y
git commit -a -m "update environs for cluster GCC 11.2"
git push
python setup.py sdist bdist_wheel
python -m twine upload dist/*
scfurl
rm -R dist/* build/* henipipe.egg-info/*


##Install pipx
## python3 -m pip install --user pipx
## python3 -m pipx ensurepath

**At SCRI do the following**

module load python
python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install --include-deps --pip-args '--trusted-host pypi.org --trusted-host files.pythonhosted.org' henipipe
pipx install --spec git+https://github.com/scfurl/henipipe --include-deps henipipe --pip-args '--trusted-host pypi.org --trusted-host files.pythonhosted.org'
pipx install --spec git+https://github.com/scfurl/henipipe@cleaner --include-deps henipipe --pip-args '--trusted-host pypi.org --trusted-host files.pythonhosted.org'


**At the FHCRC do the following...**

module load Python/3.6.7-foss-2016b-fh1
python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install --include-deps henipipe
pipx install git+https://github.com/scfurl/henipipe --include-deps
pipx uninstall henipipe


pipx uninstall henipipe
pipx install git+https://github.com/scfurl/henipipe --include-deps


## running test data
ml Python
pipx uninstall henipipe
pipx install henipipe
git clone https://github.com/scfurl/henipipe.git
cd ~/henipipe/test_data
mkdir henipipe
cd henipipe
henipipe MAKERUNSHEET -fq ../fastq
awk -F ',' '{print $1, $2}' runsheet.csv

#proceed with henipipe steps
henipipe ALIGN -t 16 -r runsheet.csv -n spike_in
henipipe SCALE -r runsheet_fixed.csv -n spike_in
henipipe SEACR -r runsheet_fixed.csv


"""