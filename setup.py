import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name='fragalysis_api',
    version='0.0.1.1',
    author='Fragment 5',
    author_email="",
    description="A package to load PDBs into fragalysis format.",
    long_description=long_description,  # README.md file as description
    long_description_content_type="text/markdown",
    url="https://github.com/xchem/fragalysis-api",
    packages=setuptools.find_packages(),
    install_requires=['biopython',
                      'numpy',
                      'pandas',
                      'pypdb'],  # Install requirements extracted from requirements.txt
    include_package_data=True,  # Allow to include other files than .py in package
    package_data={
        '': ['fragalysis_api/xcimporter/non_ligs.json']
    },  # Define which additional files should be included in package
    classifiers=[
        'Development Status :: 4 - Beta',
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    )
