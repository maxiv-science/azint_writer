from setuptools import setup, find_packages

setup(
    name="nxazint",
    description="A Python package for writing azimthual Integrated HDF5 files in NeXus format",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Meghdad Yazdi",
    author_email="meghdad.yazdi@maxiv.lu.se",
    url="https://github.com/maxiv-science/nxazint",
    packages=find_packages(),
    install_requires=[
        "h5py",
        "azint",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
