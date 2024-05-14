# setup.py

from setuptools import setup, find_packages

setup(
    name="oer-db-tools",
    version="0.1.0",
    author="Wei",
    author_email="wei.chen@uclouvain.be",
    description="Toolkit for Mongo databases dedicated to OER workflows",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/wch3n/oer-db-tools",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    install_requires=[
        "jobflow",
        "pymatgen",
    ],
)
