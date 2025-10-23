"""
Setup script for Komega Python Library

This script allows the Komega Python library to be installed as a package
using pip or other Python package managers.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

import os

from setuptools import find_packages, setup


# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), "README.md")
    if os.path.exists(readme_path):
        with open(readme_path, "r", encoding="utf-8") as f:
            return f.read()
    return "Komega Python Library - A library for solving linear systems in materials science"


# Read requirements
def read_requirements():
    requirements_path = os.path.join(os.path.dirname(__file__), "requirements.txt")
    if os.path.exists(requirements_path):
        with open(requirements_path, "r", encoding="utf-8") as f:
            return [
                line.strip() for line in f if line.strip() and not line.startswith("#")
            ]
    return ["numpy>=1.19.0", "scipy>=1.7.0"]


setup(
    name="komega-python",
    version="2.0.0",
    author="Kazuyoshi Yoshimi",
    author_email="k-yoshimi@issp.u-tokyo.ac.jp",
    description="A Python library for solving linear systems in materials science",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/k-yoshimi/Komega",
    project_urls={
        "Bug Tracker": "https://github.com/k-yoshimi/Komega/issues",
        "Documentation": "https://github.com/k-yoshimi/Komega/tree/python/python_komega",
        "Source Code": "https://github.com/k-yoshimi/Komega/tree/python/python_komega",
    },
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.10.0",
            "flake8>=3.8.0",
            "black>=21.0.0",
            "isort>=5.0.0",
            "mypy>=0.910",
        ],
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
        "mpi": [
            "mpi4py>=3.0.0",
        ],
        "profiling": [
            "line-profiler>=3.0.0",
            "memory-profiler>=0.60.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "komega-test=run_tests:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="linear-systems, iterative-solvers, materials-science, fortran, bicg, cg, cocg",
)
