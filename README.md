# KAPOW - Kapow is A Pythonic OpenMM Wrapper

## Introduction

Kapow is an implementation of an adapter class designed to provide a
more Pythonic interface to the OpenMM molecular simulation
toolkit. Taking advantage of the dynamic properties of the Python
programming language this module mostly aims to circumvent the use of
getters and setters in OpenMM objects, replacing them with the use of
properties and more sophisticated descriptors as required.

Apologies for the project name. In case you're wondering its a
[recursive acronymn](https://en.wikipedia.org/wiki/Recursive_acronym)
(like GNU). Whilst common in the open source world, it is my firm
belief they should be more widely used in scientific software.

## Installation

At present this is very basic. The simtk.openmm package must already
be installed so configure your Python environment as required,
e.g. load any desired virtual environments. Download the source code
of the project and from it's root directory run: 
> python setup.py install

If you encounter permission problems, an installation into your home
directory may also be requested with:
> python setup.py install \-\-user

Test the installation by changing into the tests directory and
running:
> python test\_adapter.py

## Getting Started

An introductory ramble is available as demo.ipynb in the project's
root directory.

## Caveats

This is an experimental early alpha release subject to rapid changes
and guarantee of a fixed API.

Kapow has not yet been widely tested against different OpenMM
versions. This release is tailored to OpenMM 7.2.1, other versions may
not pass all the tests or even import successfully.

