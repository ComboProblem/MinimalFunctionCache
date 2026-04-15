Minimal Function Cache
~~~~~~~~~~~~~~~~~~~~~~
``minimalFunctionCache`` is a optional package for ``cutgeneratingfunctiology`` which contains necessary data for using pre-computed cell descriptions of the space of continuous minimal functions with at most :math:`k` breakpoints where :math:`k=7`. This repository is in a beta state.

This package is based on the dissertation of Acadia Larsen. 

Notes
~~~~~

- Data is written in .cvs files.
- pip installable from cloned repository.
- Tools to (re)generate data using a HPC.
- Practical evidence for polyedral computational speed ups provided by ``pplite`` as compared to ``ppl``.

Documentation
~~~~~~~~~~~~~
Documentation is minimal as this repository is intended to be a purely optional data source for ``cutgeneratingfunctionology``.

Installation
~~~~~~~~~~~~

- Clone the GitHub repository https://github.com/ComboProblem/cutgeneratingfunctionology/tree/MinFunStable.git and https://github.com/ComboProblem/MinimalFunctionCache.git::

    git clone https://github.com/ComboProblem/cutgeneratingfunctionology/tree/MinFunStable.git
    git clone https://github.com/ComboProblem/MinimalFunctionCache.git
    cd cutgeneratingfunctionology

- Create a virtual environment::

    python3 -m venv venv-cgf
    source venv-cgf/bin/activate

- Install the cutgeneratingfunctionology package using pip::

    pip install ".[passagemath]"

- Install the MinimalFunctionCache using pip::

    cd MinimalFunctionCache
    pip install .

- Install pplite::

    pip install pplitepy

- Start using the function cache in sagemath::

    sage
    sage: from cutgenerationfunctionology.igp import PiMinContConatiner
    sage: Pi5 = PiMinContConatiner(5) # loaded with the function cache!
    sage: Pi5_cell_description = [cell for cell in Pi5.get_semialgebraic_sets()]

License 
~~~~~~~
The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation.
