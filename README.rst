Minimal Function Cache
~~~~~~~~~~~~~~~~~~~~~~

Minimal Function Cache is a repository containing two packages python packages ``minimalFunctionCache`` and ``parametricCutGen``.

``minimalFunctionCache`` is a optional package for ``cutgeneratingfunctiology`` which contains precomputed cell descptions of minimal functions with at most :math: `k` breakpoints where :math: `k=6`. This repository is in a beta state.

``parametricCutGen`` is a package which implemments a single row optiomal cut selection for Mixed Integer Programs over the domain (and restricted domains) of continuous minimal functions with at most :math: `k` breakpoints. This repository is in an alpha state.

These packages are based on the disseration of Acadia Larsen. 

Included Software and Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``minimalFunctionCache`` includes precomputed cell descriptions of minimal functions, code for generating the function cache on a HPC, and a comparison of ``ppl`` and ``pplite`` as polyhedral backends.  

- ``parametricCutGen`` is implemented to work with ``pyscipopt`` the python wrapper of ``SCIP``. This package includes a solver for cut selection problems over the domain (with possibly some constraints) of continuous minimal functions, an interface to ``pyscipopt`` and reproduceable parametric experinments for testing cut selection.


Goals and Non-Goals
~~~~~~~~~~~~~~~~~~~

 - Illustrate concept of explicit optimal cut selection as a proof of concept for MIP solvers. 

 - Reproducibliblity of experimental data.

 - Mathematical correctness of cut generation up to some :math: `(M,\epsilon)` parameters. 

 - Demonstrate use of ``passsagemath`` in application; in particular illustrate application of cutting edge mathematics to application of MIPs.
   
Non-goals
~~~~~~~~~

 - Performance and code optimization; the cut generation technique presented is a proof of concept and is intended to function (for practical problem) with an excess of compute.

 - Documenation; documentation and testing is minimal. What is written is what is requried to reproduce the results for an expert user. 

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

    cd ..
    sage
    sage: from cutgenerationfunctionology.igp import PiMinContConatiner
    sage: MinFun_with_at_most_5_breakpoints = PiMinContConatiner(5) # loaded with the function cache!

License 
~~~~~~~
The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation.
