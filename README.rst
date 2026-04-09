Minimal Function Cache
~~~~~~~~~~~~~~~~~~~~~~

Minimal Function Cache is a repository for a proposed cut generation method in the thesis of Acadia Larsen which uses an explicit method select optimal cuts over a space of valid intersection cuts for a given MIP relaxation basis. 

The repository is in an alpha state. User be ware.

Included Software and Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Experimental Cut and Branch Solver using explicitly optimized cuts over the space of minimal functions with a finite number of breakpoints. 

- Code for generating representative elements of minimal functions for a semialgebraic parametric description of minimal functions in the Gomory Johnson Group Problem using high performance computing cluster.

- A cache of computed representative elements.

- Reproducible computational experiments.

Goals and Non-Goals
~~~~~~~~~~~~~~~~~~~

 - Illustrate concept of explicit optimal cut selection.

 - Reproducible on local and HPC machines.

 - Mathematical correctness of cut generation.

 - Interface with current optimization software ``pyscipopt``.

 - Demonstrate use of ``passsagemath`` in application.
   
Non-goals
~~~~~~~~~

 - Performance and code optimization; the cut generation technique presented is a proof of concept and is intended to function (for practical problem) with an excess of compute.

 - Documentation light. 

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
