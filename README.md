# Minimal Function Cache

Minimal Function Cache is a repository for a proposed cut generation method in the thesis of Acadia Larsen which uses an explicit method select optimal cuts over a space of valid intersection cuts for a given MIP relaxation basis. 

The repository is in an alpha state. 

## Included Software and Tools

- Experimental Cut and Branch Solver using explicitly optimized cuts over the space of minimal functions with a finite number of breakpoints. 

- Code for generating representative elements of minimal functions for a semialgebraic parametric description of minimal functions in the Gomory Johnson Group Problem using high performance computing cluster.

- A cache of computed representative elements.

- Reproducible computational experiments.

# Goals and Non-Goals

## Goals
 - Illustrate concept of explicit optimal cut selection. 
 - Reproducible on local and HPC machines.
 - Mathematical correctness of cut generation.
 - Interface with current optimization software ``pyscipopt``.
 - Demonstrate use of ``passsagemath`` in application.
   
## Non-goals
 - Performance and code optimization; the cut generation technique presented is a proof of concept and is intended to function (for practical problem) with an excess of compute.
 - Documentation light, code follows thesis. 

# License 
The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation.
