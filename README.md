# Minimal Function Cache

Minimal Funciotn Cache is a repository for a proposed cut generation method in the thesis of Acadia Larsen which uses an explicit method select optimal cuts over a space of valid intersection cuts for a given MIP relaxation basis. 

The repository is in an alpha state. 

## Included Software and Tools

- Experimental Cut and Branch Solver using explicitly optimized cuts over the space of minimal functions with a finite number of breakpoints. 

- Code for generating represenative elements of minimal functions for a semialgebraic parametric description of minimal functions in the Gomory Johnson Group Problem using high preformance computing cluster.

- A cache of computed repersenative elements.

- Reproducable computational experiments.

# Goals and Non-Goals

## Goals
 - Illustrate concept of expicit optimal cut selection. 
 - Reproducablity on local and HPC machines.
 - Mathematical correcntness of cut generation.
 - Interface with current optimization software ``pyscipopt``.
 - Demonstrate use of ``passsagemath`` in application.
   
## Non-goals
 - Preformance and code optimization; the cut generation technique presented is a proof of concept and is intended to function (for practical problem) with an excess of compute.
 - Documentaiton light, code follows thesis. 

# Installation
TBD

# Lisence 
The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation.

The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation.
