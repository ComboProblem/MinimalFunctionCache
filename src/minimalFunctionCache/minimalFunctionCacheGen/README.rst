Minimal Function Cache Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instructions for adding to the function cache.

Requirements
~~~~~~~~~~~~

- cluster with ``SLURM``.

- cluster with ``apptainer`` (https://apptainer.org/).

- a terminal text editor (``vim``, ``nano``, ``emacs``, ect).

- ``git``

Generating the Cache
~~~~~~~~~~~~~~~~~~~~

- Clone the MinimalFunctionCache repository using git.::

    git clone https://github.com/ComboProblem/MinimalFunctionCache.git

- Edit the cache job parameters script using editor of choice.::

    vim ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh

- Compile and run the job generation and submission script.::

    chmod +x ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cluster_cache_gen.sh
    ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cluster_cache_gen.sh

- When finished; you can use git to make a PR to submit the additionally generated function cache or move to desired repository. An example to make a PR to the main MinimalFunctionCache repository.::

    cd MinimalFunctionCache
    git add .
    git commit -m "This is an example commit message"
    git push
