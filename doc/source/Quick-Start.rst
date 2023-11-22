Quick Start
-----------

Install ``chemcomp`` (see :ref:`Installation`)

Activate environment:

.. code-block:: bash

   conda activate chemcomp


create a run directory (e.g., ``chemcomp_runs``) somewhere, where you want to run and store your simulations.

.. code-block:: bash

   mkdir chemcomp_runs
   cd chemcomp_runs


Configure your run and get it running.

.. code-block:: bash

   # a single run config:
   chemcomp_main -c <path to config file>

   # probing multiple parameters:
   chemcomp_pipeline -c <path to config file> -j <path to job file> -d 0 -o 0

Head over to see  :ref:`Configuration` to understand the run scripts!

Running chemcomp will create the following folders (if not existing) in the current working directory:

.. code-block:: bash

   chemcomp_run
   ├── config              # directory where the config files of the mesh will be stored
   ├── output              # output files
   ├── zip_output          # Storage of the ziped output files

If you want to know how to open and plot these files, head over to :ref:`Output`!

Examples
^^^^^^^^

You can find a few examples in the ``jobs`` and ``config`` folder.
These examples will help you to understand, how ``chemcomp`` has been used to create some figures in the original publications.

These examples include:

* ``jobs/paper_1_disk.yaml``
    Job file, used to generate the disk models in the first chemcomp paper (Schneider & Bitsch 2021a)
* ``jobs/paper_1_grid_thorngren.yaml``
    Job file, used to generate the grid of planets to compare the heavy element content in the first chemcomp paper (Schneider & Bitsch 2021a)
* ``jobs/paper_1_grid_thorngren.yaml``
    Job file, used to generate the growthtracks of planets in the first chemcomp paper (Schneider & Bitsch 2021a)
* ``jobs/paper_2_Jupiter.yaml``
    Job file, used to generate Jupiter like planets for the second chemcomp paper (Schneider & Bitsch 2021b)