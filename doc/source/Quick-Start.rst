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
* ``jobs/mstar01_C60.yaml``
    Job file, used to generate a disc around a 0.1 M_Sun star (Mah, Bitsch et al. 2023). Job files for other stellar masses are also included, check out the ``jobs`` folder! These job files are to be used with the corresponding config files because all the disk parameters are defined in the config files.
* ``jobs/mstar01_C60_moreCO.yaml``
    Job file, used to generate a disc around a 0.1 M_Sun star but with more CO and less CH4 in the disc (Mah, Bitsch et al. 2023). Use this with the corresponding config file.
* ``jobs/mstar01_C60_001.yaml``
    Job file, used to generate a less massive disc around a 0.1 M_Sun star (Mah, Bitsch et al. 2023). Use this with the corresponding config file.
