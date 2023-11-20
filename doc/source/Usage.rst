Usage
-----

Install ``chemcomp`` (see :ref:`Installation`)

Activate environment:

.. code-block:: bash

   conda activate chemcomp


create a run directory (e.g., chemcomp_runs) somewhere, where you want to run and store your simulations.

.. code-block:: bash

   mkdir chemcomp_runs
   cd chemcomp_runs


Configure your run and get it running.

.. code-block:: bash

   # a single run config:
   chemcomp_main -c <path to config file>

   # probing multiple parameters:
   chemcomp_pipeline -c <path to config file> -j <path to job file> -d 0 -o 0

   # probing multiple parameters with slurm:
   chemcomp_pipeline_slurm -c <path to config file> -j <path to job file> -d 0 -o 0 -s <slurm file>

Head over to see  :ref:`Overview of files` to understand the run scripts!

Running chemcomp will create the following folders (if not existing) in the current working directory:

.. code-block:: bash

   chemcomp_run
   ├── config              # directory where the config files of the mesh will be stored
   ├── output              # output files
   ├── zip_output          # Storage of the ziped output files

If you want to know how to open and plot these files, head over to :ref:`Output`!
