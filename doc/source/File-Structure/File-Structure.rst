Overview of files
-----------------

.. code-block :: bash

   chemcomp
   ├── analysis              # sample python scripts for post-processing
   ├── chemcomp              # main directory of the code
   │   ├── accretion_models  # modules for accretion
   │   ├── disks             # disk modules
   │   │   ├── _0pacity      # opacity files
   │   ├── helpers           # helper (e.g. wrapper for multiprocessing, etc)
   │   │   └── units         # units and other constants
   │   └── planets           # planet modules
   ├── config                # config files
   ├── jobs                  # job files
   ├── output                # output files
   └── zip_output            # ziped output files


There are several ways how you can use ``chemcomp``.
Upon installation of ``chemcomp``, several scripts are added to your python environment.
We will discuss these in more detail below.


.. toctree::
   :maxdepth: 1
   :caption: Run configuration:

   Config-File
   Job-File
   Slurm-File


chemcomp_main
^^^^^^^^^^^^^

Wrapper for running a single configuration. Execute with:

.. code-block :: bash

   chemcomp_main [-h/--help][-c/--config_file]


* ``-c`` specifies the path to the config file.
* ``-h`` shows the help dialog.

Parameters that are set in the config file can be found :ref:`here <config.yaml>`.

chemcomp_pipeline
^^^^^^^^^^^^^^^^^

Wrapper for running multiple configurations. Uses ``dask-jobqueue``. Execute with:

.. code-block :: bash

   chemcomp_pipeline [-h/--help][-c/--config_file][-d/--delete][-j/--job_file][-o/overwrite]

* ``-c`` specifies the path to the config file.
* ``-o`` specifies if you want to overwrite existing simulations. Use ``-o 1`` if you want to overwrite all your already finished runs.
* ``-h`` shows the help dialog.
* ``-d`` specifies wether to delete (use ``-d 1``) or not delete (use ``-d 0``) the files in ``output`` after zipping.
* ``-j`` specifies the path to the job file.

Explanations about the config file can be found :ref:`here <config.yaml>`.
Explanations of the job file can be found :ref:`here <job.yaml>`.

:Note: Your execution did not finish? Don't worry. You can continue your runs with ``-o 0``, which will only runs on those parameters again, that did not yet finish.

chemcomp_pipeline_slurm
^^^^^^^^^^^^^^^^^^^^^^^

Unlike ``chemcomp_pipeline``, ``chemcomp_pipeline_slurm`` is build for execution on a slurm cluster.
This run script takes identical arguments, but adds an additional flag for a file with slurm arguments:

.. code-block :: bash

   chemcomp_pipeline_slurm [-h/--help][-c/--config_file][-d/--delete][-j/--job_file][-o/overwrite][-s/--slurm_file]

* ``-s`` specifies the path to the slurm file.

Explanations about the slurm file can be found :ref:`here <slurm.yaml>`.


The chemcomp directory
^^^^^^^^^^^^^^^^^^^^^^

The ``chemcomp/chemcomp`` directory is the directory where the main code resides.

.. code-block :: bash

   chemcomp/chemcomp
   ├── __init__.py
   ├── accretion_models
   │   ├── __init__.py
   │   ├── _accretion_class.py     # framework for accretion models
   │   ├── gas.py                  # gas accretion model
   │   ├── pebbles.py              # pebble accretion model
   │   └── planetesimals.py        # planetesimal accretion model
   ├── disks
   │   ├── _0pacity
   │   │   ├── __init__.py
   │   │   ├── meanopac1_5050.dat  # opacity lookup table
   │   │   ├── meanopac5.dat       # opacity lookup table
   │   │   └── opacity_models.py   # delivers the opacities
   │   ├── __init__.py
   │   ├── _chemistry.py           # Bitsch2020 composition model
   │   ├── _chemistry_new.py       # SchneiderBitsch2020 composition model
   │   ├── _disk_class.py          # general disk model (parent of individual disks)
   │   ├── aaron.py                # analytical visc heating. Do not use.
   │   ├── bert.py                 # weird disk from bert
   │   ├── ida.py                  # weird Ida2016 disk
   │   ├── kees.py                 # visc heating that works
   │   ├── micha.py                # Michas disk. Do not use.
   │   ├── mmsn.py                 # simple mmsn disk with Hayashi H/R.
   │   └── twopop.py               # disk used in twopop.
   ├── helpers
   │   ├── __init__.py
   │   ├── analysis_helper.py      # framework for post processing
   │   ├── debugger_helper.py      # some useful functions for debugging during runtime
   │   ├── import_config.py        # functions that are used for the import of the parameters in the config file.
   │   ├── main_helpers.py         # functions that are used for main.py and pipeline.py.
   │   ├── solvediffonedee.py      # functions used for visc disk evolution
   │   ├── tridag.py               # functions used for visc disk evolution
   │   └── units
   │       ├── __init__.py
   │       ├── berts_units.py      # Berts weird code units. Not used.
   │       └── chemistry_const.py  # constants used for the _chemistry_new and _chemistry
   └── planets
       ├── __init__.py
       ├── _planet_class.py        # main planet module (parent).
       ├── bertmigration.py        # P11 + K18 migration + feedback torques.
       ├── migration_1.py          # only Lindblad torque.
       ├── no_accretion.py         # no accretion, no migration -> only disk
       └── simple_planet.py        # no migration


The code is structured in four main modules:

+-------------------------------------+-------------------------------------+------------------------------------------------------------------------------+
|              Module                 |   directory of parent class         |                                   Used for                                   |
+=====================================+=====================================+==============================================================================+
| :ref:`Disk <Disk Module>`           |         ``_disk_class``             | Deals with the disk related physics (e.g. dust growth, visc evolution, etc.) |
+-------------------------------------+-------------------------------------+------------------------------------------------------------------------------+
| :ref:`Chemistry <Chemistry Module>` | ``_chemistry_new`` / ``_chemistry`` |                     Chemical compositions used in ``Disk``.                  |
+-------------------------------------+-------------------------------------+------------------------------------------------------------------------------+
| :ref:`Planet <Planet Module>`       |        ``_planet_class``            |                Controls the other modules. Grows and migrates.               |
+-------------------------------------+-------------------------------------+------------------------------------------------------------------------------+
| :ref:`Accretion <Accretion Module>` |       ``_accretion_class``          |               Used in ``Planet``. Calculates the accretion rates.            |
+-------------------------------------+-------------------------------------+------------------------------------------------------------------------------+


Other directories
"""""""""""""""""

The paths for ``output``, ``zip_output``, ``config`` and ``jobs`` can be adjusted in ``chemcomp/chemcomp/helpers/__init__.py``

+----------------+---------------------+--------------------------------------------------------+
| Path           | Variable            | Use of path                                            |
+================+=====================+========================================================+
| ``output``     | ``OUTPUT_PATH``     | storage path for the output files                      |
+----------------+---------------------+--------------------------------------------------------+
| ``zip_output`` | ``ZIP_OUTPUT_PATH`` | storage path for the zipped output files               |
+----------------+---------------------+--------------------------------------------------------+
| ``jobs``       | ``JOB_PATH``        | path in which job configurations are stored            |
+----------------+---------------------+--------------------------------------------------------+
| ``config``     | ``CONFIG_PATH``     | path in which configuration/parameter files are stored |
+----------------+---------------------+--------------------------------------------------------+

The complete set of parameters in the config file is explained :ref:`here <config.yaml>`.
Explanations of the job file can be found :ref:`here <Job.yaml>`.
