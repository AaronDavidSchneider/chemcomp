slurm.yaml
^^^^^^^^^^

The ``chemcomp_pipeline_slurm`` script uses ``dask-jobqueue`` to distribute the emberessing parallel workload of the grid simulations into multiple slurm jobs.
The script creates and destroys individual workers that are distributed over the cluster.
This way you can utilize the full capacity of the cluster for short periods, while only blocking as little as nescessary. You don't need to use long walltimes for individual jobs, since individual jobs are recreated when the walltime expires.

The ``slurm.yaml`` file is used to configure the individual slurm jobs that are used to distribute the parameter grid of individual ``chemcomp`` simulations to a variety of slurm jobs.

.. code-block :: yaml

   max_jobs:  384
   min_jobs:  1
   cores: 2
   memory: '10GB'
   queue: 'astro_short'
   walltime: '24:00:00'
   max_duration_of_worker: '11h'
   interface: 'bond0'


Explanation
"""""""""""

The general philosphy of the ``chemcomp_pipeline_slurm`` script is to use many individual slurm jobs that are distributed to the cluster.
This can be achieved by using few cores per job (e.g., 2), short walltimes per job and many individual jobs (high value for ``max_jobs``).

+----------------------------+--------------------------------------------------------------------------------------------------------+
| Parameter                  | Usage                                                                                                  |
+============================+========================================================================================================+
| ``max_jobs``               | Maximum number of individual slurm jobs used for execution                                             |
+----------------------------+--------------------------------------------------------------------------------------------------------+
| ``min_jobs``               | Minimum number of individual slurm jobs used for execution (should be one for most cases).             |
+----------------------------+--------------------------------------------------------------------------------------------------------+
| ``cores``                  | Number of cores per slurm job. Choose low number, to fit better into the cluster.                      |
+----------------------------+--------------------------------------------------------------------------------------------------------+
| ``memory``                 | memory per individual slurm job.                                                                       |
+----------------------------+--------------------------------------------------------------------------------------------------------+
| ``queue``                  | partition/queue used for individual slurm jobs. Use a short one, ``dask`` creates new jobs, if needed. |
+----------------------------+--------------------------------------------------------------------------------------------------------+
| ``walltime``               | walltime associated with the queue                                                                     |
+----------------------------+--------------------------------------------------------------------------------------------------------+
| ``max_duration_of_worker`` | time after which the individual jobs are canceled. ``dask`` will simply create new ones.               |
+----------------------------+--------------------------------------------------------------------------------------------------------+
| ``interface``              | Interface between nodes, use ``ib0`` if your cluster has infiniband                                    |
+----------------------------+--------------------------------------------------------------------------------------------------------+

:Note: Depending on your cluster, you might need to wrap the ``chemcomp_pipeline_slurm`` script in a  slurm job. This is especcially the case, if your login nodes are not allowed to communicate with the compute nodes.

:Note: Your execution did not finish? Don't worry. You can continue your runs with ``-o 0``, which will only runs on those parameters again, that did not yet finish.
