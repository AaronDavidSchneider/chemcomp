Disk only runs
""""""""""""""
If you want to use ``chemcomp`` to run a disk model and don't care about a planet, you may use the following settings:

In the ``config_planet`` section:

* Set ``model`` to ``NoAccretion``. This will introduce a passive planet that ends the simulation, once planted into the disk.
* Set ``t_0`` to the value for how long you want to evolve the disk, before the planet ends the simulation.

In the ``defaults`` section:

* Set ``save_disk`` to ``True`` to make sure the disk is stored
* Set ``save_interval`` and ``save_disk_interval`` to a reasonable interval for which you require the disk quantities to be stored
