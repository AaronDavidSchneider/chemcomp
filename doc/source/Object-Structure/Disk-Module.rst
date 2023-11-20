Disk Module
^^^^^^^^^^^^^^^^

The Disk module has all the physics of the protoplanetary disk included. It is possible to create your own disk (for examples of disks see ``chemcomp/disks`` folder).

Important Functions
###################

* **update()**

``update()`` is called in ``planet.update()`` and is the function that manages the disk timestep. ``update()`` updates the time and then calls ``update_parameters()`` which is defined in the child (per default static disk). For a good example see ``DiskLabDisk`` in ``kees.py``.

* **compute_viscous_evolution()**

Calculates the complete time evolution of the disk for this timestep. See functions there in.
The whole disk evolution is wrapped up here.