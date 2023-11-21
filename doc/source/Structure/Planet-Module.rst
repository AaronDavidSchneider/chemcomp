Planet Module
^^^^^^^^^^^^^

Supervisor module that handles the growth of the planet (see ``grow_mass()``) and updates all other modules. It is possible to overwrite methods (for examples of planets see ``chemcomp/planets`` folder). The most advanced planet is certainly ``BertPlanet``, since it uses all kinds of migration techniques.

The `Planet`` module also handles the outputing which is handled by the ``DataObject`` (see :ref:`Output`) class.

Important Functions
###################

* **grow_mass()**

IMPORTANT: This is the Supervisor! From here everything else is directed. This should be the only method that loops over the time.

* **update_all()**

Handles the interpolation of disk quantities and the migration of the planet.

* **update()**

all the interpolation steps are done here (see also [Operating Principle](Operating-Principle)).

* **update_a_p()**

Can be used to calculate the new position of the planet based on torques calculated in ``update_torque()``. This function is called during ``update_all()``.

* **evolve_disk_init()**

Is used before the planet is put into the disk.

* **calc_timestep()**

calculates the timestep.

* **print_params()**

does the outputing and the occasional printing.
