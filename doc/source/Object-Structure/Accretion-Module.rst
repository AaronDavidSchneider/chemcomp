Accretion Module
^^^^^^^^^^^^^^^^

Used in the flavours of pebble, planetesimal and gas accretion. Calculates the accretion rates.

Important Functions
"""""""""""""""""""

* **calc_m_dot()**

is called by the planet at each timestep to update the recalculate the accretion rates.
This method is overloaded by the child classes (pebbles, planetesimals, gas).

* **_init_params()**

is called by ``init_accretion()`` which is called at the final step during the planets initialisation. This function is especially useful to communicate parameter between the modules (see e.g. `PebbleAccretion`).

* **remove_mass_from_flux()**

is called at the final stage of the timestep. The removal of mass from the disk should be handled in this method.

Important Attributes
""""""""""""""""""""

+------------------+--------+-------------------------------------------------------------------------------+
| attribute        | format | meaning                                                                       |
+==================+========+===============================================================================+
| ``m_dot``        | float  | total sum of the accretion rate                                               |
+------------------+--------+-------------------------------------------------------------------------------+
| ``m_dot_a``      | float  | accretion rate onto the atmosphere                                            |
+------------------+--------+-------------------------------------------------------------------------------+
| ``m_dot_c``      | float  | accretion rate onto the core                                                  |
+------------------+--------+-------------------------------------------------------------------------------+
| ``m_dot_a_chem`` | 2xM    | accretion rate onto the atmosphere disaggregated into the chemical compositon |
+------------------+--------+-------------------------------------------------------------------------------+
| ``m_dot_c_chem`` | 2xM    | accretion rate onto the core disaggregated into the chemical compositon       |
+------------------+--------+-------------------------------------------------------------------------------+

where M is the number ob molecules in ``Chemistry``.




