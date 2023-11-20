Chemistry Module
^^^^^^^^^^^^^^^^

Calculates the chemical composition of the disk (gas and solids).
There are two different ways to calculate the composition:

1. static: ``get_composition()``
2. dynamically: ``get_gas_composition()`` and ``get_solid_composition()``

Static means that the composition of the disk is computed based on the disk temperature. Dynamically means that the composition is computed based on the knowledge of the molecular composition (used if evaporation & condensation is used). 

Compositional vectors (e.g. ``sigma_g_components``) are numpy arrays with dimensions ``N x 2 x M``, where ``N`` is the number of gridcells and ``M`` is the number of included molecular species (19 in case of SB20, 10 in case of BB2020).
Compositional vectors can be accessed in the following way:

.. code-block :: python

   sigma_g_components[i, type, idx]
   # i: Number of the gridcell
   # type: 0 (elements), 1 (molecules)
   # idx: idx of the element / molecule

This is best explained in the function ``calc_composition()``.

Important Functions
"""""""""""""""""""
* **calc_composition()**

takes molecular mass fractions and enumerates the mass fractions of elements based on their molecular composition.

* **get_composition()**

wrapper for the calculation of the static disk chemistry.

* **get_gas_composition()**

wrapper for the calculation of the dynamic gas disk composition based on the molecular composition.

* **get_solid_composition()**

wrapper for the calculation of the dynamic solid disk composition based on the molecular composition.

* **abundances()**

calculates the mass fractions for the static composition.

* **get_position_of_ice()**

calculates the position of icelines.


important attributes
""""""""""""""""""""

+------------+--------+-----------------------------------------------------------------------------------------------------------+
| attribute  | format | meaning                                                                                                   |
+============+========+===========================================================================================================+
| mu         | float  | mean molecular mass                                                                                       |
+------------+--------+-----------------------------------------------------------------------------------------------------------+
| CV         | float  | heat capacity                                                                                             |
+------------+--------+-----------------------------------------------------------------------------------------------------------+
| gamma      | float  | adiabatic index                                                                                           |
+------------+--------+-----------------------------------------------------------------------------------------------------------+
| mask_array | 2xM    | 1 if active dimension, 0 if unactive. This is needed since there are more molecular species than elements |
+------------+--------+-----------------------------------------------------------------------------------------------------------+





