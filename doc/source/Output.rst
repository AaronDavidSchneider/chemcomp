Output
------

Outputing is handled in ``chemcomp/chemcomp/planets/_planet_class.py`` with the class ``DataObject``.

Every run creates a ``.h5`` file with the following structure:

.. code-block :: bash

   .
   ├── disk              # all the disk quantities
       ├── [quantities]
   ├── planet            # all the planet quantities
       ├── [quantities]
       ├── units
   ├── acc               # accretion related quantities
       ├── [quantities]


The ``print_params()`` inside the Planet class is the function that deals with IO operations. This function is called in the mainloop (``grow_mass()``).

The ``print_params()`` function calls the save function of the ``DataObject`` that belongs to the ``Planet``. All the saving arithmetics (how often we want outputs, etc.) are handled therein.

It is important to understand that the ``DataObject`` object owns a reference to the ``Planet`` (and therefore also to the ``Disk`` object and the accretion objects). The save function of the ``DataObject`` has therefore full access (at all times) to the current disk and planet quantities.

Get familiar with the output
""""""""""""""""""""""""""""

``.h5`` files are self documented files!

#Browse through the content using:
#
#.. code-block :: python
#
#   import h5py
#   file = "yourfilename.h5"
#   with h5py.File(file, "r") as f:
#       print(dict(f["planet"].keys()))
#       print(dict(f["acc"].keys()))
#       print(dict(f["disk"].keys()))
#

Units for planets can be retrieved by

.. code-block :: python

   with h5py.File(file, "r") as f:
     units = dict(f["planet"]["units"].attrs)

.. toctree::
   :maxdepth: 2
   :caption: Plotting Tutorial:

   notebooks/basic_plotting/basic_plotting


Further Examples
""""""""""""""""
   * See `import_data()` in `chemcomp/chemcomp/helper/analysis_helper` for an example of how to load planet data
   * See `chemcomp/analysis/paper_1/disk.py` for an example of how to load disk data
