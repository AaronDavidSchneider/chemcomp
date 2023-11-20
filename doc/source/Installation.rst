Installation
------------

Preperation
^^^^^^^^^^^

1. Get access to the code (contact `Aaron Schneider <mailto:aarondavid.schnieder@kuleuven.be>`_ for details)
2. Create a github account
3. crate an ssh-key and add it to your github account (see `here <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_)
4. Install `anaconda <https://www.anaconda.com/products/individual>`_
5. Familiarise with python

Installation
^^^^^^^^^^^^

Clone the project

.. code-block:: bash

   git clone git@github.com:chemcomppla/chemcomp.git

move into directory:

.. code-block:: bash

   cd chemcomp


Install `conda <https://www.anaconda.com/products/individual>`_. virtual environment:

.. code-block:: bash

   conda create -n chemcomp astropy numpy scipy numba matplotlib pyyaml PyTables h5py
   conda activate chemcomp


Install ``chemcomp``:

.. code-block:: bash

   pip install -e .
