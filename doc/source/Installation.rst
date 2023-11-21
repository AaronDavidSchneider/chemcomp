Installation
------------

Optional Preperationsteps
^^^^^^^^^^^^^^^^^^^^^^^^^

1. Create a github account
2. crate an ssh-key and add it to your github account (see `here <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_)
3. Familiarise with python
4. Install `anaconda <https://www.anaconda.com/products/individual>`_

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

   conda create -n chemcomp -f ci/environment-3.11.yml
   conda activate chemcomp


Install ``chemcomp``:

.. code-block:: bash

   pip install -e .


Installation from pypi
^^^^^^^^^^^^^^^^^^^^^^

Instead of cloning the repository, you can also install ``chemcomp`` from pypi

.. code-block:: bash

   pip install chemcomp
