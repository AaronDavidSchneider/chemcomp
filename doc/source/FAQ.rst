FAQ
---

How do I contribute?
""""""""""""""""""""
If you want to contribute feel free to create a pull request! This way others can profit from your contribution as well!

How do I do a disk only run?
""""""""""""""""""""""""""""
Checkout the ``config/config_paper_1.yaml`` together with ``jobs/paper_1_disk.yaml`` for an example

What happens, if I don't set a paramter
"""""""""""""""""""""""""""""""""""""""
A default value will be used in most cases. Make sure to set all parameters, so that you can be sure about its value!

Where do I find a list of Parameters?
"""""""""""""""""""""""""""""""""""""
Head over to the Configuration -> config.yaml

I want to debug but my runs take too long
"""""""""""""""""""""""""""""""""""""""""
Just use a higher value for the timestep (``DEF_TIMESTEP``), like 500 years. 
But make sure to change it back when you want to do science!