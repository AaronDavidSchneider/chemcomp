config.yaml
^^^^^^^^^^^
.. Note:: The code is build to use cgs units! Either use cgs values or set unit conversion with `scale` in `chemcomp/helpers/import_config.py`

Example of unit conversion
""""""""""""""""""""""""""

In ``config.yaml``:

.. code-block :: yaml

   config_disk:
     M_STAR:   1.0
     ALPHA:    5.0e-4


+------------+-----------------+-------------------------+-----------------------------------+
|Variable    | Unit            | ``import_config.py``    | ``__init__()`` of module ``Disk`` |
+============+=================+=========================+===================================+
| ``ALPHA``  | unit less value | untouched               | untouched                         |
+------------+-----------------+-------------------------+-----------------------------------+
| ``M_STAR`` | solar masses    | becomes a quantity      | ``M_STAR.cgs.value``              |
|            |                 | (with unit solar masses)|                                   |
|            |                 | by use of ``scale``     |                                   |
+------------+-----------------+-------------------------+-----------------------------------+

The same applies for all modules

The journey of parameters
"""""""""""""""""""""""""

There are generally two ways to bring a value from the configuration file into a module:

1. Use ``kwargs`` (example: ``a_0`` in ``__init__()`` of ``Disk``)
2. Use function arguments in ``__init__()`` (example: ``M_STAR`` in ``__init__()`` of ``Disk``)

The sections
""""""""""""

A config file contains the following sections:

+-----------------------------------+----------------------------------------------------+
| Section                           | Used by ``__init__()`` of                          |
+===================================+====================================================+
| ``config_disk``                   | ``Disk`` module (specify with parameter ``model``) |
+-----------------------------------+----------------------------------------------------+
| ``config_planet``                 | ``Planet`` (specify with parameter ``model``)      |
+-----------------------------------+----------------------------------------------------+
| ``config_pebble_accretion``       | ``PebbleAccretion``                                |
+-----------------------------------+----------------------------------------------------+
| ``config_planetesimal_accretion`` | ``PlanetesimalAccretion``                          |
+-----------------------------------+----------------------------------------------------+
| ``config_gas_accretion``          | ``GasAccretion``                                   |
+-----------------------------------+----------------------------------------------------+
| ``chemistry``                     | ``Chemistry``                                      |
+-----------------------------------+----------------------------------------------------+
| ``defaults``                      | ``Disk``                                           |
+-----------------------------------+----------------------------------------------------+
| ``output``                        | ``Planet``                                         |
+-----------------------------------+----------------------------------------------------+

List of parameters
""""""""""""""""""

You can change all parameters in the ``config.py`` file. Most parameters are optional and have some default value.
All parameters are set in the respective ``__init__()`` methods of the corresponding object.
Default values are either set in the function header or in a call to ``eval_kwargs()`` which extracts the arguments from the corresponding yaml section.
This is an example of a config file. It should contain all important parameters for the configuration of a DiskLabDisk and a BertPlanet.
If you use/write your own Disk- or Planetclass or simply use a different one you likely need different parameters!
In case of doubt check the ``__init__()`` methods of the respective classes!


.. code-block :: yaml

  # default config file

  # Take care of units here! Otherwise adapt in import_config.scale_units!
  # Physical parameters
  config_disk:
    model:                        DiskLabDisk    # Set the name of the Disk model that you want to use. This Disk model should be imported in the ``__init__.py` in the disk directory!
    M_STAR:                       1.0        # mass of central star in sun masses
    ALPHA:                        5.0e-4     # Alpha viscosity
    ALPHAHEIGHT:                  1.0e-4     # alpha_z vertical stirring parameter
    M0:                           0.128      # Diskmass in solar masses
    R0:                           137        # Disk radius in AU
    DTG_total:                    0.02       # total solid to gas ratio
    DTG_pla:                      0.000      # please use 0, thanks. Planetesimal model needs to be changed
    static:                       False      # static disk means no gas/dust evolution
    evaporation:                  True       # use evaporation and condensation (according to the model of Schneider & Bitsch 2021a)
    static_stokes:                False      # use a static stokes number for the large grain population (set stokes number using STOKES in config_pebble_accretion)
    tau_disk:                     1.0e4      # Disk dispersal time after begin_photoevap (in years)
    begin_photoevap:              3.0        # time before we start the diskdispersal in Myr
    temp_evol:                    False      # Dont evolve the temperature!
    evap_width:                   1.0e-3     # width of the evaporation front in AU (see Schneider & Bitsch 2021a)

  chemistry:
    FeH: 0.0          # iron Fraction relative to hydrogen, log relative to the sun
    SH: 0.0           # like above for sulfur. Note that this is possible for every abundancy
    use_FeH: False    # use FeH proxy fit according to Bitsch, Battestini 2020, not compatible with the complete model
    C_frac: 0.2       # Fraction of C/H in carbon grains (will be reduced from CH4)
    # CH4_frac: 0.0     # Set the fraction of C/H in CH4.
    # CO_frac: 0.45     # Set the fraction of C/H in CO.
    # CO2_frac: 0.1     # Set the fraction of C/H in CO2.

  config_planet:
    model: BertPlanet                 # name of Planet used. Should be imported in the ``__init__.py`` file of the planets folder
    matter_removal: True              # remove accreted matter
    use_heat_torque: True             # Use the heating torque
    use_dynamical_torque: True        # Use the dynamical torque
    migration: True                   # migrate the planet. Can be set to false for insitu formation
    M0_fact: 1                        # factor on the Lambrecht2012 transitionmass for the initial seed mass
    #M_c: 0.0045                      # set random initial core mass
    #M_a: 0.0005                      # set random initial envelope mass
    a_p: 30.0                         # set random initial core mass
    t_0: 0.05                        # time before planets start growing, total time: time_disk_0+t_0
    rho_c: 5.5                     # density of the core (g/cm**3)
    r_in: 0.2                  # position in the disk where planetary growth is stopped
    keep_peb_iso: True         # dont change the pebiso mass once the planet has reached it (avoids a second phase of pebble accretion)
    use_pebiso_diffusion: False  # ignore the diffusion part of the pebble isolation mass
    pebiso_start: False          # plant the seed with a mass of the local pebbleisolation mass

  config_pebble_accretion:
    #STOKES:                       0.01     # static Stokes Number, only use with caution, only use in combination with disk config static_stokes
    u_frag:                       5.0       # fragmentation velocity in m/s. If you dont set it the code will use the Drazkowska2017 fragmentation velocity of 1 m/s for dry and 10 m/s for icy pebbles
    epsilon_p:                    0.5       # sticking efficiency
    #H_p_over_H:                   0.1      # static pebble scale hight. Not used by default
    #twoD:                         True     # Not used by default
    #REGIME:                       Hill     # Not used by default

  config_gas_accretion:
    kappa_env:                    0.05      # envelope opacity
    f_machida:                    1         # Machida efficiency
    f_disk_max:                   1.0       # Maximum of the disk accretionrate

  config_planetesimal_accretion:          # Outdated, set efficiency to 0 and DTG_pla to 0!
    R_pla:                        50      # Planetesimal radius in km
    rho_pla:                      1       # density of a single planetesimal (g/cm^3) = 1000 kg/m^3
    stirring:                     1.0e-4  # planetesimal
    efficiency:                   0.00    # Set the efficiency of planetesimal formation
    pla_method:                   MPIA    # can be 'MPIA' or 'Drazkowska' in order to use the formation prescription of Lenz2019 or Drazkowska2017

  # modelling parameters
  defaults:
    DEF_R_IN:                     0.1      # inner r boundary (in AU)
    DEF_R_OUT:                    1000     # outer r boundary (in AU)
    DEF_GRIDSIZE:                 500      # radial gridsize
    DEF_LIN_SPACING:              False    # Spacing of radial grid
    DEF_TIMESTEP:                 10       # custom timestep in years (use 10 years to be safe)

  output:
    name:                        Bert     # name of output file. Will be overwritten if you use a job.yaml
    save_disk:                   True     # output the disk or dont save the disk (saving disk quantities is expensive in terms of storage)
    save_interval:               5000     # snapshot interval for output, time in years
    save_disk_interval:          20       # interval (relative to save_interval) at which disk quantities should be snapshoted/saved
    plot_sigma_live:             False    # Some function to do live plots of certain quantities (see corresponding functions in the DataObject class)
    # acc_files:                          # Outdated, shouldnt be used
    #   - "pebble"
    #   - "planetesimal"
    #   - "gas"
