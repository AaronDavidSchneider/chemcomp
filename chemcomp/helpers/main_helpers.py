import os
import zipfile
from datetime import datetime

import numpy as np
from slugify import slugify
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
from yaml import load, dump
import h5py

import multiprocessing

import chemcomp.disks
import chemcomp.planets
from chemcomp.accretion_models import Accretion
from chemcomp.accretion_models import GasAccretion as GasAcc
from chemcomp.accretion_models import PebbleAccretion as PebAcc
from chemcomp.helpers import import_config, CONFIG_PATH, OUTPUT_PATH, ZIP_OUTPUT_PATH


def run_pipeline(job, delete, continue_job, config):
    files, names = create_pipeline(job, config)

    files_filtered = [f for f in files if file_checks_out(f, continue_job)]

    with multiprocessing.Pool() as pool:
        for item in files_filtered:
            pool.apply_async(run_setup, args=[item])
        pool.close()
        pool.join()

    zip_all(names, delete=delete)


def file_checks_out(config_file, continue_job):
    if not continue_job:
        return True

    output = import_config(config_file).get("output",{})
    name = output.get("name", "planet")
    filename = output.get("directory", os.path.join(OUTPUT_PATH, "{}.h5".format(name)))

    if os.path.isfile(filename):
        try:
            with h5py.File(filename, "r") as f:
                attr = dict(f["/"].attrs)
                if bool(attr.get("finished", 0)):
                    # simulation did finish, we do not need to run the simulation!
                    return False
                else:
                    return True
        except Exception:
            return True
    else:
        # file does not exist, we need to run the simulation!
        return True


def create_pipeline(job, config):
    template = config
    name = job.get("output_name", "")
    save_disk = job.get("save_disk", False)
    save_interval = job.get("save_interval", None)
    double_disks = job.get("double_disks", False)

    parameter_list = []
    parameter_info = []
    for i in range(1, 10):
        data = job.get(f"par_{i}", None)
        if data is None:
            break
        arange = data.get("arange", None)
        logspace = data.get("logspace", None)
        linspace = data.get("linspace", None)
        vals = data.get("vals", None)

        if arange is not None:
            parameter_list.append(np.arange(*arange))
        elif logspace is not None:
            parameter_list.append(np.logspace(*logspace))
        elif linspace is not None:
            parameter_list.append(np.linspace(*linspace))
        elif vals is not None:
            parameter_list.append(np.array(vals))
        else:
            print("please use vals, arange, linspace or logspace")

        parameter_info.append(job.get(f"par_{i}"))

    parameter_mesh = np.meshgrid(*parameter_list)
    values = np.array([mesh.flatten() for mesh in parameter_mesh])

    if double_disks:
        files, names = double_disks_wrapper(template, values, parameter_info, name=name,
                                            save_disk=save_disk, save_interval=save_interval)
    else:
        files = create_yaml_test_range(template, values, parameter_info, name=name, save_disk=save_disk,
                                       save_interval=save_interval)
        names = [name]

    return files, names


def run_setup(config_file):
    """load config and start modelling"""
    if not os.path.isdir(OUTPUT_PATH):
        print(f'creating Folder {OUTPUT_PATH}')
        os.mkdir(OUTPUT_PATH)

    # Load config
    config = import_config(config_file)
    defaults = config.get("defaults", {})
    config_disk = config.get("config_disk", {})
    config_pebble_accretion = config.get("config_pebble_accretion", {})
    config_gas_accretion = config.get("config_gas_accretion", {})
    config_planet = config.get("config_planet", {})
    output = config.get("output", {})
    chemistry_conf = config.get("chemistry", {})
    # initialize Disk and planets
    chemistry = chemcomp.disks.Chemistry(chemistry_conf)

    disk_model = config_disk.get("model", None)
    planet_model = config_planet.get("model", None)
    if disk_model is None:
        disk_model = "DiskLabDisk"
        print("No diskmodel specified, using DiskLabDisk")
    if planet_model is None:
        planet_model = "BertPlanet"
        print("No planetmodel specified, using BertPlanet")

    disk = getattr(chemcomp.disks, disk_model)(defaults, chemistry=chemistry, **config_disk)

    accretion_models = get_acc_models(
        config_pebble_accretion,
        config_gas_accretion,
    )
    planet = getattr(chemcomp.planets, planet_model)(
        output=output, disk=disk, accretion_models=accretion_models, **config_planet
    )
    print("start with planet {}".format(output.get("name")))
    planet.grow_mass()
    print("finished planet {}".format(output.get("name")))

    del accretion_models
    del chemistry
    del disk
    del planet


def zip_all(names=[], delete=False):
    if not os.path.isdir(ZIP_OUTPUT_PATH):
        print(f'creating Folder {ZIP_OUTPUT_PATH}')
        os.mkdir(ZIP_OUTPUT_PATH)

    zipf = zipfile.ZipFile(
        os.path.join(ZIP_OUTPUT_PATH,
                     "{}_{}.zip".format(os.path.commonprefix(names), datetime.now().strftime("%Y%m%d-%H%M%S"))),
        "w",
        zipfile.ZIP_DEFLATED,
    )
    flist = []

    for name in names:
        with open(os.path.join(OUTPUT_PATH, "slugs_{}.dat".format(name))) as f:
            out = f.read().split("\n")
        for o in out:
            flist.append(os.path.join(OUTPUT_PATH, "{}.h5".format(o)))
        flist.append(os.path.join(OUTPUT_PATH, "plot_names_{}.dat".format(name)))
        flist.append(os.path.join(OUTPUT_PATH, "slugs_{}.dat".format(name)))

    flist = list(dict.fromkeys(flist))
    for f in flist:
        common_path = os.path.commonpath([f, ZIP_OUTPUT_PATH])
        zipf.write(f, os.path.relpath(f, common_path))

    zipf.close()

    if delete:
        for f in flist:
            os.remove(f)

    return


def get_acc_models(conf_pebble, conf_gas):
    def is_defined(conf, acc):
        if conf == {}:
            return Accretion(conf)
        else:
            return acc(conf)

    return [
        is_defined(conf_pebble, PebAcc),
        is_defined(conf_gas, GasAcc),
    ]


def write_yaml(template, parameter_info, values, template_name, name="", save_disk=False, save_interval=None):
    if isinstance(template, str):
        with open(os.path.join(CONFIG_PATH, "{}.yaml".format(template))) as f:
            config_template = load(f, Loader=Loader)
    else:
        config_template = template.copy()

    params = [p.get("name", None) for p in parameter_info]

    for i in range(len(parameter_info)):
        section = parameter_info[i].get("section", None)
        sub_section = parameter_info[i].get("sub_section", None)
        param = parameter_info[i].get("name", None)

        config = config_template[section]
        if sub_section is not None:
            config = config[sub_section]

        config[param] = values[i]

    slug = slugify("{}{}".format(name, tuple(values)).replace("-", "neg"))
    plot_name = ", ".join(
        ["{}:{}".format(params[i], values[i]) for i in range(len(values))]
    )
    new_file_name = os.path.join(CONFIG_PATH, "p_{}/{}.yaml".format(name, slug))

    config_template["output"]["name"] = slug
    config_template["output"]["save_disk"] = save_disk
    if save_interval is not None:
        config_template["output"]["save_interval"] = save_interval

    with open(new_file_name, "w") as f:
        dump(config_template, f)  # , Dumper=Dumper)

    return new_file_name, slug, plot_name


def create_yaml_test_range(template, values, parameter_info, name="", save_disk=False, save_interval=None):
    files, slugs, plot_names = [], [], []

    if not isinstance(template, str):
        template_name = template["__FILE_NAME__"]
    else:
        template_name = template

    if not os.path.isdir(CONFIG_PATH):
        print(f'creating Folder {CONFIG_PATH}')
        os.mkdir(CONFIG_PATH)

    path = os.path.join(CONFIG_PATH, "p_{}".format(name))

    create_config_and_output()

    if not os.path.exists(path):
        os.mkdir(path)

    for i in range(len(values.T)):
        file, slug, plot_name = write_yaml(template, parameter_info, values.T[i],
                                           template_name=template_name, name=name,
                                           save_disk=save_disk, save_interval=save_interval)
        files.append(file)
        slugs.append(slug)
        plot_names.append(plot_name)

    with open(os.path.join(OUTPUT_PATH, "slugs_{}.dat".format(name)), "w") as slugfile:
        slugfile.write("\n".join(slugs))

    with open(
            os.path.join(OUTPUT_PATH, "plot_names_{}.dat".format(name)), "w"
    ) as plot_name_file:
        plot_name_file.write("\n".join(plot_names))

    return files


def double_disks_wrapper(template, values, parameter_info, name, save_disk, save_interval):
    """
    change template and create four different disks: With Pla, With evap, With Pla+evap, plain
    Note: currently only two disks. Thank you referee!
    """
    files = []

    with open(os.path.join(CONFIG_PATH, "{}.yaml".format(template))) as f:
        config_template = load(f, Loader=Loader)
        config_template["__FILE_NAME__"] = template

    # None:
    template_None = config_template.copy()
    template_None["config_disk"]["DTG_pla"] = 0.0
    template_None["config_disk"]["evaporation"] = False
    files.append(create_yaml_test_range(template_None, values, parameter_info, name, save_disk, save_interval))

    # evap:
    template_evap = config_template.copy()
    template_evap["config_disk"]["DTG_pla"] = 0.0
    template_evap["config_disk"]["evaporation"] = True
    name_evap = f"{name}_evap"
    files.append(create_yaml_test_range(template_evap, values, parameter_info, name_evap, save_disk, save_interval))

    return [item for sublist in files for item in sublist], [name_evap, name]


def create_config_and_output():
    if not os.path.exists(CONFIG_PATH):
        os.mkdir(CONFIG_PATH)
    if not os.path.exists(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)
    if not os.path.exists(ZIP_OUTPUT_PATH):
        os.mkdir(ZIP_OUTPUT_PATH)
