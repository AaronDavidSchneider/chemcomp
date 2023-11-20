import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='chemcomp',
    version='v1.0',
    packages=setuptools.find_packages(),
    include_package_data=True,
    url='https://github.com/aarondavidschneider/chemcomp',
    license='MIT',
    author='Aaron David Schneider',
    author_email='aaron.schneider@nbi.ku.dk',
    description='planet formation model including pebbles and gas accretion',
    long_description=long_description,
    long_description_content_type="text/markdown",
    scripts=['scripts/chemcomp_main', 'scripts/chemcomp_pipeline'],
    install_requires=[
        "astropy>=4.0",
        "scipy",
        "numpy",
        "h5py",
        "tables",
        "python-slugify",
        "matplotlib",
        "pyyaml",
        "parse",
        "adjustText",
    ]
)
