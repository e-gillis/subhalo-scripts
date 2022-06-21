Package for postprocessing subhalos. Imports a class that can be initialized to
save load and plot subhalo characteristics.


# Package Structure Adapted from:
# PyPackageTemplate
### Fergus Horrobin


Use this template to create a basic python package.

## Steps

1. Change the package metadata in `setup.py`. If more complicated behaviour is needed, additional stuff can be added to `setup.py` per setuptools documentation.

2. Rename the parent directory and package directory (here named `PyPackageTemplate`) to the desired name of your package (ie you will write `import PackageName` where `PackageName` is what you choose here).

3. Edit this README such that it pertains to your package.

4. Add your scripts in the subfolder for your package.

5. Import the functions you want to be able to use from your scripts in `package_folder/__init__.py`.

    For example, if you create a file called `test.py` and in write a function `f(x)` which you want to export then in `__init__.py` write:

    ```from test import f```

    Then you can import your pakcage as `import PackageName` and do `PackageName.f(2)` or you could do `from PackageName import f` and use `f(2)`.

6. If you want to install your package for local development use (such that changes to your package will become available without re-instralling), run the following under the Python environment you want to install it in while in this directory:

    ```pip install -e .```

    Here the `-e` refers to installing in development mode and `.` is telling pip to install using the `setup.py` contained here.
