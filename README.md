# SolidFEM
A library of FEM algorithms for elasticity

Open the Visual Studio 2019 solution in order to build the library and the application. Make sure you have the "Desktop development with C++" and "Python development" components installed. Use vcpkg to install the assimp library.

In order to build the Python extensions using Visual Studio 2019 please follow the instructions for installing pybind11 using vcpkg from https://pybind11.readthedocs.io/en/stable/installing.html
This will probably force the Python version to 3.8, therefore you will need to use Anaconda or Miniconda until Visual Studio supports this version.
I recommend installing Miniconda and using the "Python 3.8 (64 bit)" environment inside Visual Studio. Also please set the PYTHONPATH env var to the root folder of this Python installation.
Be very careful with recent versions of vcpkg as they use Python 3.9. I advise to roll back before that (e.g., revision e803bf11).
Please use the pytorch channel when installing packages in conda to avoid conflicts, notably the matplotlib and scipy packages.
Find installation instructions on the torch-scatter page: https://github.com/rusty1s/pytorch_scatter

For Intel MKL support please go to the official website and download the SDK. You may have to copy mkl_sequential.dll to x64/Release or add it to a system path.
