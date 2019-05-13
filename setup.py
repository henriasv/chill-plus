from setuptools import setup, Extension, find_packages
import numpy 
from Cython.Build import cythonize

chill_lib = Extension("chill_plus", 
    ["chill_lib/chill_lib.pyx","chill_lib/chill_lib.cpp", "vec3/vec3.cpp"], 
    include_dirs=["chill_lib", "vec3", numpy.get_include()],
    libraries = ["boost_math_c99"],
    library_dirs = ["/usr/lib/x86_64-linux-gnu"],
    extra_compile_args=["-std=c++11"])

setup(
    name="chill_plus",
    version="0.1",
    description = "The chill plus algorithm",
    author="Henrik Andersen Sveinsson",
    author_email="henrik.sveinsson@me.com",
    ext_modules=cythonize([chill_lib])
)
