from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("chill_lib",
    sources = ["chill_lib.pyx","chiller.cpp", "vec3/vec3.cpp"],
    include_dirs=["vec3"],
    libraries = ["boost_math_c99"],
    library_dirs = ["/usr/lib/x86_64-linux-gnu","/home/henriasv/ovito-2.6.2-x86_64/lib/python3.4/config-3.4m"],
    language="c++",
    extra_compile_args=["-std=c++11"]
    )
]

setup(
    ext_modules = cythonize(extensions)
)
