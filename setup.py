from setuptools import setup, find_packages

# Numba typically requires an older version of numpy. So put first.
install_requires = [
                        "numba",
                        "pandas",
                        "numpy",
                        "joblib",
                        "matplotlib",
                        "hxform @ git+https://github.com/rweigel/hxform.git#egg=hxform",
                        "swmfio @ git+https://github.com/GaryQ-physics/swmfio.git#egg=swmfio"
                    ]

setup(
    name='magnetopost',
    version='0.1.0',
    author='Gary Quaresima and Bob Weigel',
    author_email='garyquaresima@gmail.com',
    packages=find_packages(),
    description='Post processing of magnetosphere MHD simulation data',
    install_requires=install_requires
)
