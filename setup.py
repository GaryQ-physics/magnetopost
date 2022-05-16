from setuptools import setup, find_packages

install_requires = [
                        "numpy",
                        "pandas",
                        "numba",
                        "joblib",
                        "matplotlib",
                        "hxform @ git+https://github.com/rweigel/hxform.git#egg=hxform",
                        "swmf_file_reader @ git+https://github.com/GaryQ-physics/swmf_file_reader.git#egg=swmf_file_reader"
                    ]

setup(
    name='magnetopost',
    version='0.0.2',
    author='Gary Quaresima and Bob Weigel',
    author_email='garyquaresima@gmail.com',
    packages=find_packages(),
    description='Post processing of magnetosphere MHD simulation data',
    install_requires=install_requires
)
