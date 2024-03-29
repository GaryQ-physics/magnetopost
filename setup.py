from setuptools import setup, find_packages

install_requires = ["numpy","numba","joblib","tqdm"]
# also in git repos swmf_file_reader, hxform

setup(
    name='magnetopost',
    version='0.0.1.1',
    author='Gary Quaresima',
    author_email='garyquaresima@gmail.com',
    packages=find_packages(),
    description='post processing multiple space weather simulation models (focused on swmf)',
    install_requires=install_requires
     )
