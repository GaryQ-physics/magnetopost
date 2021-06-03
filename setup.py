from setuptools import setup, find_packages

install_requires = ["numpy","pandas"]

setup(
    name='magnetopost',
    version='0.0.0',
    author='Gary Quaresima',
    author_email='garyquaresima@gmail.com',
    packages=find_packages(),
    description='post processing multiple space weather simulation models (focused on swmf)',
    install_requires=install_requires
     )
