# Overview

A set of post--processing scripts written in support of the paper [Blake et al., 2021, Recreating the Horizontal Magnetic Field at Colaba During the Carrington Event With Geospace Simulations](https://doi.org/10.1029/2020SW002585) and a draft follow-up paper.

Given directory containing output files from a SWMF run<sup>\*</sup>, the magnetic field at given locations in the magnetosphere and ionosphere are computed by

1. computing the surface integral over the boundary between the near--Earth and far--Earth region (denoted by $\oint_{\mathcal{I}}$),
2. adding the Biot--Savart integral ($\int_{\mathcal{M}}\text{Biot--Savart}$) and the Coulomb integral ($\int_{\mathcal{M}}\text{Coulomb}$) over the far--Earth region, $\mathcal{M}$, to the a surface integral over the outer simulation surface ($\oint_{\mathcal{O}}$), and
3. computing the Biot--Savart integral, $\int_{\mathcal{M}}\text{Biot--Savart}$, only,

The resulting field from these methods is then added to the Biot--Savart integral of the currents in near--Earth region to give the total $\mathbf{B}$ on Earth's surface.

<sup>\*</sup> Native SWMF output files or CCMC .cdf; support for output from other global magnetosphere models may be added.

# Use

Requires Python 3. Last tested on Python 3.8.8.

See the [runs](https://github.com/GaryQ-physics/magnetopost/tree/main/runs) directory for sample configuration scripts and output figures.

## User

```
pip install 'git+https://github.com/GaryQ-physics/magnetopost.git' --upgrade
```

## Developer

```
git clone https://github.com/GaryQ-physics/magnetopost.git
cd magnetopost
pip install --editable .
```


# Acknowledgments

This work was in part supported by NASA Grant 80NSSC20K0589 "Physics-based modeling of the magnetosphere-ionosphere system under Carrington-scale solar driving: response modes, missing physics and uncertainty estimates", PI: Antti Pulkkinen and the subaward "Ground Magnetic Field Perturbations under Extreme Space Weather Conditions", PI: R.S. Weigel.

