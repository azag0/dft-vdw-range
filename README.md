# Range of electron correlation in density functional theory

### Requirements

All Python code in the repository requires Python 3.6. The following Python packages are required at various stages (all available on PyPI): Numpy, Scipy, Docopt, Pandas, Seaborn.

To calculate raw data from scratch:

-   [FHI-aims](https://aimsclub.fhi-berlin.mpg.de) for DFT calculations
-   The [MBD code](https://github.com/azag0/mbd) for MBD calculations
-   The [D3 code](http://www.thch.uni-bonn.de/tc/index.php?section=downloads&subsection=getd3) for D3 calculations
-   [Quantum Espresso](http://www.quantum-espresso.org) for VV10 calculations

To generate figures and the manuscript:

-   Full LaTeX installation
-   [STIX fonts](http://www.stixfonts.org)

The project uses [Caf](http://github.com/azag0/caf) to manage calculations, which is 