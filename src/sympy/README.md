This project **VinAl** is a work-in-progress Python3 module that provides a realization of [Vinberg's algorithm](https://en.wikipedia.org/wiki/Vinberg%27s_algorithm) of searching a fundamental polytope of an arythmetic group. The previous Sage version is described in the article [Vinberg's Algorithm for Hyperbolic Lattices](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=mzm&paperid=11889&option_lang=eng).

## Status
Some parts may work not as intended or not work at all. 


## Installation

### Installing without Sage
1. Install [Python3](https://www.python.org) including the package manager `pip` if they are not already present at your system.
2. Install necessary python modules for the repository directory: `pip install -r requirements.txt`.
Main dependencies here are [SymPy](https://www.sympy.org/) for matrix operations and [`pycddlib`](https://pycddlib.readthedocs.io/en/latest/), a Python wrapper for Komei Fukuda’s [cddlib](https://inf.ethz.ch/personal/fukudak/cdd_home/) used for calculating dual cones in constructing the fundamental cone. There is a known [issue](https://github.com/mcmtroffaes/pycddlib/issues/2) with `pycddlib` installation, so you may need to install the C/C++ library [GMP](https://gmplib.org/) either installing the `libgmp3-dev` module with a package manager or [manually](https://www.mersenneforum.org/showthread.php?t=23079).  
3. Compile Rafael Guglielmetti's [CoxIter](https://rgugliel.github.io/CoxIter) into a binary `./CoxIter/build/coxiter` as [instructed](https://rgugliel.github.io/CoxIter/pageInstall.html). It is used to check whether the found set of roots is complete.

### Installing in Sage
1. Install Sagemath of version at least 9.x (which is on Python 3). You may use the official docker image.
2. Install required packages in Sage: `sage --pip install -r requirements.txt`

## Usage example
```
from vinal import VinAl

M = [[0,1,0], [1,0,0], [0,0,1]]
V = VinAl(M)

roots = V.run()
print(roots)

```

## Contributing

Created and maintained by [Nikolay Bogachev](https://github.com/nvbogachev) and [Alexander Perepechko](https://github.com/aperep). Suggestions and collaboration are appreciated.

