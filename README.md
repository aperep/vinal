This project **VinAl** is a realization of [Vinberg's algorithm](https://en.wikipedia.org/wiki/Vinberg%27s_algorithm) of searching a fundamental polytope of an arythmetic group. The Sage version is described in the article [Vinberg's Algorithm for Hyperbolic Lattices](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=mzm&paperid=11889&option_lang=eng).



## Installation

1. Please refer to the [Readme](src/sage/README.md) of the stable Sagemath version.
2. There is also an experimental sympy version available, see [Readme](src/sympy/README.md).

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

