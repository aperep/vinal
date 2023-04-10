This project provides a realization of [Vinberg's algorithm](https://en.wikipedia.org/wiki/Vinberg%27s_algorithm) on [SageMath](http://www.sagemath.org/).

This is a work-in-progress.

It uses [CoxIter](https://github.com/rgugliel/CoxIter) by Rafael Guglielmetti for checking the completeness of the found set of roots.

Created and maintained by [Nikolay Bogachev](https://github.com/nvbogachev) and [Alexander Perepechko](https://github.com/aperep).

# Installation
1. Install [SageMath](http://doc.sagemath.org/html/en/installation/)
2. Add necessary python modules, e.g. `line_profiler`:
```
sage -python -m easy_install line_profiler
```
3. Compile [CoxIter](https://rgugliel.github.io/CoxIter) into a binary `./CoxIter/coxiter`.

# Usage
1. Start the SageMath interactive shell `sage`.
2. Attach the project
```
sage: attach('vinal.sage')
```
3. Create the Vinberg Algorithm instance from a square integer matrix M: 
```
sage: a = VinAl(M)
```
4. Run the roots search
```
sage: a.FindRoots()
```
