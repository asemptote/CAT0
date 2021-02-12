# Constructible spherical trigonometry in Python

While there exist a range of spherical trigonometry libraries for areas such as
geodesy, these libraries use floating point, or at best rational numbers.
Constructible numbers, on the other hand, have historical and geometric
importance, being those numbers which can be 'constructed' using a ruler and
compass. This implementation thus generalises the exact arithmetic used by
arbitrary-precision rational-number-based implementations.

The mathematical derivations underlying the algorithms in
`ConstrSphPoint.ang(..)` and `ConstrSphPoint.linear(..)` are presented in the
file `constrsphtrig.ipynb`.

## Requirements
Python 2 or 3

## Current state:
- Write unit tests
- Some usage examples
- Upload to PyPI

## Related links:
- https://github.com/bth2008/python_spheric
- https://pypi.org/project/algebraics/#description
- http://abstract.ups.edu/aata/section-constructions.html
    
## Dependencies:
- [`constructible.py`](https://github.com/leovt/constructible): To use a different underlying constructible number object, simply
        change what `constrsphtrig.sqrt` refers to after importing. (If
        `constructible.py` is not available create a file with stub variables
        `sqrt` and `Constructible`.)
