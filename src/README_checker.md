# Exact crescent checker (`src/checker.py`)

`python src/checker.py file.json [--json] [--quiet]`. The exit code is 0 when every configuration in the file is a crescent configuration, and 1 otherwise.
The file holds one spec object or a list of them. Each spec is `{"name": optional, "field": F, "points": [[x, y], ...]}`, and every coordinate is a string.

| `field` | coordinates | example |
|---|---|---|
| `"Q"` | rationals, written `"p/q"` | `["3/5", "-4"]` |
| `{"type":"sqrt","d":d}` | sympy expressions over Q in `sqrt(d)` (or in the symbol `s` = sqrt(d)). `d` must be squarefree and not 0 or 1 | `["5*sqrt(3)/2", "s/2+1"]` |
| `{"type":"minpoly","minpoly":"a**3-2","var":"a"}` | polynomials over Q in `a`, a root of the given polynomial, which must be irreducible over Q | `["1", "a**2/3"]` |
| `{"type":"algebraic"}` | any closed-form real algebraic sympy expressions (radicals, `CRootOf`, ...) | `["1", "sqrt(2)+3**(1/3)"]` |

In Python, `check(points) -> CheckResult` returns `ok`, `reasons` (the duplicate points, the collinear triples, the concyclic quadruples and the bad class sizes), `multiplicities` and `classes` (for each class: the squared distance, its count and its pairs). You can also call `check_spec(dict)` or `load_points(dict)`.

## How equality is decided exactly
- Every coordinate lies in one number field K = Q[t]/(m), where m is irreducible over Q. FLINT's `fmpq_poly.factor()` verifies the irreducibility when K is built. The rationals use m = t, and Q(sqrt d) uses m = t^2 - d.
- An element is stored as its reduced remainder mod m, a FLINT `fmpq_poly` of degree < deg m with exact rational coefficients. K is a field, so 1, t, ..., t^(deg m - 1) is a Q-basis. That makes the representative unique, so **x = y holds exactly when the coefficient vectors are identical**. Squared distances are grouped by hashing these coefficient tuples.
- Only +, - and * are used. There is no division, no square root and no floating point.
  - Collinearity is tested by whether the 2x2 determinant is 0.
  - Concyclicity is tested by whether the 4x4 lifted determinant |x y x²+y² 1| is 0. It is computed as a 3x3 determinant after translating the 4th point to the origin, which changes it only by a sign. Quadruples that contain a collinear triple are already rejected by the collinearity check, so they are skipped here.
- A zero test is a statement about K itself, so it does not depend on which real embedding (which root of m) is meant. A Galois conjugate of 0 is 0.
- Parsing:
  - `sqrt` mode: sympy rewrites `sqrt(12)` as `2*sqrt(3)`. The parser then substitutes `sqrt(d) -> s` and requires the result to be a polynomial over Q in `s`. If any other radical remains, parsing fails with an error. It never approximates.
  - `algebraic` mode: sympy's `primitive_element` and `to_number_field` compute a primitive element θ, its minimal polynomial, and each coordinate as a Q-polynomial in θ. These are exact algebraic algorithms (resultants and factorisation over Q). The one numerical step inside sympy is choosing which irreducible factor vanishes at the given expression, and sympy refines that to separation. The parser checks that the generator's minimal polynomial agrees with θ's. After parsing, all arithmetic is the FLINT arithmetic above.
  - When results matter, prefer the `sqrt` or `minpoly` modes, which involve no sympy decision procedure beyond `Poly(..., domain=QQ)`.

## Lattice search (`src/lattice_search.py`)
`python src/lattice_search.py R k... [--all] [--out file]` searches the hexagonal patch of radius R of the triangular lattice, where the lattice point (a, b) has coordinates (a + b/2, b·√3/2). All arithmetic is in integers:
- The squared distance is p²+pq+q².
- The collinearity test is the integer det of the (a, b)-differences, because the Cartesian determinant equals (√3/2)·det.
- The concyclicity test is the integer det |a b a²+ab+b² 1|, because |X Y N 1| = (√3/2)·|a b N 1| after a column operation.

Both integer criteria are cross-checked against `checker.py` in `tests/test_checker.py`. The search backtracks with a class-size feasibility prune: the class counts, sorted in descending order, must satisfy c_i ≤ k-1-(i-1). Solutions are deduplicated up to translation and the 12 lattice symmetries.
