# Exact representation calculations

Run with Python 3, from this directory or by absolute path:

```sh
python3 ih_llv_test.py
python3 t4_restricted_character.py
```

No external dependencies. The first calculation constructs Hilbert-scheme
characters through degree nine, applies the exterior-square identity and
Freudenthal recursion, and checks the 23-type decomposition and its dimension
138065600. The second recovers the same answer from the corank-one torus and
runs failure controls. Expected transcripts are in `output/`; elapsed times
may differ. All operations affecting the results use exact arithmetic.

The geometric interpretation in complex dimension eighteen is conditional on
(H3), the Hodge-structure comparison printed in Section 10. These calculations
verify finite representation data, not that geometric hypothesis or the
analytic proofs. The lower-dimensional O'Grady comparison is separate.
