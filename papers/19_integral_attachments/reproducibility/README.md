# Exact calculations

Run the following from the paper directory:

    python3 reproducibility/run_all.py

The wrapper also works from another directory. All calculations use exact
rational arithmetic and the Python standard library.

- check_blowup_model.py constructs the diagonal-blow-up algebra, integral middle
  pairing, determinant, Chern pairings and exceptional-plane coordinates.
- check_universal_indices.py calculates four Chern classes for 24 topological
  K-basis inputs and eight further virtual inputs, with integral coordinates and
  plane binomial restrictions.
- check_wall_action.py checks every degree-two, degree-four and degree-six
  generator image, the involution, normalized-index transforms and mod-eight action.
- check_all_attachments.py reconstructs the printed Gram pairing separately,
  checks the six attachment systems and the inverse-transpose boundary maps.

The four JSON outputs retain the coordinates and checks. These programs
supplement the printed proofs of the compact geometric maps, integral topology
and analytic applicability; finite arithmetic does not establish those theorems.
