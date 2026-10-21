"""certcheck: exact, independent verification of vlink certificates.

    from certcheck import check_certificate
    ok, message = check_certificate(cert)

Trusted base: vlink/field.py (exact multiquadratic-field arithmetic) and
Regina's combinatorics (isosig decoding, vertex links, gluings). Nothing else
from vlink is used: bead validity, junction angles and lengths are checked
with different methods from vlink's search (see check.py).
"""

from .check import check_certificate  # noqa: F401
