*********************************************************
Generating Functions for Special Multivariate Polynomials
*********************************************************

This Python-module computes generating functions for the number of
some special multivariate polynomials over finite fields. The
accompanying paper [GVZ13]_ also provides proofs and asymptotics with
explicit error terms.

Module
======

multi-def.sage
    run within Sage to obtain symbolic expressions (q) for the number
    of reducible/s-powerful/absolutely irreducible polynomials at
    degree n in r variables over GF(q).

Usage
=====
First, load the module in your local Sage installation::

   $ sage -q
   sage: load('multi-def.sage')

Then, compute generating functions up to the desired degree. For
example, the number of reducible, bivariate monic polynomials at
degree 10 over GF(q) is::

   sage: RR(2, z, 10).coeff(z^10)
   q^56 + 2*q^55 + 2*q^54 + 2*q^53 + 2*q^52 + 2*q^51 + 2*q^50 +
   3*q^49 + 3*q^48 + 2*q^47 - q^45 + q^43 - q^41 - 3/2*q^40 + q^39 +
   5/2*q^38 - 11/2*q^36 - 9*q^35 - 23/2*q^34 - 15*q^33 - 16*q^32 -
   9*q^31 + 19/2*q^30 + 30*q^29 + 93/2*q^28 + 46*q^27 + 31/2*q^26 -
   181/5*q^25 - 147/2*q^24 - 58*q^23 + 2*q^22 + 53*q^21 +
   509/10*q^20 + 15/2*q^19 - 25*q^18 - 41/2*q^17 - 1/2*q^16 +
   15/2*q^15 + 2*q^14 - 5/2*q^13 - 3/2*q^12 + q^11 + 17/10*q^10 +
   1/2*q^9 - 1/2*q^8 - 1/2*q^7 + 3/10*q^5 + 1/10*q^4 - 1/5*q^2 -
   1/10*q

And the number of trivariate, squareful monic polynomials at degree 5
over GF(q) is::

   sage: SS(3, 2, z, 5).coeff(z^5)
   q^55 + q^54 + q^53 + q^52 + q^51 + q^50 + q^49 + q^48 + q^47 +
   q^46 + q^45 + q^44 + q^43 + q^42 + q^41 + q^40 + q^39 + q^38 +
   q^37 + q^36 + q^35 - q^22 - 2*q^21 - 3*q^20 - 3*q^19 - 3*q^18 -
   3*q^17 - 3*q^16 - 3*q^15 - 3*q^14 - 3*q^13 - 3*q^12 - 3*q^11 -
   3*q^10 - 2*q^9 + 3*q^7 + 5*q^6 + 5*q^5 + 3*q^4 + q^3

Finally, the number of absolutely absolutely irreducible, bivariate monic
polynomials at degree 7 over GF(q) is::

   sage: AA(2, z, 7).coeff(z^7)
   q^35 + q^34 + q^33 + q^32 + q^31 + q^30 - q^28 - 2*q^27 - 2*q^26 -
   3*q^25 - 3*q^24 - 3*q^23 - 2*q^22 + 3*q^20 + 8*q^19 + 9*q^18 +
   4*q^17 - 5*q^16 - 10*q^15 - 6*q^14 + 2*q^13 + 5*q^12 + 2*q^11 -
   q^10 - q^9

Requirements
============

This code requires the free mathematical software [Sage]_ which is
available for download at http://www.sagemath.org and as cloud service
at https://cloud.sagemath.org. It has been tested under GNU/Linux with
Sage 6.4.


References
==========

.. [GVZ13] Joachim von zur Gathen, Alfredo Viola & Konstantin
  Ziegler (2013). Counting reducible, powerful, and relatively irreducible
  multivariate polynomials over finite fields. *SIAM Journal on
  Discrete Mathematics* **27**\(2):855–891. URL
  http://dx.doi.org/10.1137/110854680. Also available at
  http://arxiv.org/abs/0912.3312. Extended Abstract in Alejandro
  López-Ortiz (ed.), Proceedings of LATIN 2010, Oaxaca, Mexico, volume
  6034 of Lecture Notes in Computer Science, 243–254 (2010).

.. [Sage] W. A. Stein et al. (2014). Sage Mathematics Software (Version
  6.4). The Sage Development Team. URL http://www.sagemath.org.

Author
======

- Konstantin Ziegler (2013-12-24): initial version

License
=======

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/.
