# Outdated!

## Welcome to the ramses wiki !

The [ramses][1] code is intended to be a versatile platform to develop
applications using  Adaptive Mesh  Refinement (AMR)  for computational
astrophysics.  The current implementation allows solving the classical
and relativistic Euler equations in presence of self-gravity, magnetic
field  and radiation  field.   The  [ramses][1] code  can  be used  on
massively  parallel  architectures,  if  properly linked  to  the  MPI
library.  It  can also  be used on  single processor  machines without
MPI.  Output  data are generated  using Fortran unformatted  files.  A
suite  of post-processing  routines  is delivered  within the  present
release,  allowing  the user  to  perform  a  simple analysis  of  the
generated output files.

[1]: https://github.com/ramses-organisation/ramses
[2]: https://bitbucket.org/ohahn/music


# [Click here for the user's guide][3]

# [Click here for the automatic tests page][4]

# [Click here for a list of useful tools and external links][6]

## About this wiki

The goal of this wiki is to provide a step-by-step tutorial in
using the [ramses][1] code.  We will first describe a simple example
of its use.  More complex set-up will  be addressed in a second step.
Typical [ramses][1] users can be grouped into 3 categories:

* Beginners: It is possible to execute [ramses][1] using only run time
parameter files. The code is compiled once, and the user only modifies
the  input   file  to   perform  various   simulations.   Cosmological
simulations, for example,  can be performed quite  efficiently in this
way, using the initial conditions  generated by external packages such
as [music][2].

* Intermediate  users: For  more  complex applications,  the user  can
easily modify  a small  set of  Fortran routines  in order  to specify
specific  initial or  boundary conditions.  These routines  are called
__patches__.  The code  should be recompiled each  time these routines
are modified.

* Advanced users:  It is finally  possible to modify the  base scheme,
add new equations, or add new  routines in order to modify the default
[ramses][1]  application.   This  wiki  will describe  some  of  these
advanced features.