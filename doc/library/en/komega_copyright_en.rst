Re-distribution of this library
===============================

Contain Komega in your program
------------------------------

:math:`K\omega` library is distributed with the :ref:`lgplicense` (LGPL).
It is summarized as follows:

 * :math:`K\omega` can be freely distributed, modified, copied and pasted,
   in a private program (in the research group, co-workers, etc.).
   
 * For the released program (open-source, free, commercial software etc.):
   
    * When you contain the source-code of :math:`K\omega` (either as is and modified)
      in the distributed source code of your program,
      please distribute your program with LGPL/GPL.
      
    * If you do not include the source-code of :math:`K\omega` (just call it),
      you can freely distribute your program with any licenses.
      
    * If you distribute a binary file which is statically linked to :math:`K\omega` library,
      please use LGPL/GPL. However, if you distribute a binary file which is dynamically linked to
      :math:`K\omega` library (therefore :math:`K\omega` itself is not contained),
      you can freely distribute your binary file with any licenses.

Build Komega without Autoconf
-----------------------------

In this package, :math:`K\omega` is built with Autotools (Autoconf, Automake, Libtool).
If you do not want to use Autotools for your distributed program with :math:`K\omega` source,
you can use the following simple Makefile (please care about TAB).

.. code-block:: makefile

   F90 = gfortran
   FFLAGS = -fopenmp -g -O2 #-D__MPI -D__NO_ZDOT -D__KOMEGA_THREAD
   
   .SUFFIXES :
   .SUFFIXES : .o .F90
   
   OBJS = \
   komega_cg_c.o \
   komega_cg_r.o \
   komega_cocg.o \
   komega_bicg.o \
   komega_math.o \
   komega_vals.o
   
   all:libkomega.a
   
   libkomega.a:$(OBJS)
        ar cr libkomega.a $(OBJS)
   
   .F90.o:
        $(F90) -c $< $(FFLAGS)
   
   clean:
        rm -f *.o *.a *.mod
   
   komega_cg_c.o:komega_math.o
   komega_cg_c.o:komega_vals.o
   komega_cg_r.o:komega_math.o
   komega_cg_r.o:komega_vals.o
   komega_cocg.o:komega_math.o
   komega_cocg.o:komega_vals.o
   komega_bicg.o:komega_math.o
   komega_bicg.o:komega_vals.o
   komega_math.o:komega_vals.o

Preprocessor macros ``__MPI``, ``__NO_ZDOT``, and ``__KOMEGA_THREAD`` correspond to
``--with-mpi=yes``, ``--disable-zdot``, and ``--enable-thread`` of the options of ``configure``, respectively.
      
.. _lgplicense:
      
License
=======

GNU LESSER GENERAL PUBLIC LICENSE Version 3

*¢í 2016- The University of Tokyo. All rights reserved.*

The development of KOmega is supported by ¡ÈProject for advancement of software usability in materials science¡É of ISSP, The University of Tokyo.

(¡ö) We hope that you cite the following paper and repository when you publish the results using K¦Ø.
Papaer: : `¡ÈK¦Ø ¡½ Open-source library for the shifted Krylov subspace method of the form (zI¡ÝH)x=b¡É, Takeo Hoshi, Mitsuaki Kawamura, Kazuyoshi Yoshimi, Yuichi Motoyama, Takahiro Misawa, Youhei Yamaji, Synge Todo, Naoki Kawashima, Tomohiro Sogabe, Computer Physics Communications, Volume 258, January 2021, 107536.<https://www.sciencedirect.com/science/article/pii/S0010465520302551>`_
