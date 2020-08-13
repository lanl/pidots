PIDOTS: Parallel Integral Discrete Ordinates Transport Solver

Updated 8/13/2020 by zerr

PIDOTS is a software application for studying an alternative iterative and
parallel solution technique for linear, deterministic neutral-particle radiation
transport. The code solves a very generic, simplified—i.e., dumb—form of the
system of relevant system of equations. PIDOTS is less capable than the open
source proxy application SNAP. Given limitations introduced by the equations as
well as limitations in data and meshing, PIDOTS is incapable of performing real
physics simulations or addressing actual physics questions for transport.
Instead of using transport mesh sweeps that we see in many production codes (and
SNAP), PIDOTS uses a Jacobi or Jacobi-like iterative scheme that updates the
flux information between adjacent, decoupled subdomains. By completely
decoupling subdomains, PIDOTS permits massively large-scale decomposition (i.e.,
tens of thousands of processes or more) of the spatial mesh, making it an
interesting alternative to mesh sweeps, worth continued studying and
development.

PIDOTS was developed as a research product at Penn State, North Carolina State
University, and Los Alamos National Laboratory, based on original research
performed at Oak Ridge National Laboratory.


Principal Author:

Joe Zerr, Los Alamos National Laboratory, GitHub Username: zerr


Contributors:

Raffi Yessayan, North Carolina State University

Yousry Azmy, North Carolina State University

