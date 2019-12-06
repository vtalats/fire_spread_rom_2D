# Fire Spread ROM (2D)
This project provides a 2D fire spread model suggested by [1].  The full-order solution is
generated using finite differences.  The ROM is implemented using POD, POD/DEIM, and POD/IRKA.  This
project provides extensive test code to compare the results of each of these techniques.

# Running the model
## run_adr.m
This allows you to run the full model, save the matrices, and perform POD and POD/DEIM ROMs using
those matrices.

## adr\_2D\_bd.m
Solves the full model using backward differences.

## build\_adr\_mats2D.m
Builds the matrices used for the backward difference solution. Does not actually solve the system.


_See individual files for additional documentation_

**Alan Lattimer, (c)2015-2019**

## Updates

* December 2019
  * Minor cleanup of the files in the directory.
  * Updated REDAME.

## Other Information

* **twod\_grid\_Q.m** and **twod\_mesh\_Q.m** developed by Alan Lattimer (c)2015-2109 based on code
developed by Jeff Borggaard.
* **\@Msgcls** is a wrapper class for level based logging. Alan Lattimer (c)2013-2019.

