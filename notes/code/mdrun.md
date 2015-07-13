Source code
===========

This note applies to one of the following file
  * programs/mdrun/mdrun.cpp for v5.0 or later
  * kernel/mdrun.c for v4.6 or eariler versions



The `main()` function and its aliases
=================================

* For v4.5 or earlier, this is the actual main function of `mdrun`.
* For v4.6, it is called `cmain()`.
* For v5.0, it is called `gmx_mdrun()`.



Analyses of `main()`
======================
Calls `mdrunner()`.
