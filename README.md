# OBGMX
A generator of GROMACS topologies using Openababel's UFF force field.

OBGMX is a set of patches to [openbabel-2.3.2](https://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/openbabel-2.3.2.tar.gz/download) (`ob-gg-2.3.2.patch`), which enable dealing with periodic systems, and a program (`obgmx.cpp`) that uses openbabel's UFF implementation to generate [GROMACS](https://www.gromacs.org/) topologies.

I used to run a web-based interface to produce topologies at `http://software-lisc.fbk.eu/obgmx/`, but unfortunately this service has been discontinued for reasons beyond my control.
Additionally, my research interested had drifted away from molecular simulations so **I am not able to maintain this code anymore**.

This release comes with a statically compiled executable, made under [Ubuntu 16.04 LTS](https://releases.ubuntu.com/16.04/)/x86-64 that you can use right away provided you copy the `share` directory under `/usr/local/ob-gg-232`. At the time of this writing, the executable works fine in Ubuntu 22.04.3 LTS/x86-64.

If you want to compile the program yourself, these are the steps

0. Set up an Ubuntu 16.04 LTS environment (a virtual machine will do just fine).
1. Download and uncompress [openbabel-2.3.2](https://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/openbabel-2.3.2.tar.gz/download)
2. Apply `ob-gg-2.3.2.patch`
3. Compile and install OpenBabel
4. Compile `obgmx.cpp` linking the patched openbabel library built in the previous step
5. Enjoy!

Calling the `obgmx` executable without options, will write a list of possible options. Default should work just fine, but for the various switches I refer to the  [original paper](https://doi.org/10.1002/jcc.23049).

In principle, OBGMX can generate topologies for any molecular format supported by openbabel-2.3.2. In practice, I usually used the `xyz` format for molecules and the `pdb` format for periodic systems.

---

If you find OBGMX useful in your research, please [cite the original paper](https://doi.org/10.1002/jcc.23049).

---
