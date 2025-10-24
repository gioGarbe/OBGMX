# OBGMX
A generator of GROMACS topologies using Openbabel's UFF force field.

OBGMX is a set of patches to [openbabel-2.3.2](https://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/openbabel-2.3.2.tar.gz/download) (`ob-gg-2.3.2.patch`), which enable dealing with periodic systems, and a program (`obgmx.cpp`) that uses openbabel's UFF implementation to generate [GROMACS](https://www.gromacs.org/) topologies.

The OBGMX project began with the goal of enabling molecular-dynamics simulations of flexible
metal-organic frameworks. Since then, my research interest has drifted away from molecular
simulations so **I can no longer maintain the code on a regular basis**.

Nevertheless, I felt that I should somehow celebrate the [2025 Nobel prize in
Chemistry](https://www.nobelprize.org/prizes/chemistry/2025/summary/), so I enlisted the help of an LLM
and came up with a version of OBGMX that compiles on modern distributions, such as Ubuntu 24.04.3
LTS or macOS Tahoe.

The shell script `install_obgmx.sh` takes care of:

* download openbabel-2.3.2 (using `curl`)
* apply the patch (`ob-2.3.2-gg-2025.patch`)
* build the patched library (needs `cmake` and a C++ compiler)
* install it (might require administrative privileges, obtained via `sudo`)
* build the `obgmx` executable

**NOTE:** this procedure will create a lot of warnings, but it should complete just fine.

Running the `obgmx` executable without any arguments prints a list of available options. Default
values should work just fine, but I refer to the [original
paper](https://doi.org/10.1002/jcc.23049) for the meaning of the various switches.  
Notice that this release adds the possibility of computing _partial charges_ using various models.

In principle, OBGMX can generate topologies for any molecular format supported by
openbabel-2.3.2. In practice, I usually used the `xyz` format for molecules and the `pdb` format for
periodic systems.  
OBGMX produces and .itp file and a .top file, almost ready to be used with GROMACS.

---

If you find OBGMX useful in your research, please [cite the original paper](https://doi.org/10.1002/jcc.23049).

---

## Historical notes

I used to run a web-based interface to produce topologies at `http://software-lisc.fbk.eu/obgmx/`, but unfortunately this service has been discontinued for reasons beyond my control.

OBGMX was originally released as an executable, `obgmx-release`, statically compiled under [Ubuntu
16.04 LTS](https://releases.ubuntu.com/16.04/)/x86-64 that you could use right away provided you
copy the `share` directory under `/usr/local/ob-gg-232`. The executable seemed to work fine in
Ubuntu 22.04.3 LTS/x86-64.  
This version is now distributed in the `Ubuntu 16.04 LTS` directory.

If you want to compile the program yourself, these are the steps

0. Set up an Ubuntu 16.04 LTS environment (a virtual machine will do just fine).
1. Download and uncompress [openbabel-2.3.2](https://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/openbabel-2.3.2.tar.gz/download)
2. Apply `ob-gg-2.3.2.patch`
3. Compile and install OpenBabel
4. Compile `obgmx.cpp` linking the patched openbabel library built in the previous step
5. Enjoy!

<!--
Calling the `obgmx-release` executable without options, will write a list of possible options. Default should work just fine, but for the various switches I refer to the  [original paper](https://doi.org/10.1002/jcc.23049).

In principle, OBGMX can generate topologies for any molecular format supported by openbabel-2.3.2. In practice, I usually used the `xyz` format for molecules and the `pdb` format for periodic systems.
-->

