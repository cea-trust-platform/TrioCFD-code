
# TrioCFD

**TrioCFD** (previously named "Trio_U") is the Computational Fluid Dynamics (CFD) code
based on the TRUST platform ("TRUST" with Front_Tracking, Radiation, ALE for fluid/structure interactions and Turbulence LES & RANS models).

This software is OpenSource (BSD license).



# **How to install TrioCFD-1.8.4 version ?**

### If TRUST-1.8.4 is not already installed, [please follow TRUST install instructions](https://github.com/cea-trust-platform/trust-code#readme).

### Once TRUST installed, install TrioCFD-1.8.4 with:

    $> git clone https://github.com/cea-trust-platform/TrioCFD-code.git TrioCFD-1.8.4
    $> cd TrioCFD-1.8.4
    $> source PathToTRUST-1.8.4/env_TRUST.sh
    $> baltik_build_configure -execute
    $> make optim debug

# **How to install TrioCFD's development version ?**
**for developers and those interested in new features only.**

**Warning: "next" branch may not compile or some tests fail if important developments merged**

### If TRUST-next is not already installed, [please follow TRUST install instructions](https://github.com/cea-trust-platform/trust-code/tree/next#readme).

    $> git clone https://github.com/cea-trust-platform/TrioCFD-code.git TrioCFD-next
    $> cd TrioCFD-next
    $> source PathToTRUST-next/env_TRUST.sh
    $> baltik_build_configure -execute
    $> make optim debug

# **How to start ?**

    $> source ./env_TrioCFD.sh

To check:

    $> make check_all_optim


# **How to run a TrioCFD preinstalled version**

If you are a CEA worker, it is possible to use a TMA preinstalled version of TrioCFD. Here is the list of the machines and the paths that can be used to source the TRUST environment

- for a local machine at CEA Saclay

      $> source /home/triou/env_TrioCFD_1.8.4.sh

      or

      $> source /home/trust_trio-public/env_TrioCFD-1.8.4.sh

- for CEA Saclay cluster (orcus):

      $> source /home/trust_trio/env_TrioCFD-1.8.4.sh

- for CCRT-TGCC supercomputers (cobalt & irene-ccrt):

      $> source /ccc/cont002/home/den/triou/env_TrioCFD-1.8.4.sh

- for CCRT-TGCC supercomputers (irene-amd int64):

      $> source /ccc/cont002/home/den/triou/env_TrioCFD-1.8.4-int64.sh

- for CINES supercomputer (occigen):

      $> source /panfs/panasas/softs/applications/trio_u/TrioCFD/env_TrioCFD-1.8.4.sh

