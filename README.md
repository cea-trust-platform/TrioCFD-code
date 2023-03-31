
# TrioCFD

**TrioCFD** (previously named "Trio_U") is the Computational Fluid Dynamics (CFD) code
based on the TRUST platform ("TRUST" with Front_Tracking, Radiation, ALE for fluid/structure interactions and Turbulence LES & RANS models).

This software is OpenSource (BSD license).



# **How to install TrioCFD-1.9.1 version ?**

### If TRUST-1.9.1 is not already installed, [please follow TRUST install instructions](https://github.com/cea-trust-platform/trust-code#readme).

### Once TRUST installed, install TrioCFD-1.9.1 using one of these methods:

### **First method**

    $> git clone https://github.com/cea-trust-platform/TrioCFD-code.git TrioCFD-1.9.1
    $> cd TrioCFD-1.9.1
    $> source PathToTRUST-1.9.1/env_TRUST.sh
    $> baltik_build_configure -execute
    $> make optim debug

### **Second method**

    $> wget ftp://ftp.cea.fr/pub/TRUST/TrioCFD/versions/v1.9.1/TrioCFD-1.9.1.tar.gz
    $> tar xzf TrioCFD-1.9.1.tar.gz
    $> mv TrioCFD TrioCFD-1.9.1
    $> cd TrioCFD-1.9.1
    $> source PathToTRUST-1.9.1/env_TRUST.sh
    $> baltik_build_configure -execute
    $> make optim debug

# **How to install TrioCFD's development version ?**
**for developers and those interested in new features only.**

**Warning: "next" branch may not compile or some tests fail if important developments merged**

### If TRUST-next is not already installed, [please follow TRUST install instructions](https://github.com/cea-trust-platform/trust-code/tree/next#readme).

    $> git clone https://github.com/cea-trust-platform/TrioCFD-code.git TrioCFD-next
    $> cd TrioCFD-next
    $> git checkout next
    $> source PathToTRUST-next/env_TRUST.sh
    $> baltik_build_configure -execute
    $> make optim debug

# **How to start ?**

    $> source ./env_TrioCFD.sh

To check:

    $> make check_all_optim

To see documentation:

    $> triocfd -index
