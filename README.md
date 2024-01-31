
# TrioCFD

**TrioCFD** (previously named "Trio_U") is the Computational Fluid Dynamics (CFD) code
based on the TRUST platform.
There are different physical modules such as:
- Turbulence LES & RANS models,
- Front-Tracking,
- Radiation,
- ALE for fluid/structure interactions.

This software is OpenSource (BSD license).



# **How to install TrioCFD-1.9.3 version ?**

### If TRUST-1.9.3 is not already installed, [please follow TRUST install instructions](https://github.com/cea-trust-platform/trust-code#readme).

### Once TRUST installed, install TrioCFD-1.9.3 using one of these methods:

### **First method**
```bash
git clone https://github.com/cea-trust-platform/TrioCFD-code.git TrioCFD-1.9.3
cd TrioCFD-1.9.3
source PathToTRUST-1.9.3/env_TRUST.sh
baltik_build_configure -execute
make optim debug
```

### **Second method**
```bash
wget ftp://ftp.cea.fr/pub/TRUST/TrioCFD/versions/v1.9.3/TrioCFD-1.9.3.tar.gz
tar xzf TrioCFD-1.9.3.tar.gz
mv TrioCFD TrioCFD-1.9.3
cd TrioCFD-1.9.3
source PathToTRUST-1.9.3/env_TRUST.sh
baltik_build_configure -execute
make optim debug
```

# **How to install TrioCFD's development version ?**
**for developers and those interested in new features only.**

**Warning: "next" branch may not compile or some tests fail if important developments merged**

### If TRUST-next is not already installed, [please follow TRUST install instructions](https://github.com/cea-trust-platform/trust-code/tree/next#readme).
```bash
git clone https://github.com/cea-trust-platform/TrioCFD-code.git TrioCFD-next
cd TrioCFD-next
git checkout next
source PathToTRUST-next/env_TRUST.sh
baltik_build_configure -execute
make optim debug
```
# **How to start ?**
```bash
source ./env_TrioCFD.sh
```

To check:
```bash
# All non-regression test cases:
make check_optim or make ctest_optim
make check_debug or make ctest_debug

# A given non-regression test list
make check_optim TESTLIST="./share/testList/{nameOfTheTestList}"

# All validation report (warning it may take many days !):
make check_validation

# A given list of validation reports
make check_validation TESTLIST="./share/testList/{nameOfTheValidationList}
```

To see documentation:
```bash
triocfd -index
```
