# H-bonded Amino acid dimers with RDKit
Let's construct hydrogen-bond preorganized amino acid dimers in a physically meaningful linear geometry.
We will use RDKit for molecular construction + NumPy for geometry manipulation.

We will:
* Generate all 20 × 20 natural amino acid dimers
* Build each amino acid as a neutral backbone fragment (NH₂–CH(R)–COOH)
* Arrange them linearly
* Orient them such that:
 One H from the NH₂ of amino acid 1 is positioned at ~2.0 Å from the carbonyl oxygen (C=O) of amino acid 2
* Save as XYZ files
* CPU-only
* Fully deterministic

**Physical Target Geometry**
* Typical hydrogen bond: H···O distance ≈ 1.8–2.2 Å
* N–H···O angle ≈ 160–180°

**We’ll enforce: ==>  Script: aa_dimer_0.py**
* H···O = 2.0 Å
* Linear N–H···O alignment
  
**What This Does**
For each of the 400 dimers:
* Optimizes both amino acids individually
* Identifies:
i) N atom
ii) One H attached to N
ii) Carbonyl oxygen
iii) Moves second amino acid so: H···O = 2.0 Å and N–H···O approximately linear
iv) Saves XYZ
  
**This Does NOT**
* Optimize the dimer
* Enforce perfect angle
* Consider steric clashes
* Consider protonation state

**We add stricter chemical reasonings now  ==>  Script: aa_dimer_1.py**
* N–H···O = 180°
* H···O = 2.0 Å
* Rotational scan around hydrogen bond axis
* Best steric separation

*Next we add CPU parallelisation on top of it to boost up performance and make it production ready.*
**==>  Script: aa_dimer_1_pal.py**



