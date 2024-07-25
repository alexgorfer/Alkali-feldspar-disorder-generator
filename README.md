# Alkali feldspar disorder generator
This standalone Python script constructs alkali feldspar supercells with precisely defined Al-Si and Na-K disorder to use in atomistic simulations. It depends on the ASE package and can optionally utilize the sqsgenerator package. This script was introduced in the manuscript:

A. Gorfer, D. Heuser, R. Abart, and C. Dellago, <br>
"Thermodynamics of alkali feldspar solid solutions with varying <br>
Al-Si order: Atomistic simulations using a neural network potential", <br>
[arXiv:2407.17452](https://doi.org/10.48550/arXiv.2407.17452) [cond-mat.mtrl-sci] (2024).

## Input settings
```
--supercell_dim   # The multiplicity of the supercell in a b c (e.g. 3 3 3) 
--T1O_frac        # Fraction of aluminum at T1O-sites (0.0 to 1.0)
--T1M_frac        # Fraction of aluminum at T1M-sites (0.0 to 1.0)
--T2O_frac        # Fraction of aluminum at T2O-sites (0.0 to 1.0)
--T2M_frac        # Fraction of aluminum at T2M-sites (0.0 to 1.0)
                  # Note: T1O_f + T1M_f + T2O_f + T2M_f must equal 1.0
--Na_fraction     # Fraction of sodium in the system  (0.0 to 1.0)
--sqsGen_NaK      # Use sqsgenerator to distribute sodium/potassium (optional)
```
## Examples
Let us recreate the alkali-feldspar archetypes used in the manuscript above. Note that using the original manuscript's supercell dimensions (8 6 8) with --sqsGen_NaK requires substantial memory resources. We will use a smaller supercell here:
### Al ordered
```
python alkali_feldspar_disorder_generator.py \
    --supercell_dim 3 3 3 \
    --T1O_frac 1.0 \
    --T1M_frac 0.0 \
    --T2O_frac 0.0 \
    --T2M_frac 0.0 \
    --Na_fraction 0.5 \
    --sqsGen_NaK
```
### Al T1-disordered
```
python alkali_feldspar_disorder_generator.py \
    --supercell_dim 3 3 3 \
    --T1O_frac 0.5 \
    --T1M_frac 0.5 \
    --T2O_frac 0.0 \
    --T2M_frac 0.0 \
    --Na_fraction 0.5 \
    --sqsGen_NaK
```
### Al disordered
```
python alkali_feldspar_disorder_generator.py \
    --supercell_dim 3 3 3 \
    --T1O_frac 0.25 \
    --T1M_frac 0.25 \
    --T2O_frac 0.25 \
    --T2M_frac 0.25 \
    --Na_fraction 0.5 \
    --sqsGen_NaK
```
