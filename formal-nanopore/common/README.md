# Formal 50 ns high-salt DNA pore workflows

The four `formal-50ns` systems use the median-Rg representative structure from
the final 10 ns of each formal DNA-only 3 M KCl trajectory. The 5' end is
placed above the cis entrance, while the placement search independently picks
the lowest clash-free Z position and azimuth for each pore.

The original pore boxes did not contain enough cis-side solvent for a folded
40-mer. The formal boxes therefore extend only along Z. The electric field is
rescaled to retain a 1.0 V transmembrane potential.

Preparation is 50,000 minimization steps, 0.5 ns heating, 1 ns strongly
restrained equilibration, 1 ns weakly restrained equilibration, and 2 ns with
DNA free. Production is 5 ns at a 1 fs timestep.
