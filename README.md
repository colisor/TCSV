# TCSV — Variational Semantic Field Theory

**Power-law spectral structure of synaptic connectivity matrices  
is conserved across cortical areas, layers, and species**

Gesivaldo Santos  
Department of Biological Science, UESB, Bahia, Brazil

## Results summary

| Metric | BBP L2/3 | MICrONS V1 |
|--------|----------|------------|
| N neurons | 8,173 | 50,948 |
| Spectral attractors | 20 | 22 |
| Exponent alpha | -0.450 +/- 0.074 | -0.422 +/- 0.062 |
| R2 | 0.990 | 0.984 |
| level spacing r | 0.535 (GOE) | 0.654 > GOE |
| Lyapunov exponent | -0.224 (stable) | -0.082 (near-critical) |
| p null model | -- | < 0.001 |
| delta AIC (power law) | 25.6 | 23.0 |

## Repository structure

    TCSV/
    notebooks/
        layers_analysis.py      # Laminar spectral analysis
    results/                    # All figures (PNG)
    tcsv_paper_final.tex        # Manuscript (REVTeX4-2, PRE)
    README.md

## Data access

All data used is publicly available:

BBP:
  aws s3 cp s3://openbluebrain/Simulatable_Circuit/ConnectivityMatrices/ . --no-sign-request --recursive

MICrONS:
  aws s3 cp s3://openbluebrain/Simulatable_Circuit/ConnectivityMatrices/MICrONS_V1/connectivity_matrix.h5 . --no-sign-request

Electrophysiology:
  aws s3 ls s3://openbluebrain/Experimental_Data/Electrophysiological_recordings/ --no-sign-request

## Dependencies

  pip install numpy scipy matplotlib h5py pynwb

## Citation

Santos, G. (2025). Power-law spectral structure of synaptic connectivity
matrices is conserved across cortical areas, layers, and species.
Physical Review E (submitted). arXiv:[to be added]
