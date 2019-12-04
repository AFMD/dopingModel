# dopingModel

### Controlling energy levels and Fermi level en route to fully tailored energetics in organic semiconductors

Authors: Ross Warren<sup>1</sup>, Alberto Privitera<sup>1</sup>, Pascal Kaienburg<sup>1</sup>, Andreas E. Lauritzen<sup>1</sup>, Oliver Thimm<sup>2</sup>, Jenny Nelson<sup>3</sup> & Moritz Riede<sup>1</sup>

Affiliations:<br>
    <sup>1</sup> Clarendon Laboratory, Department of Physics, University of Oxford, Parks Road, Oxford OX1 3PU, UK.<br>
    <sup>2</sup>  IEK5-Photovoltaics, Forschungszentrum Jülich, 52425 Jülich, Germany. .<br>
    <sup>2</sup> Department of Physics, Imperial College London, Exhibition Road, London SW7 2AZ, UK.

This repository is associated with a study that has been published in Nature Communications (link tba). It was created so that our research is both transparent and reproducible.

The repo contains the data and code for:

  * PDS data (Figure 1)
  * EPR data (Figure 2)
  * Statistical model code (Figures 3, 4 and Supplementary Figures 6,7).

### Requirements

  * python3
  * cython
  * pandas

### How to use

 1. Clone or download this repository
 2. Build the statistical model cython module:

    `$ python setup.py build_ext --inplace`

 3. Run whichever script e.g.

    `$ python fig-4a-fixed-dopant-EA.py`

### Reference

If you find this code or data helpful, please reference it as: (tba)
