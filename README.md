# Coupling SDMs from Alaska and Canada, and calculating initial biomass for Atlantis GOA

This code performs the following tasks:

1. Calculate values of initial biomass for `init.nc` based on 1990 (or closest year) biomass estimates from stock assessments or from Aydin et al. (2007), expanded to account for unassessed age classes and for biomass in British Columbia.
2. Calculate estimates of biomass per box that we then use to calculate the spatial distributions S1-S4.
3. Writes out S1-S4 vectors.

Biomass in Alaska is taken from stock assessments or estimates per unit area from Aydin et al. (2007) expanded to the model domain area. Biomass in British Columbia is taken from one of the following, in this order depending on availability:
- DFO stock assessments 
- [RAM Legacy Stock Assessment Database](https://www.ramlegacy.org/)
- Assumption of equal density at depth with SE GOA

![](https://github.com/somros/Atlantis_GOA_SDM_coupling_and_initbiom/blob/06264a58309c6f09a2720ec931116ae563409d0e/assessments%20and%20sdms.png) 

Biomass is allocated to Atlantis boxes based on: SDM results where these satisfy a set of validation criteria; spatial distributions derived from alternative sources (e.g. for species not captured in the bottom trawl data); and heuristics like habitat type or distance from shore. See scripts for details.

__References__

Aydin, K. Y., Gaichas, S., Ortiz, I., Kinzey, D., & Friday, N. (2007). A comparison of the Bering Sea, Gulf of Alaska, and Aleutian Islands large marine ecosystems through food web modeling. NOAA Technical Memorandum NMFS-AFSC. No. 178, 178, 1–298.

