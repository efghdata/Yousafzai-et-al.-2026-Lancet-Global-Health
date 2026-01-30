1.ProcessIncidenceDatasets calculates the following componenets of incidence:
- denominator - estimated number of children in the catchment area accounting for clusters that were not enumerated and households that could not be reached. Computed at the site level and overall as well as stratified by age group.
- Enrollment adjustment - weight accounting for children who were not screened (but otherwise could have been according to study guidelines) and children who were not enrolled but met study inclusion except for reasons that mean they could have potentially been counted as a case.
- Healthcare seeking adjustment - weight accounting for children in the community who had diarrhea but did not seek care. Ultimately replaced by a propensity to seek care model.
- Numerator - estimated number of children who tested positive for Shigella by culture or qPCR, stratfied by a variety of subgroups of interest. For analytic reasons, the healthcare seeking adjustment was incorporated at this point in the code for the relevant incidence indicators.
- Additional data created in this file are estimates for all diarrhea incidence as well as the proportions of qPCR and culture samples meeting particular serotypes and serogroups.
  
2.ProcessSecondaryDatasets computes all other data sources needed for the analysis:
- Description of the study population (table 1).
- Description of pre-screening, screening and enrollment as well as retention during longitudinal follow-up (figure 1).
- Proportions of culture samples that were resistant or showed intermediate resistance to common antibiotics, stratified by site and Shigella species.
- The proportion of qPCR samples attributable to other enteric pathogens from the TAC card among Shigella culture positive and Shigella qPCR-attributable samples and further resticted to samples from participants who presented with moderate or severe diarrhea.

3.RunAnalysis merges adjustments, numerators and denominators together and calculates the various incidences as described in the paper. All tables and figures showing primary data are generated here (With the exception of figure 1 and supplementary figure 1 which were made in powerpoint by transcribing data from the R summary).
