---
editor_options: 
  markdown: 
    wrap: 72
---

# krill-phyto

## Analysis of krill distributions in relation to phytoplankton and other dirvers

# Preliminary steps

-   use add_depth.r to add depth from Blake's grid to the data file.
    Need to update names etc. when doing so.

-   \_01_krill-sdm.r is the main analysis file.

## Methods

(1) Run base model with different distribution to test qq plots and fit

### Selecting distribution

## avNASC models (krill distributions)

(1) Ran four models with or w/o spatio-temporal fields but without
    phytoplankton for three distributions: (a) Tweedie, (b)
    delta-lognormal, and (c) delta-gamma distributions:

**Models**

paste0(spp," \~ 0 + f_year + s(scale_day, k = 3) + s(scale_depth, k =
3)")

QQplots for delta-lognormal showed divergence for both high and low
values.

QQplots for delta-gamma showed more divergence than Tweedie in the
abundance model with more under-prediction at the low end

QQplots the Tweedie distribution showed some over prediction for high
values, but was better than delta-lognormal and delta-gamma.

*Chose* Tweedie distribution as least-worse.

### Fitting

-   Run all models with and without spatiotemporal field:
-   phyto = diat_diatDino ratio for first set.

paste0(spp," \~ 0 + f_year"), paste0(spp," \~ 0 + f_year + s(scale_day,
k = 3)"),

paste0(spp," \~ 0 + f_year + s(scale_depth, k = 3)"),

paste0(spp," \~ 0 + f_year + s(scale_day, k = 3) + s(scale_depth, k =
3)") ,

paste0(spp," \~ 0 + f_year + s(scale_day, k = 3) + s(scale_depth, k =
3) + phyto"),

paste0(spp," \~ 0 + f_year + s(scale_day, k = 3) + s(scale_depth, k =
3) + s(phyto,k=3)") )

-   Build AIC and Sanity table.
-   Including iid and phyto tends to cause sanity problems (expected)
-   Then run cross-validation.



## diat_diatDino models (diatom/dinoflagelate ratio distributions)

- file contains some negative values ~ 20; deleted because they cause fitting issues.  

# diatom only models

# dinoflagellate only models
- ran base models for QQplots with spatial and spatiotemporal fields
-- difficulty with most models. Dropping spatial field allows fit with good Sanity outputs but QQplots for most distributions pretty marginal.



# Notes to self

**2025-07-29** Testing distributions for diat_diatDino ratio. There are
several below zero number that mess up QQplots. Might delete them. 20 in total. 
