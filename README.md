# SocialMap
Analysis code for Map Making: Constructing, Combining, and Inferring on Abstract Cognitive Maps https://doi.org/10.1016/j.neuron.2020.06.030

1. Univariate analysis

1.1. Preprocessing.m
- fMRI preprocessing using SPM12

1.2. Level1.m
- 1st level analysis (GLM)
- See Model files for model specification (/model/Model_*.m)
- See Contrast files for contrast specification (/model/Cont_*.m)

1.3. Level2.m
- 2nd level analysis (Group level one sample t-test)

2. Representational similarity analysis (RSA) in the region of interest (ROI)

2.1. RSA_1ExtctRaw.m
- Extract raw signals from an ROI

2.2. RSA_2NoiseCov.m
- Noise nomalization
- Computing representational dissimilarity matrix (RDM)

2.3. RSA_3KendallT.m
- Computing rank correation (Kendall's tau) over the permutation

3. Representational similarity analysis (RSA) in the region of interest (ROI)

3.1. RSA_sL1_RDMcrs.m
- Create a RDM matrix per searchlight

3.2. RSA_sL2_AppRSA.m
- Computing rank correation (Kendall's tau) over the permutation

3.3. RSA_sL3_smoothing.m
- Creating a map from rank correation and smooting
