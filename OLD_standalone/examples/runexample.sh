#MKM
../bin/monas -SMKM -DSMKM -fGetParticleContribution -fSetMultieventStatistic 10000 -topasScorer Scorer.phsp -MKM_rDomain 0.4571564 -MKM_rNucleus 8.00 -MKM_alpha 0.1626047 -MKM_beta 0.0789754 -MKM_y0 150 -MKM_alphaX 0.313 -MKM_betaX 0.030.6 -Doses 0 10 1 #HSG

#MKM
../bin/monas -MKMSatCorr -topasScorer Scorer.phsp -MKM_rDomain 0.44 -MKM_rNucleus 3.9 -MKM_alpha 0.19 -MKM_beta 0.07 -MKM_y0 150 -MKM_alphaX 0.313 -MKM_betaX 0.030.6 -Doses 0 10 1 #HSG

#GSM2
#./exe  -GSM2 -fSetMultiEventStatistic 100000 -fGetParticleContribution -topasScorer $topas -GSM2_rDomain 0.4 -GSM2_rNucleus 8 -GSM2_alphaX 0.19 -GSM2_betaX 0.05 -GSM2_a 0.03691869  -GSM2_b 0.1822144  -GSM2_r 3.640363  -GSM2_kappa 0.0750 -GSM2_lambda 0.0000750 -Doses 0 10 1 #HSG
