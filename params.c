#Main parameters
#Type of fit (BAOISO/BAOANISO/FS/FSBAOISO/FSBAOANISO/FSalphasrecon): FS

#Include power spectrum (yes/no): yes
#FIT BAO to P0/P2/P4/P02/P024/P24/P04: P02
#k-range for P0 computation (double/double): 0.02 0.30
#k-range for P2 computation (double/double): 0.02 0.30
#k-range for P4 computation (double/double): 0.02 0.30
#FIT FS to P0/P2/P4/P02/P024/P24/P04: P024
#k-range for P0 computation (double/double): 0.02 0.15
#k-range for P2 computation (double/double): 0.02 0.15
#k-range for P4 computation (double/double): 0.02 0.15
#Number of chunks: 2

#Read input parameters for Power Spectrum
#Path to data1BAO: none
#Path to data2BAO: none
#Path to data1FS: /home/hector/eboss/lrg_v72/bispectrum/data_30jan/Power_Spectrum_comb_NGC_dataNorm_datav72_Dk8.txt
#Path to data2FS: /home/hector/eboss/lrg_v72/bispectrum/data_30jan/Power_Spectrum_comb_SGC_dataNorm_datav72_Dk8.txt
#Path to mocks1BAO/cov1: none
#Path to mocks2BAO/cov1:  none
#Path to mocks1FS: /home/hector/eboss/lrg_v7/bispectra/ezmocks/Power_Spectrum_comb_NGC_ezmocks_
#Path to mocks2FS: /home/hector/eboss/lrg_v7/bispectra/ezmocks/Power_Spectrum_comb_SGC_ezmocks_
#number of realizations: 1000
#Force to have diagonal covariance (yes/no): no
#Path to mask1 (none if no mask): /home/DATA/hector/boss_eboss_final/lrg_eboss/window/W_randomdata_eBOSS_LRGpCMASS_NGC_v72_yama_newformat.txt
#Path to mask2 (none if no mask): /home/DATA/hector/boss_eboss_final/lrg_eboss/window/W_randomdata_eBOSS_LRGpCMASS_SGC_v72_yama_newformat.txt
#Renormalize window (yes/no): yes
#Apply mask as a matrix (yes/no): yes
#Correction by (Hartlap/Sellentin-Heavens/none): Hartlap
#Covariance for FSalpharecon (fixed/varying): varying

#Read input parameters for Bispectrum
#Include bispectrum (yes/no): yes
#Fit bispectrum (B) or reduced bispectrum (Q): B
#Path to data1BAO: none
#Path to data2BAO: none
#Path to mocks1BAO: none
#Path to mocks2BAO: none
#k-range for B0 computation (double/double): 0.02 0.30
#Path to data1FS: /home/hector/eboss/lrg_v72/bispectrum/data_30jan/Bispectrum_comb_NGC_dataNorm_datav72_Dk8.txt
#Path to data2FS: /home/hector/eboss/lrg_v72/bispectrum/data_30jan/Bispectrum_comb_SGC_dataNorm_datav72_Dk8.txt
#Path to mocks1FS: /home/hector/eboss/lrg_v7/bispectra/ezmocks_Dk8/Bispectrum_comb_NGC_ezmocks_
#Path to mocks2FS: /home/hector/eboss/lrg_v7/bispectra/ezmocks_Dk8/Bispectrum_comb_SGC_ezmocks_
#k-range for B0 computation (double/double): 0.02 0.15

#Do Compression? (no/linear/nonlinear): no
#Path to compression file: compression_file_format.txt

#BAO Polynomial fit
#Order of polynomial broadband: 5
#Sigma variables (effective/para-perp): para-perp
#Sigma variables independent (yes/no): yes
#f-factor for Sigma dependence: 0.8204
#Smoothing scale in Mpc/h (only for BAOANISO): 15

#Sigma_variables
#SigmaP_0 (free/fixed/gaussian): fixed
#SigmaP_0 (mean/devstd): 5 1
#SigmaP_2 (free/fixed/gaussian): fixed
#SigmaP_2 (mean/devstd): 10 1
#SigmaP_4 (free/fixed/gaussian): fixed
#SigmaP_4 (mean/devstd): 9 3
#SigmaB_0 (free/fixed/gaussian): fixed
#SigmaB_0 (mean/devstd): 9 3
#Sigma_para (free/fixed/gaussian): fixed
#Sigma_para (mean/devstd): 7.0 3
#Sigma_perp (free/fixed/gaussian): fixed
#Sigma_perp (mean/devstd): 2.0 3

#FS priors/options
#Apply noise in the model or subtract it in the data (model/data): data
#Anoise (flat/gaussian): flat
#Anoise (mean/devstd): 1 0.30
#b2 (flat/gaussian): flat
#b2 (mean/devstd): 5 2.5
#b2s2 (flat/gaussian): flat
#b2s2 (mean/devstd): 0 5
#b3nl (flat/gaussian): flat
#b3nl (mean/devstd): 0 5

#Priors on alpha: 0.7 1.3
#Path to smooth Plin: /home/DATA/hector/boss_eboss_final/ptcool/bao/Pk_smkirkby_eboss_comb_z070_matterpower_15.txt
#Path to Olin: /home/DATA/hector/boss_eboss_final/ptcool/bao/Olinkirkby_eboss_comb_z070_matterpower_15.txt
#Type of fit (mcmc/analytic): mcmc

#For mcmc fit
#Number of threads: 40
#Use proposal covariance (yes/no): yes
#Path to proposal covariance: results/mcmcFS_output_test7b.txt
#Maximum number of mcmic accepted steps: 4000000
#mcmc sampling step (recommended 1.9): 1.9

#For BAO anaylytical fit
#Interval for alpha: 0.00001

#For FS mcmc fit
#Input file from PTcool: /home/DATA/hector/boss_eboss_final/ptcool/Perturbation_theory_eboss_comb_z070_matterpower.txt
#Type of PT-Model (linear/1L-SPT/2L-SPT/1L-RPT/2L-RPT): 1L-SPT
#Type of RSD-model (Kaiser87/Scoccimarro04/TNS10): TNS10
#Type of FoG-model (Exponential/Lorentzian/Exponential_avir): Lorentzian
#Type of model for bispectrum (tree-level/1L-SPT/GilMarin14/GEO23): GilMarin14
#Redshift input for GEO model: 0.5
#Lagrangian local bias b2s2 (yes/no/off): yes
#Lagrangian local bias b3nl (yes/no/off): yes
#RSD-fit (yes/no/shape/shape2): yes
#sigma8 free parameter (yes/no): no
#FoG free parameter (yes/no): yes
#Same FoG parameter for bispectrum (yes/no): no

#Read out parameters
#Path to output: ./results
#Identifier of output: test_window
#Plot best-fit (yes/no): yes
