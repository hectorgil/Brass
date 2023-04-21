#Main parameters
#Type of fit (BAOISO/BAOANISO/FS/FSBAOISO/FSBAOANISO/FSalphasrecon): //This is the type of fit. BAOISO and BAOANISO stands for two type of BAO-only fits. FS stands for RSD or FullShape type of fit. The FSBAOISO, FSBAOANISO consists of simultaneous fits to BAO (post-recon in this case) and FS (to prerecon catalogues). Finally FSalphasrecon uses the datavector consisting of P(k) prerecon + alpha_parallel and alpha_perp derived from the post-recon BAO type of fit

#Include power spectrum (yes/no): //if the power spectrum data is included
#FIT BAO to P0/P2/P4/P02/P024/P24/P04: //Which multipoles are going to be considered for a BAO type of fit
#k-range for P0 computation (double/double):  //kmin and kmax of the monopole range for BAO
#k-range for P2 computation (double/double): //kmin and kmax of the quadrupole range for BAO
#k-range for P4 computation (double/double):  //kmin and kmax of the hexadecapole range for BAO
#FIT FS to P0/P2/P4/P02/P024/P24/P04:  //Which multipoles are going to be considered for a FS type of fit
#k-range for P0 computation (double/double): //kmin and kmax of the monopole range for FS
#k-range for P2 computation (double/double): //kmin and kmax of the quadrupole range for FS
#k-range for P4 computation (double/double): //kmin and kmax of the hexadecapole range for FS
#Number of chunks: //number of independent chunks of data. Can be only 1 or 2

#Read input parameters for Power Spectrum
#Path to data1BAO:  //path to 1st chunk of data for BAO analysis
#Path to data2BAO:   //path to 2nd chunk of data for BAO analysis
#Path to data1FS:  //path to 1st chunk of data for FS analysis
#Path to data2FS:  //path to wnd chunk of data for FS analysis
#Path to mocks1BAO/cov1: //path to the mock file name containing all the mocks used for covariance, for the 1st chunk of BAO
#Path to mocks2BAO/cov1: //same as above but for the 2nd chunck for BAO
#Path to mocks1FS: //same as above but for the 1st chunck for FS
#Path to mocks2FS: //same as above but for the 2nd chunck for FS
#number of realizations: //number of mocks used for covariance. 0 for analytical covariance
#Path to mask1 (none if no mask): //path to the mask file for the 1st chunk
#Path to mask2 (none if no mask): //same as a above, but for the 2nd chunck
#Renormalize window (yes/no): //whether the mask should be renormalized matched to the data amplitude
#Apply mask as a matrix (yes/no): //whether mask should be applied as a matrix, or as a direct FF in each mcmc step
#Correction by (Hartlap/Sellentin-Heavens/none): //type of covariance correction to be added
#Covariance for FSalpharecon (fixed/varying): //whether the covariance for FSalphpa should vary accros realization, or is fuly fixed by the mocks. 

#Read input parameters for Bispectrum
#Include bispectrum (yes/no): //whethter bispectrum signal is used
#Fit bispectrum (B) or reduced bispectrum (Q): //whether the fit should be done on the bispectrum, or to the reduced bispectrum
#Path to data1BAO: //path of data for bispectrum for a BAO fit. 1st chunck
#Path to data2BAO: //same as above, for 2nd chunk
#Path to mocks1BAO: //path to the mocks for covariance for the BAO type of fit on the 1st chunk
#Path to mocks2BAO: //same as abaove for the 2nd chunck
#k-range for B0 computation (double/double): //kmin kmax for the BAO bispectrum fit
#Path to data1FS: //path of data for bispectrum for a FS fit. 1st chunck
#Path to data2FS: //same as above, for 2nd chunk
#Path to mocks1FS: //path to the mocks for covariance for the FS type of fit on the 1st chunk
#Path to mocks2FS: //same as abaove for the 1stt chunck
#k-range for B0 computation (double/double): //kmin kmax for the FS bispectrum fit

#Do Compression? (no/linear/nonlinear): //whether a compression should be applied (not implemented)
#Path to compression file: //file for the compresion (not implemented)

#BAO Polynomial fit
#Order of polynomial broadband: //For BAO fits, order of the broad-band polynomial to be used
#Sigma variables (effective/para-perp): //whether the BAO-damping sigmas are parametrized as parallel-perpendicular to the LOS, or as an effective values on monopole, quadrupole, and hexadecapole signals
#Sigma variables independent (yes/no): //whether the different sigmas are or not independent
#f-factor for Sigma dependence: //if they are not independent, which growth parameter relates them
#Smoothing scale in Mpc/h (only for BAOANISO): //For BAOANISO fits, the smoothing scale usded for reconstruction

#Sigma_variables
#SigmaP_0 (free/fixed/gaussian): //If using and effective sigma for the monopole, if its free, fixed or has a gaussian prior
#SigmaP_0 (mean/devstd): //if it's fixed or has a gaussian prior, here its mean (or fixed value) and error is given
#SigmaP_2 (free/fixed/gaussian): //same for effective sigma for quadrupole
#SigmaP_2 (mean/devstd): //same for effective sigma for quadrupole
#SigmaP_4 (free/fixed/gaussian): //same for effective sigma for hexadecapole
#SigmaP_4 (mean/devstd): //same for effective sigma for hexadecaple
#SigmaB_0 (free/fixed/gaussian): //same for the effective sigma for the bispectrum
#SigmaB_0 (mean/devstd): //same for the effective sigma for the bispectrum 
#Sigma_para (free/fixed/gaussian): //same for sigma_parallel type of sigma
#Sigma_para (mean/devstd): //mean (or fixed value) and error for he sigma_parallel
#Sigma_perp (free/fixed/gaussian): //same for alpha_perp
#Sigma_perp (mean/devstd): //same as above for alpha_perp

#FS priors/options
#Apply noise in the model or subtract it in the data (model/data): //if the noise correction is applied on the data (no affected by the window), or on the model (convolved with the window). If no window no diffference
#Anoise (flat/gaussian): //type of prior for Anoise (FS fits)
#Anoise (mean/devstd): //center and error for the Anoise prior if gaussian
#b2 (flat/gaussian): //same for the b2 arameter
#b2 (mean/devstd): //same for the b2 parameter
#b2s2 (flat/gaussian): //same for the bs2 parameter (if free)
#b2s2 (mean/devstd): //same for the bs2 parameter (if free)
#b3nl (flat/gaussian): //same for the b3nl parameter (if free)
#b3nl (mean/devstd): //same for the b3nl parameter (if free)

#Priors on alpha: //priors or range on alphas to be explored
#Path to smooth Plin: //path for the Plin used for the BAO fit
#Path to Olin: //path for the Plin/Pnl=Olin file used for the BAO fit
#Type of fit (mcmc/analytic): //whether the BAO parameters are marginalized analytically, or via a mcmc

#For mcmc fit
#Number of threads: //number of threads to be used
#Use proposal covariance (yes/no): //if a proposal covariance is to be used or the values from the prior.c file will be taken for the mcmc exploration
#Path to proposal covariance: //file to the proposal covariance
#Maximum number of mcmic accepted steps: //maximum number of mcmc steps to be run
#mcmc sampling step (recommended 1.9): //size of the step for mcmc. 

i#For BAO anaylytical fit
#Interval for alpha: //for analytical fits, the size of the alpha-grid to be explored

#For FS mcmc fit
#Input file from PTcool: //path to the PTCool file for the FS fits
#Type of PT-Model (linear/1L-SPT/2L-SPT/1L-RPT/2L-RPT): // type of model for tthe FS fit
#Type of RSD-model (Kaiser87/Scoccimarro04/TNS10): //type of model for the FS fit
#Type of FoG-model (Exponential/Lorentzian/Exponential_avir): //type of model for the FS fit
#Type of model for bispectrum (tree-level/1L-SPT/GilMarin14): //type of model for the FS fit for the bispecturm
#Lagrangian local bias b2s2 (yes/no/off): //if b2s2 it is fixed to the local prediction (yes), it's a free parameter (no) or it's set to 0 (off)
#Lagrangian local bias b3nl (yes/no/off): //same for b3nl parameter
#RSD-fit (yes/no/shape/shape2): //type of fit. 'no' for alphas set to 1, f set to 0. 'yes' for alpha's and f varying, but with no shape. 'shape' if m has to vary, 'shape2' for m and n variation (along with alphas and f)
#sigma8 free parameter (yes/no): //if sigma8 is to be varied
#FoG free parameter (yes/no): //if sigma_FOG is free or it's set to 0 (no)
#Same FoG parameter for bispectrum (yes/no): //if the Sigma_FOG from P and B are the same

#Read out parameters
#Path to output: //path of the output (need writting permisions)
#Identifier of output: //unique identifier for your output files. The string . is not allow and should NEVER be used as part of the identifier
#Plot best-fit (yes/no): //if file with best-fitting file is to be printed
