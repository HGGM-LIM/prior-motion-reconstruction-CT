# prior-motion-reconstruction-CT
This repository contains data, code and results for prior- and motion-based reconstruction (PRIMOR) method for respiratory gated CT presented in the paper **A novel prior- and motion-based compressed sensing methods for small-animal respiratory gated CT. JFPJ Abascal, M Abella, E Marinetto, J Pascau, and M Desco. Plos One, 2016 (in press).** DOI: 

We propose PRIMOR method that extends prior-based reconstruction (PBR) by including a model of the motion between gates. Motion is estimated using a nonrigid registration method based on hierarchical B-splines. PRIMOR solves 

![](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/PRIMOR_equation.jpg)

where Psi is the spatial gradient, which leads to TV, Phi is the wavelet transform, T is the temporal operator, F is the Radon transform, up is the prior image and v is the image variation with respect to the prior image. 

## Data 
Methods are assessed using small-animal respiratory gated-CT data. In this case, data is sorted among four respiratory gates, leading to few noisy and irregularly distributed projections. Data for the different scenarios used in the paper can be found at http://dx.doi.org/10.5281/zenodo.15685 

## Code
We provide MATLAB code for PBR and PRIMOR methods. A demo loads small-animal respiratory gated-CT data and reconstructed images and display results. To run PBR and PRIMOR methods you need the following packages: 

* IRT package: The image reconstruction toolbox is needed for forward and backprojection operators (J A Fessler, Image reconstruction toolbox [IRT], 2011, retrieved from <http://www.eecs.umich.edu/~fessler/code/index.html>).   
 
* Wavelab 850: Package for wavelet transform computation (Buckheit JB, Chen S, Donoho DL, Johnstone IM, Scargle JD. WaveLab. Reference Manual. ftp://playfair.stanford.edu/pub/wavelab/WaveLabRef.ps. 1995., retrieved from http://statweb.stanford.edu/~wavelab/).

* FFD-based registration package: For the motion estimation step we used the FFD-based registration package available from MATLAB Central (Dirk-Jan Kroon; B-spline grid,  image and point based registration; 2012, retrieved from http://www.mathworks.com/matlabcentral/fileexchange/20057-b-spline-grid-image-and-point-based-registration), 


## Summary of results ##
We provide videos of reconstructed images with FDK, prior-based reconstruction (PBR) and prior- and motion-based reconstruction (PRIMOR) for the different scenarios created by changing photon flux (low dose protocols) and number of projections (subsumpled scenarios)

-**Static protocol** (120 projections per gate, I0 where I0=4.5e4) reconstructed with PBR (left) and PRIMOR (right)

![](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_StaticProtocol_120p_I0_PBR.gif)
![](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_StaticProtocol_120p_I0_PRIMOR.gif)

-**Subsampled protocol** (40 projections, I0 where I0=4.5e4) reconstructed with PBR (left) and PRIMOR (right)

![](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_Subsampled_40p_I0_PBR.gif)
![](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_Subsampled_40p_I0_PRIMOR.gif)


-**Low dose protocol** (120 projections, I0/6 where I0=4.5e4) reconstructed with PBR (left) and PRIMOR (right)

![](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_LowDose_120p_07I0_PBR.gif)
![](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_LowDose_120p_07I0_PRIMOR.gif)


##  Repository files ##

The repository contains the following files:

### Data and MATLAB functions ###

- **Demo_PBR_PRIMOR_CT.m:** Demo that laods data and shows how to use PBR and PRIMOR methods 

- **PBR_CT.m:** PBR method

- **PRIMOR_CT.m:** PRIMOR method

- **ComputeSplineRegTwoDirectionsPoolDouble.m** Function to compute the registration between consecutive gates

- **data_I0By6_120p.mat:** Data example corresponding to a low dose protocol (120 projections, I0/6 where I0=4.5e4)

- **image_HighDose.mat:** Image of four respiratory gates for a high dose protocol (four times the standard dose) 

- **RefImage_I0By6_120p.mat:** Prior image computed from the average of all low-dose data 

- **G_linux.mat:** Projection operator from IRT software

- **TParameters_I0By6_120p.mat:** Registration result needed for PRIMOR method

- **FDK_I0By6_120p.mat:** Low-dose data reconstructed with FDK (first gate)

- **PBR_Rec_I0By6_120p.mat:** Low-dose data reconstructed with PBR

- **PRIMOR_Rec_I0By6_120p.mat:** Low-dose data reconstructed with PRIMOR

### Videos of results ###

#### High dose ####
- **ResGat_HighDose.gif:** Video of target images acquired using a high dose protocol (around four times the dose for static imaging, I0)

#### Static protocol ####

- **ResGat_StaticProtocol_120p_I0_FDK.gif:** Video of images reconstructed with FDK using a dose corresponding to the static  (120 projections per gate, I0 where I0=4.5e4)

- **ResGat_StaticProtocol_120p_I0_PBR.gif:** Video of images reconstructed with PBR using a dose corresponding to the static protocol

- **ResGat_StaticProtocol_120p_I0_PRIMOR.gif:** Video of images reconstructed with PRIMOR using a dose corresponding to the static protocol

#### Subsampled protocol ####
- **ResGat_Subsampled_40p_I0_FDK.gif:** Video of images reconstructed with FDK using a subsampled protocol (40 projections, I0 where I0=4.5e4)

- **ResGat_Subsampled_40p_I0_PBR.gif:** Video of images reconstructed with PBR using a subsampled protocol

- **ResGat_Subsampled_40p_I0_PRIMOR.gif:** Video of images reconstructed with PRIMOR using a subsampled protocol

#### Low dose protocol ####
- **ResGat_LowDose_120p_07I0_FDK.gif:** Video of images reconstructed with FDK using a low dose protocol (120 projections, I0/6 where I0=4.5e4)

- **ResGat_LowDose_120p_07I0_PBR.gif:** Video of images reconstructed with PBR using a low dose protocol

- **ResGat_LowDose_120p_07I0_PRIMOR.gif:** Video of images reconstructed with PRIMOR using a low dose protocol

If you use this code, please reference the publication JFPJ Abascal et al. A novel prior- and motion-based compressed sensing method for small-animal respiratory gated CT. Plos One, 2016 (in press). If you need to contact the author, please do so at mabella@hggm.es, juanabascal78@gmail.com

