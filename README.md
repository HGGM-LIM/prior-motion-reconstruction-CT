# prior-motion-reconstruction-CT
This repository contains data, code and results for prior- and motion-based reconstruction (PRIMOR) method for respiratory gated CT presented in the paper: 
**“A novel prior- and motion-based compressed sensing methods for small-animal respiratory gated CT”. J F P J Abascal, M Abella, E Marinetto, J Pascau, and M Desco.** 
DOI: 

Small-animal respiratory gated-CT data used in this work is available from http://dx.doi.org/10.5281/zenodo.15685. 

We provide videos of reconstructed images with FDK, prior-based reconstruction (PBR) and prior- and motion-based reconstruction (PRIMOR) for the different scenarios created by changing photon flux (low dose protocols) and number of projections (subsumpled scenarios)

## Summary of results ##

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

**The repository contains the following files:**

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

If you need to contact the author, please do so at mabella (AT) mce (DOT) hggm (DOT) es.
