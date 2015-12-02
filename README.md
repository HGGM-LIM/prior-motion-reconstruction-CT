# prior-motion-reconstruction-CT
Data, code and results for prior- and motion-based reconstruction (PRIMOR) method for respiratory gated CT

This repository contains the data files used in the paper: 
**“A novel prior- and motion-based compressed sensing methods for small-animal respiratory gated CT”. J F P J Abascal, M Abella, E Marinetto, J Pascau, and M Desco.** 
DOI: 

Small-animal respiratory gated-CT data, with four respiratory gates. We provide videos of reconstructed images with FDK, prior-based reconstruction (PBR) and prior- and motion-based reconstruction (PRIMOR). 

**SUMMARY OF RESULTS**

-**Static protocol** reconstructed with PBR (left) and PRIMOR (right)

![Image of Yaktocat](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_StaticProtocol_120p_I0_PBR.gif)
![Image of Yaktocat](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_StaticProtocol_120p_I0_PRIMOR.gif)

-**Subsampled protocol** reconstructed with PBR (left) and PRIMOR (right)

![Image of Yaktocat](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_Subsampled_40p_I0_PBR.gif)
![Image of Yaktocat](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_Subsampled_40p_I0_PRIMOR.gif)


-**Low dose protocol** reconstructed with PBR (left) and PRIMOR (right)

![Image of Yaktocat](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_LowDose_120p_07I0_PBR.gif)
![Image of Yaktocat](https://github.com/HGGM-LIM/prior-motion-reconstruction-CT/blob/master/ResGat_LowDose_120p_07I0_PRIMOR.gif)



**The repository contains the following files:**

VIDEOS OF RESULTS

HIGH DOSE
- **ResGat_HighDose.gif:** Video of target images acquired using a high dose protocol(around four times the dose for static imaging, I0)

STATIC PROTOCOL
- **ResGat_StaticProtocol_120p_I0_FDK.gif:** Video of images reconstructed with FDK using a dose corresponding to the static protocol

- **ResGat_StaticProtocol_120p_I0_PBR.gif:** Video of images reconstructed with PBR using a dose corresponding to the static protocol

- **ResGat_StaticProtocol_120p_I0_PRIMOR.gif:** Video of images reconstructed with PRIMOR using a dose corresponding to the static protocol

SUBSAMPLED PROTOCOL
- **ResGat_Subsampled_40p_I0_FDK.gif:** Video of images reconstructed with FDK using a subsampled protocol (40p, I0 where I0=4.5e4)

- **ResGat_Subsampled_40p_I0_PBR.gif:** Video of images reconstructed with PBR using a subsampled protocol

- **ResGat_Subsampled_40p_I0_PRIMOR.gif:** Video of images reconstructed with PRIMOR using a subsampled protocol

LOW DOSE PROTOCOL
- **ResGat_LowDose_120p_07I0_FDK.gif:** Video of images reconstructed with FDK using a low dose protocol (120p, I0/6 where I0=4.5e4)

- **ResGat_LowDose_120p_07I0_PBR.gif:** Video of images reconstructed with PBR using a low dose protocol

- **ResGat_LowDose_120p_07I0_PRIMOR.gif:** Video of images reconstructed with PRIMOR using a low dose protocol

If you need to contact the author, please do so at mabella (AT) mce (DOT) hggm (DOT) es.
