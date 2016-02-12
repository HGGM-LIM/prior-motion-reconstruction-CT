function [GridAll,SpacingAll,GridDAll,SpacingDAll] = ComputeSplineRegTwoDirectionsPoolDouble(u0,Options)
% [GridAll,SpacingAll,GridDAll,SpacingDAll] = ComputeSplineRegTwoDirectionsPoolDouble(u0,Options)
%
% Computes the registration between consecutive gates, forward and
% backwards to be used by PRIMOR method. 
%
% Inputs:
%
% u0        = image of all gates, size n_x x n_y x number gates
% Options   = parametes for the FFD-based registration package
% Options.Penalty     = 1e-4;
% Options.Registration= 'NonRigid';
% Options.MaxRef      = 3;
% Options.Similarity  = 'sd';
%
% Outputs: 
%
% GridAll,SpacingAll,GridDAll,SpacingDAll = mesh of grid points and
% spacing for forward and backward registration required by PRIMOR method
% 
% This version uses the double of processors for registration forward and
% backwards
%
% Requirements: 
%
% FFD-based registration package: For the motion estimation step we used
% the FFD-based registration method available in MATLAB Central (Dirk-Jan
% Kroon; B-spline grid,  image and point based registration; 2012,
% retrieved from 
% http://www.mathworks.com/matlabcentral/fileexchange/20057-b-spline-grid-image-and-point-based-registration), 
%
% Juan FPJ Abascal, Monica Abella
% Departamento de Bioingenieria e Ingenieria Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% mabella@hggm.es, juanabascal78@gmail.com, desco@hggm.es

dimIm       = size(u0);
numTime     = dimIm(3);

parfor ih = 1:2*numTime
    if ih <= numTime
        % Backward
        Imoving             = abs(u0(:,:,ih));
    
        % Up-transformation and down-transformations
        if ih == numTime
            [Ireg,Grid,Spacing,M,B,F]=image_registration(Imoving,abs(u0(:,:,1)),Options);
        elseif ih == 1
            [Ireg,Grid,Spacing,M,B,F]=image_registration(Imoving,abs(u0(:,:,ih+1)),Options);
        else
            [Ireg,Grid,Spacing,M,B,F]=image_registration(Imoving,abs(u0(:,:,ih+1)),Options);
        end % numtime

        GridAll(:,:,:,ih)   = Grid;
        SpacingAll(:,:,ih)  = Spacing;       
    else 
        % Forward transformation
        ih2                 = ih-numTime;
        Imoving             = abs(u0(:,:,ih2));
        
        if ih2 == 2*numTime
            [IregD,GridD,SpacingD,MD,BD,FD]=image_registration(Imoving,abs(u0(:,:,ih2-1)),Options);
        elseif ih2 == 1
            [IregD,GridD,SpacingD,MD,BD,FD]=image_registration(Imoving,abs(u0(:,:,end)),Options);
        else
            [IregD,GridD,SpacingD,MD,BD,FD]=image_registration(Imoving,abs(u0(:,:,ih2-1)),Options);        
        end % numtime       
        GridDAll(:,:,:,ih-numTime)  = GridD;
        SpacingDAll(:,:,ih-numTime) = SpacingD;        
    end % ih2
end % ih

