% Demo_PBR_PRIMOR_CT.m
% 
% Demo for reconstructing respiratory gated-CT data with a novel prior- and
% motion-based reconstruction (PRIMOR) method proposed in the publication 
% JFPJ Abascal, M Abella, E Marinetto, J Pascau, M Desco. A novel prior-
% and motion-based compressed sensing method for small-animal respiratory
% gated CT. Plos One, 2016 (in press). 
%
% Code downloaded from the repository 
% https://github.com/HGGM-LIM/prior-motion-reconstruction-CT
%
% PRIMOR extends prior-based reconstruction (PBR) method by including a
% model of the motion between gates. Motion is estimated using a nonrigid
% registration method based on hierarchical B-splines. A prior image is
% obtained as the average of all respiratory gates.  
%
% In this demo PRIMOR is compared to PBR using respiratory gated-CT data.
% In this case data is sorted among four respiratory gates, leading to few 
% noisy irregularly distributed projections. Video of the comparison
% between these methods can be found in the repository and data for the
% different scenarios used in the paper can be found at 
% http://dx.doi.org/10.5281/zenodo.15685
%
%
% Requirements: 
%
% The demo loads data and reconstructed images and display results. To run
% the reconstruction methods you need the following packages: 
% IRT code and Wavelab 850 package for PBR and PRIMOR methods, and 
% FFD-based registration software for PRIMOR. 
%
% IRT package: Simulated data and forward and backprojection operators have
% been computed using IRT code (J A Fessler, Image reconstruction toolbox
% [IRT], 2011, retrieved from
% <http://www.eecs.umich.edu/~fessler/code/index.html>).   
% Execute setup.m to add all directories to the path and execute
% mex_build.m to compile mex files for your system. 
% Bear in mind that these forward and backprojection operations are the
% most computationally expensive parts of the algorithm, so providing
% faster implementations will speed up the algorithm considerably. 
% 
% Wavelab 850: Package for wavelet transform computation (Buckheit JB, Chen
% S, Donoho DL, Johnstone IM, Scargle JD. WaveLab. Reference Manual.
% ftp://playfair.stanford.edu/pub/wavelab/WaveLabRef.ps. 1995., retrieved
% from http://statweb.stanford.edu/~wavelab/). Execute startup.m to ad
% paths and compile mex files.   
%
% FFD-based registration package: For the motion estimation step we used
% the FFD-based registration package available in MATLAB Central (Dirk-Jan
% Kroon; B-spline grid,  image and point based registration; 2012,
% retrieved from 
% http://www.mathworks.com/matlabcentral/fileexchange/20057-b-spline-grid-image-and-point-based-registration), 
% Add all directories to the path and run compile_c_files.m for faster
% implementation. 
%
% Parallel computation: This demo makes use of parallel computation using
% parfor loops across different gates, but this is not a requirement (if 
% you do not have the parallel computing toolbox, change parfor loops by
% for loops).  
%
%
% Faster implementations: 
% 
% Change forward operator, G, and backprojection operator, G', with your
% own GPU implementations to make the algorithm up to orders of magnitude
% faster.  
%
%
% If you use this code, please reference the paper JFPJ Abascal et al. A
% novel prior- and motion-based compressed sensing method for small-animal
% respiratory gated CT. Plos One, 2016 (in press). 
% 
% Juan FPJ Abascal, Monica Abella
% Departamento de Bioingenieria e Ingenieria Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% mabella@hggm.es, juanabascal78@gmail.com, desco@hggm.es

% -------------------------------------------------------------------------
% READ Data and forward operator G (mapping 2D image to the data)

% These have been computed using Fessler's IRT software (J.A. Fessler,
% Image reconstruction toolbox (IRT), 2011, retrieved from
% http://www.eecs.umich.edu/~fessler/code/index.html. 
%
% The following lines provides the parameters used to generate the matrix
% operator, G, using IRT:
% n_m         = 1;
% numPerProj  = 350;
% N           = [350 350];
% n_x         = N(1);
% numProj     = [360]
% totalAngle  = [360];
% binning     = 4;
% n_ang = numProj;
% det_z_count = 1; detector_count= n_x; pixel_count=n_x*n_x;
% ds          = 50e-4*binning; % pixel sixe= min pixel*binning
% dx          = ds/1.6;
% cg          = sino_geom('fan','nb',n_x, 'na', n_ang, ...
%     'ds',ds, 'orbit',-totalAngle, ...
%     'offset_s', 0.0, ...
%     'dsd', 35.2, 'dod', 13.2, 'dfs', inf);
% ig          = image_geom('nx', n_x, 'ny', n_x,'dx', dx);
% G           = Gtomo2_dscmex(cg, ig);
% geom        = fbp2(cg, ig);

% FORWARD OPERATOR
modeLinux       = 'y'; % change to 'y' if IRT software is installed
if modeLinux     == 'y'
    % It requires IRT software
    % Load Forward Operator
    load('G_linux','geom','G');   % In Linux system. If it does not work,
                                  % run the code above        
end

% READ SIMULATED DATA 
%
% HIGH DOSE DATA (ideal data). 
% Target image acquired from a high dose protocol (dosis four times
% the dosis used for static imaging)
load('image_HighDose','uTarget');
N           = size(uTarget);

% Display four respiratory gates
figure;
for it = 1:4
    imagesc(uTarget(:,:,it)); colormap gray; axis image; colorbar;
    title(['High dose, gate ' num2str(it)]);
    pause(0.3);
end

% STANDARD DOSE DATA: Standard dose for static imaging leads to irregularly
% distributed 
% Data for different scenarios can be found at http://dx.doi.org/10.5281/zenodo.15685
nameData    = 'data_I0By6_120p';         % Dose reduction by 6

load(nameData,'dataAll');
Nd          = size(dataAll);

% Undersampling pattern
RAll        = dataAll>0;
% -------------------------------------------------------------------------
% PRIOR image
%
% Prior image as the average image
Rsum        = sum(RAll,3);
dataAllAv 	= sum(dataAll,3);
dataAllAv(Rsum>0) = dataAllAv(Rsum>0)./Rsum(Rsum>0);
if 0 
    % To compute the reference image download the IRT code and run this
    % part
    uref 		= fbp2(dataAllAv, geom); uref(uref<0) = 0;
    gaussianFilter = fspecial('gaussian', [5, 5], 3); % [7 7], 5 stronger filter if the prior edges is a problem
    uref_aux = imfilter(uref, gaussianFilter, 'symmetric', 'conv');
    uref        = uref_aux;    
else
    % Load precomputed prior image
    load('RefImage_I0By6_120p','uref');
end
% -------------------------------------------------------------------------
% FBP reconstruction using IRT software
frameThis   = 1;
if 0 
    % FDK reconstruction, it requires IRT code 
    im          = fbp2(dataAll(:,:,frameThis), geom); im(im<0) = 0;
else
    % Load precomputed FDK reconstruction    
    load('FDK_I0By6_120p','im');
end
figure;
subplot(2,2,1);
imagesc(uTarget(:,:,frameThis)); title('Target gate 1 (high dose image)'); colorbar;
subplot(2,2,2);
imagesc(dataAll(:,:,frameThis)); title('Low dose data (120 proj., dose by 6)'); colorbar;
subplot(2,2,3);
imagesc(im*prod(Nd(1:2))/nnz(RAll(:,:,frameThis))); title('FDK reconstruction'); colorbar; axis image;
subplot(2,2,4);
imagesc(uref(:,:,frameThis)); title('Prior image: sum of four gates'); colorbar; colormap gray;

% -------------------------------------------------------------------------
% Create the support (circle)
[n_x_u,n_y_u,n_z_u] = size(uTarget);
[X,Y]       = meshgrid(linspace(1,n_x_u,n_x_u),linspace(1,n_x_u,n_x_u));
X           = X-n_x_u/2;
Y           = Y-n_x_u/2;
indCir      = find(sqrt(X.^2+Y.^2)<=n_x_u/2);
mask        = zeros(N(1:2));
mask(indCir)= 1;

% Apply mask to the prior image
uref        = uref.*mask;
% -------------------------------------------------------------------------
% PBR method
if 0
    % To run PBR method, run this part
    matlabpool(4); % comment if matlabpool not available
    mu          = 2;
    lambda      = 1;
    nBreg       = 25;
    alpha       = 0.4;
    beta        = 0.2;
    [uPbr,errPbr] = PBR_CT(G,dataAll,RAll,N,uref,mu,lambda,alpha,beta,nBreg,uTarget);    
    matlabpool close; 
else
    % Load PBR result already computed
    load('PBR_Rec_I0By6_120p','uPbr','errPbr');
end
% -------------------------------------------------------------------------
% PRIMOR method

% Registration step
if 0
    % To compute the registration, download the FFD registration software
    % and run this part (it takes ~2min)
    
    % Compute the registration step
    matlabpool(8); % comment if matlabpool not available
    
    % Parameters for the registration step
    Options.Penalty     = 1e-4;
    Options.Registration= 'NonRigid';
    Options.MaxRef      = 3;
    Options.Similarity  = 'sd';
    
    % Estimate the temporal operator, computed by the registration of
    % consecutive gates from previous reconstruction (we used PBR)
    u0          = uPbr;
    [GridAll,SpacingAll,GridDAll,SpacingDAll] = ComputeSplineRegTwoDirectionsPoolDouble(u0,Options);
    TParameters.GridAll         = GridAll;
    TParameters.SpacingAll      = SpacingAll;
    TParameters.GridDAll        = GridDAll;
    TParameters.SpacingDAll     = SpacingDAll;
    matlabpool close; 
else
    % Load registration result already computed
    load('TParameters_I0By6_120p','TParameters');
end

% Reconstruction step
mu          = 2;
lambda      = 1;
alpha       = 0.4;
beta        = 0.2;
gamma       = 0.5;
nBreg       = 50;
if 0
    % To run PRIMOR method download FFD registration software and IRT
    % software and execute this part  (it takes ~ 1h 20min)
    
    matlabpool(4); % comment if matlabpool not available               

    [uPrimor,errPrimor] = PRIMOR_CT(TParameters,G,dataAll,RAll,N,uref,mu,lambda,gamma,alpha,beta,nBreg,uTarget);
    % Reconstructed image and auxiliary variables are displayed for TV and
    % prior terms, for some iteration numbers. The number of nonzero
    % coefficients on the respective transformed domains are given as a
    % precentage
    matlabpool close;
else
    % Load PRIMOR result already computed
   load('PRIMOR_Rec_I0By6_120p','uPrimor','errPrimor');
end
% -------------------------------------------------------------------------
% Images for ideal image (high dose) and respiratory gated data with
% six-fold dose reduction reconstructed with FDK, PBR and PRIMOR methods
figure;
subplot(2,2,1);
imagesc(uTarget(:,:,1)); axis image; axis off; colormap gray;
title('FDK, high dose');
ca = caxis;
subplot(2,2,2);
imagesc(im(:,:,1)*prod(Nd(1:2))/nnz(RAll(:,:,frameThis))); axis image; 
axis off; colormap gray;
caxis(ca); title('FDK, low dose');
subplot(2,2,3);
imagesc(uPbr(:,:,1)); axis image; axis off; colormap gray;
caxis(ca); title('PBR, low dose');
subplot(2,2,4);
imagesc(uPrimor(:,:,1)); axis image; axis off; colormap gray;
caxis(ca); title('PRIMOR, low dose');

% Zoom image
xZoom        = 120:280;
yZoom        = 70:240;

figure;
subplot(2,2,1);
imagesc(uTarget(xZoom,yZoom,1)); axis image; axis off; colormap gray;
title('FDK, high dose');
ca = caxis;
subplot(2,2,2);
imagesc(im(xZoom,yZoom,1)*prod(Nd(1:2))/nnz(RAll(:,:,frameThis))); axis image; axis off; colormap gray;
caxis(ca); title('FDK, low dose');
subplot(2,2,3);
imagesc(uPbr(xZoom,yZoom,1)); axis image; axis off; colormap gray;
caxis(ca); title('PBR, low dose');
subplot(2,2,4);
imagesc(uPrimor(xZoom,yZoom,1)); axis image; axis off; colormap gray;
caxis(ca); title('PRIMOR, low dose');

% Convergence: solution error vs iteration number
figure; plot(mean(errPbr,2));
hold on; 
plot(mean(errPrimor,2),'r');
legend('PBR','PRIMOR'); 
xlabel('Iteration number'); ylabel('Solution error');

% -------------------------------------------------------------------------
%