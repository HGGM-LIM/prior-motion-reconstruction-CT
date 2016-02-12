
function [u,errAll] = PRIMOR_CT(TParameters,G,f,R,N,uref,mu,lambda,gamma,alpha,beta,nBreg,varargin)
% u = PBR_CT(G,f,R,N,uref,mu,lambda,alpha,beta,nBreg)
% [u,errAll] = PBR_CT(G,f,R,N,uref,mu,lambda,alpha,beta,nBreg,uTarget)
%
% Inputs: 
%
% TParameters = result from the registration between gates (see
% ComputeSplineRegTwoDirectionsPoolDouble.m). It is a structure with fields
% TParameters.GridAll, TParameters.SpacingAll, TParameters.GridDAll,
% TParameters.SpacingDAll 
% G         = Projection operator where G' is the backprojection operator
% f         = data, N x numPro x time
% R         = undersampling pattern, same size as f. Entries are 1 when
% data is sampled and zero otherwise
% N         = image size n_x x n_y x n_t (number of pixels in spatial and
% temporal dimensions)
% uref      = prior image (usually the FBP of averaged data), size can be
% n_x x n_y x n_t or n_x x n_y (the latter will be replicated across time)
% mu        = 1, weight of data constraint. Increase for faster
% convergence, decrease for noisy data (1 usually works fine)
% lambda    = 1
% beta      = 0.2, TV sparsity parameter
% alpha     = 0.4, prior functional sparsity parameter,
% increase for imposing similarity to the prior image
% nBreg     = Bregman iterations. This is the most critical parameter,
% chose regarding the noise in the data (lower number of iterations for
% noisier data). It can be selected by checking the convergence (errAll)
% varargin  = {uTarget} = Gold standard to computer the error and assess
% convergence
%
% Outputs: 
%
% u         = reconstructed image of size n_x x n_y x n_t 
% errAll    = relative solution error norm at each iteration, size nBreg x
% number of gates
%
% Prior- and motion-based reconstruction (PRIMOR) method is efficiently
% solved using the Split Bregman formulation. It assumes that a good
% quality (free of artefacts) prior image is available and that there is an
% estimate of the temporal operator that relates gates, which is computed
% by a registration between consecutive gates from a previous
% reconstruction. The linear system in the image reconstruciton step is
% solved using  Gauss-Newton Krylov method.
%
%
% Requirements: 
%
% To run PRIMOR download the IRT code and FFD-based registration software.
% This demo makes use of parallel computation using parfor loops across
% different gates, but this is not a requirement (if you do not have the
% parallel computing toolbox, change parfor loops by for loops).  
%
% % IRT package: Simulated data and forward and backprojection operators have
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
%
% Faster implementations: 
% 
% Change forward operator, G, and backprojection
% operator, G', with your own GPU implementations to make the algorithm
% up to orders of magnitude faster. 
%
% If you use this code, please reference the paper JFPJ Abascal et al. A
% novel prior- and motion-based compressed sensing method for small-animal
% respiratory gated CT. Plos One, 2016 (in press). 
% 
% Juan FPJ Abascal, Monica Abella
% Departamento de Bioingenieria e Ingenieria Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% mabella@hggm.es, juanabascal78@gmail.com, desco@hggm.es

h       = figure;
hw      = waitbar(0);

% Wavelab WT parameters
L = 3;
qmf         = MakeONFilter('Symmlet',8);
NPadMax     = 2^ceil(log(N(1))/log(2));
NPad        = round((NPadMax-N(1))/2);
indPad      = NPad+1:NPadMax-NPad;

tolKrylov   = 1e-2;  % Krylov convergence criterion, the smaller the value
                     % the higher the precission for solving the linear system     

dimIm       = N;
rows        = N(1);
cols        = N(2);
numTime     = N(3);
Nf          = size(f);

if size(uref,3) <= N(3)
    uref        = repmat(uref,[1 1 N(3)]);
end

R           = reshape(R,Nf);
f           = f.*R;

% Normalize data
normFactor  = getNormalizationFactor(f,f);
f           = normFactor*f;
uref        = double(uref)*normFactor;

if nargin >= 13
    uTarget     = varargin{1};
    errAll      = zeros(nBreg,N(3));
    uTarget     = double(uTarget)*normFactor;
end % nargin

% Reserve memory for the auxillary variables
f0      = f;
v       = zeros(rows,cols,numTime);
u       = zeros(rows,cols,numTime);
x       = zeros(rows,cols,numTime);
y       = zeros(rows,cols,numTime);

bx      = zeros(rows,cols,numTime);
by      = zeros(rows,cols,numTime);

dx      = zeros(rows,cols,numTime);
dy      = zeros(rows,cols,numTime);

p       = zeros(rows,cols,numTime);
bp      = zeros(rows,cols,numTime);

vWT     = zeros(NPadMax,NPadMax,numTime);
rhs_wt  = zeros(rows,cols,numTime);
dv      = zeros(NPadMax,NPadMax,numTime);
bv      = zeros(NPadMax,NPadMax,numTime);

murf    = zeros(rows,cols,numTime);
FFuref  = zeros(rows,cols,numTime);

% Forward and backward spline-based transformations
% Forward (Grid,Spacing) must be the transformation of moving frame i to
% match the next frame, and backward (GridD,SpacingD) transform frame i to
% match the frame i-1
GridAll     = TParameters.GridAll;
SpacingAll  = TParameters.SpacingAll;
GridDAll    = TParameters.GridDAll;
SpacingDAll = TParameters.SpacingDAll;

dpref   = OperatorL(uref);
dxref   = Dx(uref);
dyref   = Dy(uref);

% Backprojection
parfor ip = 1:numTime
    murf(:,:,ip)    = double(mu*(G'*f(:,:,ip)));
    FFuref(:,:,ip)  = mu*(G'*(R(:,:,ip).*(G*uref(:,:,ip))));
end

%  Do the reconstruction
for outer = 1:nBreg;
    figure(hw); waitbar(outer/nBreg);
    parfor iw = 1:numTime
        rhs_wt_temp     = IWT2_PO(dv(:,:,iw)-bv(:,:,iw),L,qmf);
        rhs_wt(:,:,iw)  = rhs_wt_temp(indPad,indPad);
    end % iw
    rhs_wt  = lambda*rhs_wt;
    
    rhs_p   = lambda*OperatorLt(p-dpref-bp);
    rhs_tv  = lambda*(Dxt(x-dxref-bx)+Dyt(y-dyref-by));
    rhs     = murf-FFuref+rhs_wt+rhs_p+rhs_tv;
    v       = reshape(krylov(rhs(:)),N);
    u       = uref + v;
    
    % Derivatives
    dx      = Dx(u);
    dy      = Dy(u);
    dp      = OperatorL(u);
    
    % WT of image difference
    parfor iw = 1:numTime
        vWT(:,:,iw)  = FWT2_PO(padarray(v(:,:,iw),[NPad NPad]),L,qmf);
    end % iw
    
    % update x, y, p
    [x,y]   = shrink2(dx+bx, dy+by,alpha/lambda);
    dv      = shrink1(vWT+bv, gamma/lambda);
    p       = shrink1(dp+bp, beta/lambda);
    
    % update bregman parameters
    bv      = bv+vWT-dv;
    bx      = bx+dx-x;
    by      = by+dy-y;
    bp      = bp+dp-p;
    
    % Bregman iteration for the data constraint
    parfor iw = 1:numTime
        fForw         = double(G*u(:,:,iw)).*R(:,:,iw);
        f(:,:,iw)     = f(:,:,iw) + f0(:,:,iw)-fForw;
        murf(:,:,iw)  = double(mu*(G'*f(:,:,iw)));  
    end
    
    % Solution error norm
    if nargin >= 13
        for iw = 1:numTime
            errAll(outer,iw) = norm(reshape(uTarget(:,:,iw)-u(:,:,iw),[],1))...
                /norm(reshape(uTarget(:,:,iw),[],1));
        end
    end
        
    if any([outer ==1, outer == 5, rem(outer,10)==0])
        % Display image and auxiliary variables for TV and prior terms. The
        % number of nonzero coefficients on the respective transformed
        % domains are given as a precentage
        subplot(2,2,1);
        imagesc(u(:,:,1)/(normFactor)); title(['Iter ' num2str(outer)]); colorbar; axis image;
        subplot(2,2,2);
        imagesc(x(:,:,1)); axis image; title(['xw, sparsity % ' num2str(100*nnz(x(:,:,1))/(prod(N(1:2))))]);
        subplot(2,2,3);
        tmp     = dv(:,:,1);
        tmp     = IWT2_PO(tmp,L,qmf);
        tmp     = tmp(indPad,indPad);
        imagesc(tmp); colorbar; axis image; title(['w, sparsity % ' num2str(100*nnz(dv(:,:,1))/(NPadMax*NPadMax))]);
        if nargin >= 13
            subplot(2,2,4); plot(errAll(1:outer,:)); axis tight; title(['Sol. error' ]);
        end
        colormap gray;
        drawnow;
    end % rem
            
end % outer
close(hw);

% undo the normalization so that results are scaled properly
u = u/normFactor;

% =====================================================================
    function uL     = OperatorL(uThis)
        % Apply the transform operator Lu_ijk=u_ijk-u_i'j'k-1 that
        % associates a pixel ij in the frame k with a unique pixel (using
        % Nearest Neighbour interpolation) in the frame k, the
        % interpolation matrix is obtained from the spline-registratio
        % toolbox.
        % NN(:,:,1) transforms a pixel ij in frame 1 to another i'j' in
        % frame 2
        %
        % This is as Tu_k = u_k-T_{k-1}(u_{k-1})
        % Each frame is substracted the forward tranformtation of the
        % previous frame
        
        uTAll       = InterpolateBasedSplines(uThis,GridAll,SpacingAll);
        uL          = zeros(dimIm);
        
        % Compare to previous frame
        for it = 1:dimIm(3)
            if it == 1
                LuThis         = uTAll(:,:,end);
            else
                LuThis         = uTAll(:,:,it-1);
            end
            
            uL(:,:,it)         = uThis(:,:,it)-LuThis;
        end
    end

    function uLt    = OperatorLt(uThis)
        %         % Apply the transpose of the transform operator
        %         % Lu_ijk=u_ijk-u_i'j'k-1
        %
        % This is as Ttu_k = u_k-Du_k+1
        % Each frame is substracted the backward tranformtation of the next
        % frame
        
        uTAll       = InterpolateBasedSplines(uThis,GridDAll,SpacingDAll);
        uLt         = zeros(dimIm);
        
        % Compare to previous frame
        for it = 1:dimIm(3)
            if it == dimIm(3)
                LuThis      = uTAll(:,:,1);
            else
                LuThis      = uTAll(:,:,it+1);
            end
            
            uLt(:,:,it)     = uThis(:,:,it)-LuThis;
        end
    end

    function uTAll = InterpolateBasedSplines(uThis,GridAll,SpacingAll)
        % Transform each image corresponding to each frame, following the
        % spline-based transformation given by GridAll and SpacingAll. This
        % should transform each frame to be similar to the next one
        
        uTAll           = zeros(dimIm);
        for ih = 1:numTime
            % Tranform u:
            IT             = uThis(:,:,ih);
            uTAll(:,:,ih)  = bspline_transform(GridAll(:,:,:,ih),IT,SpacingAll(:,:,ih));
        end % ih
        
    end

    function normFactor = getNormalizationFactor(R,f)
        
        normFactor = 1/norm(f(:)/size(R==1,1));
        
    end

    function d = Dx(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,2:cols,:) = u(:,2:cols,:)-u(:,1:cols-1,:);
        d(:,1,:) = u(:,1,:)-u(:,cols,:);
    end

    function d = Dxt(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,1:cols-1,:) = u(:,1:cols-1,:)-u(:,2:cols,:);
        d(:,cols,:) = u(:,cols,:)-u(:,1,:);
    end

    function d = Dy(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(2:rows,:,:) = u(2:rows,:,:)-u(1:rows-1,:,:);
        d(1,:,:) = u(1,:,:)-u(rows,:,:);
    end

    function d = Dyt(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(1:rows-1,:,:) = u(1:rows-1,:,:)-u(2:rows,:,:);
        d(rows,:,:) = u(rows,:,:)-u(1,:,:);
    end

    function [xs,ys] = shrink2(x,y,lambda)
        s = sqrt(x.*conj(x)+y.*conj(y));
        ss = s-lambda;
        ss = ss.*(ss>0);
        s = s+(s<lambda);
        ss = ss./s;
        xs = ss.*x;
        ys = ss.*y;
    end

    function xs = shrink1(x,lambda)
        s = abs(x);
        xs = sign(x).*max(s-lambda,0);
    end

    % Krylov solver subroutine to solve the linear system
    % X = GMRES(A,B,RESTART,TOL,MAXIT,M)
    % bicgstab(A,b,tol,maxit)
    function dx = krylov(r)
        %dx = gmres (@jtjx, r, 30, tolKrylov, 100);
        [dx,flag] = bicgstab(@jtjx, r, tolKrylov, 100);
    end

    % Callback function for matrix-vector product (called by krylov)
    function b = jtjx(sol)
        solMat  = reshape(sol,N);
        
        % TV part
        btv     = lambda*(Dyt(Dy(solMat)) + Dxt(Dx(solMat)));
        
        % Temporal operator part
        bP      = lambda*OperatorLt(OperatorL(solMat));
        
        % Data constraint part
        bF      = zeros(N);
        parfor iq = 1:N(3)
            tmp            = R(:,:,iq).*(G*solMat(:,:,iq));
            bF(:,:,iq)     = mu*(double(G'*tmp));
        end
        
        % Prior part
        bwt     = lambda*sol;
        
        b       = bwt(:) + btv(:) + bP(:) + bF(:);
    end

end

%
