
function [u,errAll] = PBR_CT(G,f,R,N,uref,mu,lambda,alpha,beta,nBreg,varargin)
% u = PBR_CT(G,f,R,N,uref,mu,lambda,alpha,beta,nBreg)
% [u,errAll] = PBR_CT(G,f,R,N,uref,mu,lambda,alpha,beta,nBreg,uTarget)
%
% Inputs: 
%
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
%
% Outputs: 
%
% u         = reconstructed image of size n_x x n_y x n_t, it reconstruct
% all gates simultaneously 
% errAll    = relative solution error norm at each iteration, size nBreg x
% number of gates
%
% 
% Prior-based reconstruction (PBR) method is efficiently
% solved using the Split Bregman formulation. It assumes that a good
% quality (free of artefacts) prior image is available. The linear system
% in the image reconstruciton step is solved using  Gauss-Newton Krylov
% method (see the Krylov tolerance threshold below). 
%
%
% Requirements: 
% 
% To run PBR you need to download the IRT code and Wavelab
% 850 packages. This demo makes use of parallel computation using parfor loops 
% across different gates, but this is not a requirement (if you do not have
% the parallel computing toolbox, change parfor loops by for loops).  
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
%
% Faster implementations: 
%
% Change forward operator, G, and backprojection operator, G', with your
% own GPU implementations to make the algorithm up to orders of magnitude
% faster.  
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

if nargin >= 11
    uTarget     = varargin{1};
    errAll      = zeros(nBreg,N(3));
    uTarget     = double(uTarget)*normFactor;    
end % nargin

% Reserve memory for the auxillary variables
f0      = f;
u       = zeros(rows,cols,numTime);
x       = zeros(rows,cols,numTime);
y       = zeros(rows,cols,numTime);

bx      = zeros(rows,cols,numTime);
by      = zeros(rows,cols,numTime);

dx      = zeros(rows,cols,numTime);
dy      = zeros(rows,cols,numTime);

vWT     = zeros(NPadMax,NPadMax,numTime);
rhs_wt  = zeros(rows,cols,numTime);
dv      = zeros(NPadMax,NPadMax,numTime);
bv      = zeros(NPadMax,NPadMax,numTime);

murf    = zeros(rows,cols,numTime);
FFuref  = zeros(rows,cols,numTime);

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
    
    % Build RHS
    parfor iw = 1:numTime
        rhs_wt_temp     = IWT2_PO(dv(:,:,iw)-bv(:,:,iw),L,qmf);
        rhs_wt(:,:,iw)  = rhs_wt_temp(indPad,indPad);
    end % iw
    
    rhs     = murf-FFuref+lambda*rhs_wt+lambda*Dxt(x-dxref-bx)+...
        lambda*Dyt(y-dyref-by);
    
    % Reconstruct image difference
    v       = reshape(krylov(rhs(:)),N);
    
    % Images
    u       = uref + v;
    
    % Derivatives
    dx      = Dx(u);
    dy      = Dy(u);
    
    % WT of image difference
    parfor iw = 1:numTime
        vWT(:,:,iw)  = FWT2_PO(padarray(v(:,:,iw),[NPad NPad]),L,qmf);
    end % iw
    
    % update x, y, dv
    [x,y]   = shrink2(dx+bx,dy+by,beta/lambda);
    dv      = shrink1(vWT+bv,alpha/lambda);
    
    % update bregman parameters
    bv      = bv+vWT-dv;
    bx      = bx+dx-x;
    by      = by+dy-y;
    
    % Bregman iteration for the data constraint
    parfor iw = 1:numTime
        fForw         = double(G*u(:,:,iw)).*R(:,:,iw);
        f(:,:,iw)     = f(:,:,iw) + f0(:,:,iw)-fForw;
        murf(:,:,iw)  = double(mu*(G'*f(:,:,iw)));        
    end
    
     % Solution error norm
    if nargin >= 11
        for iw = 1:numTime
            errAll(outer,iw) = norm(reshape(uTarget(:,:,iw)-u(:,:,iw),[],1))...
                /norm(reshape(uTarget(:,:,iw),[],1));
        end
    end
    
    if any([outer ==1, outer == 10, rem(outer,10)==0])
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
        if nargin >= 11
            subplot(2,2,4); plot(errAll(1:outer,:)); axis tight; title(['Sol. error' ]);
        end
        colormap gray;
        drawnow;
    end % rem
end % outer
close(hw);

% undo the normalization so that results are scaled properly
u = u/(normFactor);

% =====================================================================
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
    function dx = krylov(r)
        %dx = gmres (@jtjx, r, 30, tolKrylov, 100);
        [dx,flag] = bicgstab(@jtjx, r, tolKrylov, 100);
    end

    % Callback function for matrix-vector product (called by krylov)
    function b = jtjx(sol)
        solMat  = reshape(sol,N);
        
        % TV part
        btv     = lambda*(Dyt(Dy(solMat)) + Dxt(Dx(solMat)));
        
        % Data constraint part
        bF      = zeros(N);
        
        % Jacobian u part
        parfor iq = 1:N(3)
            tmp            = R(:,:,iq).*(G*solMat(:,:,iq));
            bF(:,:,iq)     = mu*(double(G'*tmp));
        end
        
        % Prior part
        bwt     = lambda*sol;
        
        b       = bwt(:) + btv(:) + bF(:);
    end

end

%
