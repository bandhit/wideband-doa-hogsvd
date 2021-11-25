function [comp,startval,ssqres,expvar,scaling,tuckcongr,t3core] = nwaydecomp_spacetime(dat,ncomp,F,varargin)

% NWAYDECOMP_SPACETIME is a low-level function for nd_nwaydecomposition and is used
% to perform a decomposition of Fourier coefficients over space, frequency, epochs, and tapers,
% according to the SPACE-time model presented in:
%
%    van der Meij R, Jacobs J, Maris E (2015). Uncovering phase-coupled oscillatory networks in
%        electrophysiological data. Human Brain Mapping
%
% For details on the model specification and the algorithm implemented below, see the above publication.
% The following is a very succint description.
%
% This function decomposes a channel-by-frequency-by-epoch-by-taper (JxKxLxM) 4-way array into N components.
% Each component consists of a:
% A: spatial amplitude map (Jx1)
% B: frequency profile (Kx1)
% C: epoch profile (Lx1)
% Sigma: spatial time-delay maps(Jx1)
% D: Between-component coherency (the regular transpose of D in the above paper). (also denoted as phi here)
%
% Matrices A,B,C,S are of size JxN, KxN, LxN, JXN resp.
% Matrix D is of size KxNxN (if Dmode = 'kdepcomplex') or NxN (if Dmode = 'identity')
%
% Use as
%   [comp,ssqres,expvar,scaling,tuckcongr,t3core] = nwaydecomp_spacetime(dat,ncomp,F,...)
%
% Input:
%   dat   = 4-way array channel-by-frequency-by-epoch-by-taper of Fourier coefficients to be decomposed
%   ncomp = number indicating number of components
%   F     = 1xnfreq vector indicating frequencies in Hz present in 2nd dimension
%
% Output:
%         comp = 1x5 cell-array containing the estimated component parameters (A, B, C, Sigma, D)
%     startval = initialization values of the algorithm
%       ssqres = sums of squares of the residuals
%       expvar = percentage explained variance of the input data by the model
%      scaling = componont-specific scaling coefficients 
%                A,B,C have a component-specific vector of norm = 1. Sigma has a component-specific mean of half 
%                the 'circularity point' (see paper above), weigthed by A. (The phase of D is shifted accordingly if complex-valued).
%    tuckcongr = tuckers congruence coefficents between components, high values mean high correlation between components, which is a sign of
%                a degenerate model
%       t3core = a Tucker3 model core-array (vectorized). If the above model holds well, the t3core corresponds to a 4-dimensional core-array 
%                of zeros with only 1 on its super-diagonal. When the off-diagonal terms become non-zero, it can be assumed components are 
%                fitting structure that doesn't follow the model, i.e. noise. The t3core is used to calculate the core consistency diagnostic.
%
%
% Additional options should be specified in key-value pairs and can be
%   'niter'        = maximum number of iterations (default = 2500)
%   'convcrit'     = convergence criterion (default = 1e-6)
%   'startval'     = previously computed start-values
%   'dispprefix'   = prefix added to all disp-calls, handy when function is used in many loops after each other
%   'optimmode'    = Optimization mode of Sigma; 'singlecomppairals' or 'singlecomp' (default is set depending on Dmode), 
%                    optimize using a 1D grid search and Newton-Raphson/Gradient descent
%   'precision'    = number, precision used in many steps in the algorithm, but most importantly the smallest sigma difference
%   'degencrit'    = number, critical value at which to regard correlations between components as too high
%   'Dmode'        = 'identity', 'kdepcomplex', type of D to estimate/use
%   'holdparam'    = 1x5 vector of 0s and 1s indicating whether certain parameter sets are not updated in each ALS-iteration
%   'randomseed'   = scalar, seed to use for Matlab's random number generator
%
%


%
% Copyright (C) 2012-present, Roemer van der Meij, roemervandermeij AT gmail DOT com
%
% This file is part of Nwaydecomp, see https://github.com/roemervandermeij/nwaydecomp
%
%    Nwaydecomp is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    Nwaydecomp is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Nwaydecomp. If not, see <http://www.gnu.org/licenses/>.

%  TO DO: precision should really be split up in tolerance and (sigma) precision, keeping in mind their dependence
%  TO DO: precision conditions should be defined relative to the data, currently a workaround is implemented
%  TO DO: merge some of the subfunctions of SPACE models into externals
%  TO DO: the main stop conditions could be more principled wrt ssqres
%  TO DO: adjust the model formulation so conj(dat) is no longer necessary

% start execution timer
stopwatch = tic;

% Get the optional input arguments
keyvalcheck(varargin, 'optional', {'niter','convcrit','startval','dispprefix','precision','optimmode','degencrit','Dmode','holdparam','randomseed'});
niter        = keyval('niter', varargin);         if isempty(niter),        niter         = 2500;                   end
convcrit     = keyval('convcrit', varargin);      if isempty(convcrit),     convcrit      = 1e-8;                   end
startval     = keyval('startval', varargin); 
dispprefix   = keyval('dispprefix', varargin); 
precision    = keyval('precision', varargin);     if isempty(precision),    precision     = eps*1e7;                end
randomseed   = keyval('randomseed', varargin);    if isempty(randomseed),   randomseed    = round(sum(clock.*1e6)); end
degencrit    = keyval('degencrit', varargin);     if isempty(degencrit),    degencrit     = .9;                     end
Dmode        = keyval('Dmode', varargin);         if isempty(Dmode),        Dmode         = 'identity';             end
holdparam    = keyval('holdparam', varargin);     if isempty(holdparam),    holdparam     = zeros(1,5);             end
holdparam    = logical(holdparam);
optimmode    = keyval('optimmode', varargin);
if isempty(optimmode) && strcmp(Dmode,'kdepcomplex')
  optimmode = 'singlecomppairals';
elseif isempty(optimmode) && strcmp(Dmode,'identity')
  optimmode = 'singlecomp';
end
if strcmp(optimmode,'singlecomp') && ~strcmp(Dmode,'identity')
  error('optimmode = ''singlecomp'' only posssible when Dmode = ''identity''')
end

% load data from disk if string is given
if ischar(dat)
  filevars = whos('-file', dat);
  if numel(filevars)>1 || ~any(strcmp(filevars.class,{'single','double'}))
    error('filename contains either a single variable or not an array of singles or doubles')
  end
  datvarname = filevars.name;
  filecontent = load(dat);
  dat = filecontent.(datvarname);
  clear filecontent
end

% Dmode check
switch Dmode
  case 'identity'
  case 'kdepcomplex'
  otherwise
    error('improper Dmode')
end

% holdparam check
if any(size(holdparam)~=[1 5])
  error('improper ''holdparam''')
end

% set random seed (using clock as default)
rng(randomseed)

% set size and number of modes and compflg (and default compmodes)
nmode   = ndims(dat);
smode   = size(dat);
smodey  = [smode([1 2 3]) ncomp];
compflg = ~isreal(dat);

% display input data
dimstring = num2str(smode(1));
for imode = 2:nmode
  dimstring = [dimstring 'x' num2str(smode(imode))];
end
disp([dispprefix 'data is complex array with dimensions ' dimstring])
disp([dispprefix 'a SPACE-TIME model with ' num2str(ncomp) ' components will be estimated '])
% disp([dispprefix 'number of data-points per (relevant) parameter = ' num2str(prod(smodey) ./ (sum([smodey([1 2 3]) smodey(1) smodey(2)*smodey(4)])*ncomp),'%-2.1f')])
% disp([dispprefix 'number of data-points per (relevant) parameter = ' num2str((((smodey(1).^2)/2)*prod(smodey([2 3]))) ./ (sum([smodey([1 2 3]) smodey(1) smodey(2)*smodey(4)])*ncomp),'%-2.1f')])
disp([dispprefix 'number of data-points per (relevant) parameter = ' num2str(prod(smode) ./ (sum([smodey([1 2 3]) smodey(1) smodey(2)*smode(3)*smode(4)])*ncomp),'%-2.1f')])
disp([dispprefix 'maximum number of iterations = ' num2str(niter)])
disp([dispprefix 'convergence criterion = ' num2str(convcrit)])
disp([dispprefix 'precision used = ' num2str(precision)])
disp([dispprefix 'random seed = ' num2str(randomseed)])

% throw errors for inproper input
if ~compflg || (nmode ~= 4)
  error('input needs to be complex 4-way chan_time_freq_tap array')
end

% throw errors for inproper freq
if smode(2) ~= numel(F)
  error('F has improper number of elements')
end

% workaround for precision coditions not being defined relative to the data
% compute a running average of the exponent of the data
datexpfac = 0;
for ik = 1:smode(2)
  for il = 1:smode(3)
    % select dat
    currdatQ = double(permute(dat(:,ik,il,:),[1 4 2 3]));  % ensure double precision
    currdatQ(:,isnan(currdatQ(1,:))) = []; 
    csd = currdatQ*currdatQ';
    % get biggest negative exponent
    datexpfac = datexpfac + (mean(log10(abs(csd(csd~=0)))) ./ prod(smode([2 3])));
  end
end
datexpfac = min([floor(datexpfac./2) 0]);

% prepare datforQ, and use eigdecomp if applicable
datforQ    = cell(smode([2 3]));
eigflg     = false;
ssqdat     = 0;
zeropadflg = false;
for ik = 1:smode(2)
  for il = 1:smode(3)
    % select dat
    currdatQ = permute(dat(:,ik,il,:),[1 4 2 3]);
    currdatQ(:,isnan(currdatQ(1,:))) = [];
    
    % check whether speedup via eig can be done
    if size(currdatQ,1) < (size(currdatQ,2))
      % Reduce memory load and computation time by replacing each chan_tap matrix by the
      % Eigenvectors of its chan_chan cross-products weighted by sqrt(Eigenvalues).
      % This is possible because (1) SPACE only uses the cross-products of the chan_taper matrices
      % (i.e. the frequency- and trial-specific CSD) and (2) the Eigendecomposition of a symmetric
      % matrix A is A = VLV'.
      % As such, VL^.5 has the same cross-products as the original chan_tap matrix.
      currdatQ  = double(currdatQ); % ensure double precision
      csd       = currdatQ*currdatQ';
      [V L]     = eig(csd);
      L         = diag(L);
      tol       = max(size(csd))*eps(max(L)); % compute tol using matlabs default
      zeroL     = L<tol;
      eigweigth = V(:,~zeroL)*diag(sqrt(L(~zeroL)));
      eigweigth = single(eigweigth); % ensure single precision
      eigflg    = true;
      % save as currdatQ
      currdatQ = eigweigth;
    end
    
    % apply workaround for precision
    currdatQ = currdatQ ./ 10^datexpfac;
    
    % zero pad taper dimension
    if size(currdatQ,2)<ncomp
      zeropadflg = true;
      currdatQ = cat(2,currdatQ, zeros(size(currdatQ,1),ncomp-size(currdatQ,2)));
    end
    
    % take explicit conjugate, because the model expresses time-delays in a column vector (the shortest summary possible of a complicated story)
    currdatQ = conj(currdatQ);
    
    % save currdatforQ
    datforQ{ik,il} = currdatQ;
    % calc ssqdat as running sum
    ssqdat = ssqdat + sum(abs(double(currdatQ(:))).^2);
  end
end
if eigflg
  disp([dispprefix 'number of tapers higher than channels for some frequencies: using eigendecomposition to accelerate'])
end
if zeropadflg
  disp([dispprefix 'number of tapers is lower than number of components for one or more frequencies/trials, adding zero columns as filling tapers'])
end


% ensure F is a column vector, and throw error when not all integers
F = F(:);
if any(rem((round(F*10)),1)>eps*1e6)
  warning('frequencies are not integers, circularity point can become extremely large')
end
% find circularity point of sigma
% circulairy point is the smallest number number for which all frequencies have an integer number of cycles
% this number is equal to 1 over the greatest common divisor of vector F
% the greatest common divisor (GCD) of F is found by finding the GCD(GCD(F1,F2),F3), etc
% this is only defined for integers, so first scale up and and round
sigmacirc = round(F(1) .* 1e6);
for ik = 2:smode(2)
  sigmacirc = gcd(sigmacirc,round(F(ik) .* 1e6));
end
sigmacirc = 1/(sigmacirc ./ 1e6);
disp([dispprefix 'sigma circularity point = ' num2str(sigmacirc,'%-1.3f') 's'])

% set degenitercount
degenitercount = NaN;

% produce random start values if no start-values were provided
if isempty(startval)
  
  % use random start values
  comp = cell(1,nmode+1);
  % mode 1 'chan' A
  comp{1} = rand([smode(1) ncomp]);
  % mode 2 'freq' B
  comp{2} = rand([smode(2) ncomp]);
  % mode 3 'time' C
  comp{3} = rand([smode(3) ncomp]);
  % time-constants Sigma
  comp{4} = rand([smode(1) ncomp]) .* sigmacirc;
  % coherence matrix phi
  switch Dmode
    case 'identity'
      comp{5} = eye(ncomp);
    case 'kdepcomplex'
      comp{5} = complex((rand(smode(2),ncomp,ncomp)*2)-1,(rand(smode(2),ncomp,ncomp)*2)-1);
  end
  
  % norm A B and phi
  comp = normalizecomp(comp,'normA',[],smodey,[],[],0,dispprefix,Dmode);
  comp = normalizecomp(comp,'normB',[],smodey,[],[],0,dispprefix,Dmode);
  switch Dmode
    case 'identity'
      % nothing necessary
    case 'kdepcomplex'
      comp = normalizecomp(comp,'normD',[],smodey,[],[],0,dispprefix,Dmode);
  end
  
  % save startval
  startval = comp;
  
else
  comp = startval;
  
  % use random start values for empty cells
  % mode 1 'chan' A
  if isempty(comp{1})
    comp{1} = rand([smode(1) ncomp]);
  end
  % mode 2 'freq' B
  if isempty(comp{2})
    comp{2} = rand([smode(2) ncomp]);
  end
  % mode 3 'time' C
  if isempty(comp{3})
    comp{3} = rand([smode(3) ncomp]);
  end
  % time-constants Sigma
  if isempty(comp{4})
    comp{4} = rand([smode(1) ncomp]) .* sigmacirc;
  end
  % coherence matrix phi
  if isempty(comp{5})
    switch Dmode
      case 'identity'
        comp{5} = eye(ncomp);
      case 'kdepcomplex'
        comp{5} = complex((rand(smode(2),ncomp,ncomp)*2)-1,(rand(smode(2),ncomp,ncomp)*2)-1);
    end
  end
  
  % norm A B and phi
  comp = normalizecomp(comp,'normA',[],smodey,[],[],0,dispprefix,Dmode);
  comp = normalizecomp(comp,'normB',[],smodey,[],[],0,dispprefix,Dmode);
  switch Dmode
    case 'identity'
      % nothing necessary
    case 'kdepcomplex'
      comp = normalizecomp(comp,'normD',[],smodey,[],[],0,dispprefix,Dmode);
  end
  
  % save startval
  startval = comp;
  
end

% clear dat for memory issues
clear dat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   SPACE-time ALS   START  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The below code calculates an ALS estimate of the SPACE-time model
%
%
%
%
%
%


% set some  important variables
ssqres     = ssqdat;
prevssqres = 2*ssqres;
iter       = 0;

% create indices for subfunction calcmagmodely
mmiBind = repmat((1:smodey(2)).',[prod(smodey([3 4])) 1]);
mmiCind = repmat(reshape(repmat(1:smodey(3),[prod(smodey(2)) 1]),[prod(smodey([2 3])) 1]),[smodey(4) 1]);
mmiDind = mmiBind +  (reshape(repmat(1:smodey(4),[prod(smodey([2 3])) 1]),[prod(smodey([2 3 4])) 1])-1) .* smodey(2);

% start main loop of algorithm
disp([dispprefix 'starting ALS algorithm'])
while (abs((ssqres - prevssqres) / prevssqres) > convcrit) && (iter < niter) && (ssqres > precision) && (abs(ssqres - prevssqres) > precision) && (abs(ssqres ./ ssqdat) > precision)
  
  
  % Count iter
  iter = iter + 1;
  
  
  % Perform linear search every X iterations if relative ssqres increase is smaller than 10% (from the perspective of the previous iteration)
  if rem((iter-1),3) == 0 && ((prevssqres / ssqres) <= 1.10) && (((iter-1)^1/3) >= 2) && (iter < (niter-2))
    [comp,ssqres] = linsearch(datforQ,comp,prevcomp,F,smode,iter-1,prevssqres,ssqres,dispprefix,Dmode,holdparam,sigmacirc); % subfunction for performing linear search
    % normalize all signs, norms
    comp = normalizecomp(comp,'signA',[],smodey,[],[],0,dispprefix,Dmode);
    comp = normalizecomp(comp,'signB',[],smodey,[],[],0,dispprefix,Dmode);
    comp = normalizecomp(comp,'signC',[],smodey,[],[],0,dispprefix,Dmode);
    [comp,Pkl] = normalizecomp(comp,'allsignBC',[],smodey,[],Pkl,0,dispprefix,Dmode);
    comp = normalizecomp(comp,'normA',[],smodey,[],[],0,dispprefix,Dmode);
    comp = normalizecomp(comp,'normB',[],smodey,[],[],0,dispprefix,Dmode);
  end
  
  % Set previous stuff (important for e.g. linear search in next iteration)
  prevcomp = comp;
  prevssqres  = ssqres;
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate P by using Q, and calculate Y
  Pkl = calcPkl(comp,F,smode,datforQ,Dmode);
  
  % calculate Y
  Ydat = calcYdat(datforQ,Pkl,smodey);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  
  if holdparam(4)==0
    %%%%%%%%%%%%%
    % Update time-constant Sigma
    S = comp{4};
    % set ingredients for loss function evaluation
    % sigma circulairy point is defined above, outside of main loop
    % find largest active frequency
    Fmax = F(find(sum(abs(comp{2})>eps*1e6,2),1,'last'));
    estsigmax  = sigmacirc;
    estsigstep = .0625/Fmax;
    estsigall  = estsigstep:estsigstep:(estsigmax-estsigstep);
    estsigall(estsigall==0) = [];
    nestsig    = numel(estsigall);
    % set angleins
    angleinlfprt1 = bsxfun(@times,estsigall*(2*pi),F);
    % set up optimization bookkeeping
    sigupdfail  = 0;
    sigupdsame  = 0;
    sigupdsuc   = zeros(smode(1),ncomp);
    cpnewtloops = zeros(1,smode(1));
    cplsloops   = zeros(1,smode(1));
    for ij = 1:smode(1)
      % get current sigma
      sigma = S(ij,:);
      
      
      %%%%%%%%%%
      % Calculate constant ingredients for updating sigma
      % select ZdatY
      currdatY   = permute(Ydat(ij,:,:,:),[2 3 4 1]);
      ZdatYabs   = abs(currdatY);
      ZdatYangle = angle(currdatY);
      % compute model with magnitudes only
      magmodel = cell(1,ncomp);
      for icomp = 1:ncomp
        magmodel{icomp} = calcmagmodely(comp,smodey,icomp,Dmode,ij,mmiBind,mmiCind,mmiDind);
      end
      
      % calculate alpha, beta, phi
      alpha = NaN(smode(2),ncomp);
      beta  = NaN(smode(2),ncomp);
      theta = NaN(smode(2),ncomp);
      switch Dmode
        case 'identity'
          anglein = -ZdatYangle + pi; % first part the same for every channel
          cosin   = cos(anglein);
          sinin   = sin(anglein);
        case 'kdepcomplex'
          % computation done inside comp loop
      end
      for icomp = 1:ncomp
        % compute model with magnitudes only
        Zmodel  = magmodel{icomp};
        switch Dmode
          case 'identity'
            % computation done outside comp loop
          case 'kdepcomplex'
            anglein = bsxfun(@minus,reshape(angle(comp{5}(:,:,icomp)),[smode(2) 1 ncomp]), ZdatYangle) + pi; % first part the same for every channel
            cosin   = cos(anglein);
            sinin   = sin(anglein);
        end
        weight  = 2 .* ZdatYabs .* Zmodel; % first part the same for every main iteration
        alpha(:,icomp) = sum(sum(Zmodel.^2,3),2); % first part the same for every main iteration
        beta(:,icomp)  = sqrt(sum(sum(weight .* cosin,3),2).^2 + sum(sum(weight .* sinin,3),2).^2);
        theta(:,icomp) = atan2(sum(sum(weight .* sinin,3),2),sum(sum(weight .* cosin,3),2));
      end
      alpha = sum(sum(sum(ZdatYabs .^2,3),2) + sum(alpha,2),1);
      
      % calculate gamma and eta
      if ncomp ~= 1
        ncpair   = (ncomp*(ncomp-1))/2;
        cpairind = [];
        for icomp = 1:ncomp-1,
          cpairind = [cpairind; repmat(icomp,[ncomp-icomp,1]), ((icomp+1):ncomp)'];
        end
        gamma = zeros(smode(2),ncpair);
        eta   = zeros(smode(2),ncpair);
        switch Dmode
          case 'identity'
            % interaction term is always zero, so no need to calculate gamma/eta
          case 'kdepcomplex'
            for icpair = 1:ncpair
              % select component pair
              c1 = cpairind(icpair,1);
              c2 = cpairind(icpair,2);
              Zmodelc1 = magmodel{c1};
              Zmodelc2 = magmodel{c2};
              % anglein  = reshape(angle(comp{5}(:,c1) ./ comp{5}(:,c2)),[1 1 ncomp]); % first part the same for every channel (left here for possible 'nonkdepcomplex'
              anglein  = reshape(angle(comp{5}(:,:,c1) ./ comp{5}(:,:,c2)),[smode(2) 1 ncomp]); % first part the same for every channel
              cosin    = cos(anglein);
              sinin    = sin(anglein);
              weight   = 2 .* Zmodelc1 .* Zmodelc2;
              gamma(:,icpair) = sqrt(sum(sum(bsxfun(@times, weight, cosin),3),2).^2 + sum(sum(bsxfun(@times, weight, sinin),3),2).^2);
              eta(:,icpair)   = atan2(sum(sum(bsxfun(@times, weight, sinin),3),2),sum(sum(bsxfun(@times, weight, cosin),3),2));
            end
            % set NaNs in gamma and eta to zero
            gamma(isnan(gamma)) = 0;
            eta(isnan(eta))     = 0;
            % sort component pairs based on their highest interaction term, gamma
            [dum sortind] = sort(sum(gamma),'descend');
            gamma = gamma(:,sortind);
            eta   = eta(:,sortind);
            cpairind = cpairind(sortind,:);
        end
        %%%%%%%%%%
        
        
        % find sigma component vector by updating sigma alternating over component pairs
        % step 1, perform 2D grid search
        % step 2, find 2-component minumum using newton-raphson, starting from result of step 1
        % form first part of the loss function by expanding, for each component
        complfpart1 = cell(1,ncomp);
        for icomp = 1:ncomp
          complfpart1{icomp} = sum(bsxfun(@times,beta(:,icomp),cos(bsxfun(@plus,angleinlfprt1,theta(:,icomp)))),1);
        end
        switch optimmode
          
          case 'singlecomp'
            % set up optimization bookkeeping
            newtloopspair = zeros(1,ncomp);
            lsloopspair   = zeros(1,ncomp);
            
            % below approach is essentially the same as the approach for ncomp == 1 down below
            for icomp = 1:ncomp
              % consider the other parameters as constants
              constind = 1:ncomp;
              constind(icomp) = [];
              % get current lossfunval
              oldlossfunval = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
              
              
              % step 1: perform 1D grid search
              % calculate loss function ingredients dependent on active component
              lfpart1 = complfpart1{icomp};
              % all interaction terms are zero
              
              % calculate loss function ingredients dependent only on inactive components (i.e. constants)
              lfpartconst = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,cpairind,constind,Dmode);
              % sum above to final 1DGS loss function
              lf1D = lfpart1+lfpartconst;
              % find closest guess of sigma
              [dum minind] = min(lf1D);
              sigmanew = sigma;
              sigmanew(icomp) = estsigall(minind);
              
              % step 2: use gradient descent and newton-raphson to find 1D minimum
              % set controls for optimization
              maxoptimiter  = 50;
              [sigmanew,newtcount,lscount] = optimizesigma(icomp,sigmanew,alpha,beta,theta,gamma,eta,F,sigmacirc,cpairind,precision,maxoptimiter,Dmode);
              % save counts
              newtloopspair(icomp) = newtcount;
              lsloopspair(icomp)   = lscount;
              
              % only use new estimate when it decreases lossfun
              newlossfunval = sigmalossfun(sigmanew,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
              if newlossfunval < oldlossfunval
                sigma = sigmanew;
                sigupdsuc(ij,icomp) = sigupdsuc(ij,icomp) + 1;
              elseif (((newlossfunval-oldlossfunval)/oldlossfunval)<eps*1e6) && all(~(abs(sigmanew-sigma)>precision*10))
                sigupdsame = sigupdsame + 1;
              else
                sigupdfail = sigupdfail + 1;
              end
              
            end % icomp
            % process counts
            cpnewtloops(ij) = mean(newtloopspair);
            cplsloops(ij)   = mean(lsloopspair);
            
          case 'singlecomppairals'
            % The code below I left for legacy purposes. It is the old implementation of 'singlecomp'. I.e., it's the same procedure as below,
            % but only touches a single component once. This can still be a speedy strategy if Dmode = 'kdepcomplex', with small interaction terms
            % isupdated = false(1,ncomp);
            % for icpair = 1:ncpair
            %   for icpaircomp = 1:2
            %     currcomp = cpairind(icpair,icpaircomp);
            %     % check whether component has already been updated
            %     if isupdated(currcomp)
            %       continue
            %     end
            % % flag component as updated
            % isupdated(icomp) = true;
            
            % set up optimization bookkeeping
            newtloopspair = zeros(ncpair,2);
            lsloopspair   = zeros(ncpair,2);
            for icpair = 1:ncpair
              for icpaircomp = 1:2
                currcomp = cpairind(icpair,icpaircomp);
                % consider the other parameters as constants
                constind = 1:ncomp;
                constind(currcomp) = [];
                % get current lossfunval
                oldlossfunval = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
                
                
                % step 1: perform 1D grid search
                % calculate loss function ingredients dependent on active component
                lfpart1 = complfpart1{currcomp};
                if icpaircomp == 1
                  lfpart2 = sum(bsxfun(@times,gamma(:,icpair),cos(bsxfun(@plus,bsxfun(@minus,angleinlfprt1,2*pi*F*sigma(cpairind(icpair,2))),eta(:,icpair)))),1);
                elseif icpaircomp == 2
                  lfpart2 = sum(bsxfun(@times,gamma(:,icpair),cos(bsxfun(@plus,bsxfun(@minus,2*pi*F*sigma(cpairind(icpair,1)),angleinlfprt1),eta(:,icpair)))),1);
                end
                % calculate loss function ingredients dependent on active components
                partactind = find(sum(cpairind(:,1)==currcomp,2));
                partactind(partactind==icpair) = [];
                lfpart3 = zeros(1,nestsig);
                %sum(sum(bsxfun(@times,reshape(gamma(:,partactind).',[numel(partactind) 1 smode(2)]),cos(bsxfun(@plus,bsxfun(@minus,reshape(angleinlfprt1.',[1 nestsig smode(2)]),reshape((2*pi*F*sigma(cpairind(partactind,2))).',[numel(partactind) 1 smode(2)])),reshape(eta(:,partactind).',[numel(partactind) 1 smode(2)])))),3),1);
                if ~isempty(partactind)
                  for ipartactind = 1:numel(partactind)
                    currind = partactind(ipartactind);
                    lfpart3 = lfpart3 + sum(bsxfun(@times,gamma(:,currind),cos(bsxfun(@plus,bsxfun(@minus,angleinlfprt1,2*pi*F*sigma(cpairind(currind,2))),eta(:,currind)))));
                  end
                end
                partactind = find(sum(cpairind(:,2)==currcomp,2));
                partactind(partactind==icpair) = [];
                if ~isempty(partactind)
                  for ipartactind = 1:numel(partactind)
                    currind = partactind(ipartactind);
                    lfpart3 = lfpart3 + sum(bsxfun(@times,gamma(:,currind),cos(bsxfun(@plus,bsxfun(@minus,2*pi*F*sigma(cpairind(currind,1)),angleinlfprt1),eta(:,currind)))));
                  end
                end
                
                % calculate loss function ingredients dependent only on inactive components (i.e. constants)
                lfpartconst = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,cpairind,constind,Dmode);
                % sum above to final 1DGS loss function
                lf1D = lfpart1+lfpart2+lfpart3+lfpartconst;
                % find closest guess of sigma
                [dum minind] = min(lf1D);
                sigmanew = sigma;
                sigmanew(currcomp) = estsigall(minind);
                
                % step 2: use gradient descent and newton-raphson to find 1D minimum
                % set controls for optimization
                maxoptimiter  = 50;
                [sigmanew,newtcount,lscount] = optimizesigma(currcomp,sigmanew,alpha,beta,theta,gamma,eta,F,sigmacirc,cpairind,precision,maxoptimiter,Dmode);
                % save counts
                newtloopspair(icpair,icpaircomp) = newtcount;
                lsloopspair(icpair,icpaircomp)   = lscount;
                
                % only use new estimate when it decreases lossfun
                newlossfunval = sigmalossfun(sigmanew,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
                if newlossfunval < oldlossfunval
                  sigma = sigmanew;
                  sigupdsuc(ij,currcomp) = sigupdsuc(ij,currcomp) + 1;
                elseif (((newlossfunval-oldlossfunval)/oldlossfunval)<eps*1e6) && all(~(abs(sigmanew-sigma)>precision*10))
                  sigupdsame = sigupdsame + 1;
                else
                  sigupdfail = sigupdfail + 1;
                end
                
              end % icpaircomp
            end % icpair
            % process counts
            cpnewtloops(ij) = mean(mean(newtloopspair));
            cplsloops(ij)   = mean(mean(lsloopspair));
            
          otherwise
            error('invalid optimmode')
        end
        
        % save sigma
        S(ij,:) = sigma;
        
      elseif ncomp == 1
        % find sigma
        % step 1, perform 1D grid search
        % step 2, find minumum using newton-raphson, starting from result of step 1
        % get current lossfunval
        oldlossfunval = sigmalossfun(sigma,alpha,beta,theta,[],[],F,[],1,Dmode);
        
        % perform 1D grid search
        % compute loss function
        lf1D = alpha + sum(bsxfun(@times,beta,cos(bsxfun(@plus,angleinlfprt1,theta))),1);
        % find closest guess of sigma
        [dum minind] = min(lf1D);
        sigmanew = estsigall(minind);
        
        % step 2: use gradient descent and newton-raphson to find 1D minimum
        % update best guesses for all sigma by another optimization run
        % set controls for optimization
        maxoptimiter  = 50;
        [sigmanew,newtcount,lscount] = optimizesigma(1,sigmanew,alpha,beta,theta,[],[],F,sigmacirc,[ ],precision,maxoptimiter,Dmode);
        % save counts
        cpnewtloops(ij) = newtcount;
        cplsloops(ij)   = lscount;
        
        % only use new estimate when it decreases lossfun
        newlossfunval = sigmalossfun(sigmanew,alpha,beta,theta,[],[],F,[],1,Dmode);
        if newlossfunval < oldlossfunval
          sigma = sigmanew;
          sigupdsuc(ij) = sigupdsuc(ij) + 1;
        elseif (((newlossfunval-oldlossfunval)/oldlossfunval)<eps*1e6) && ~(abs(sigmanew-sigma)>precision*10)
          sigupdsame = sigupdsame + 1;
        else
          sigupdfail = sigupdfail + 1;
        end
        
        % save sigma
        S(ij) = sigma;
      end
    end % ij
    % save S
    comp{4} = S;
    %%%%%%%%%%%%%
  end
  
  
  
  
  
  if holdparam(1)==0
    %%%%%%%%%%%%%
    % Update A
    % by doing regular PARAFAC ALS, but loop over J
    A = zeros(smode(1),ncomp);
    for ij = 1:smode(1)
      
      % Calculate Z
      currB = comp{2};
      currC = comp{3};
      currS = exp(1i*2*pi*repmat(F,[1 ncomp]) .* repmat(comp{4}(ij,:),[smode(2) 1]));
      switch Dmode
        case 'identity'
          currD = comp{5};
          Z = kr(currD,currC,currB.*currS);
        case 'kdepcomplex'
          currD = reshape(permute(comp{5},[3 1 2]),[ncomp prod(smodey([2 4]))]).';
          Bind  = repmat((1:smodey(2)).',[prod(smodey([3 4])) 1]);
          Dind  = Bind +  (reshape(repmat(1:smodey(4),[prod(smodey([2 3])) 1]),[prod(smodey([2 3 4])) 1])-1) .* smodey(2);
          Z = repmat(kr(currC,currB.*currS),[smodey(4) 1]) .* currD(Dind,:);
      end
      
      % make currdat
      permorder   = [1 setdiff(1:4,1)];
      reshapesize = [1 prod(smodey([2 3 4]))];
      Zdat = reshape(permute(Ydat(ij,:,:,:),permorder),reshapesize);
      
      % calculate currA
      %A(ij,:) = fastnnls(real(Z'*Z),real(Zdat * conj(Z)).');
      A(ij,:) = real(Zdat * conj(Z)) * pinv(real(Z'*Z)).';
    end
    % save A
    comp{1} = A;
    % normalize signA and normA
    comp = normalizecomp(comp,'signA',[],smodey,[],[],0,dispprefix,Dmode);
    comp = normalizecomp(comp,'normA',[],smodey,[],[],0,dispprefix,Dmode);
    %%%%%%%%%%%%%
  end
  
  
  
  
  
  if holdparam(5)==0
    %%%%%%%%%%%%%
    % Update D
    switch Dmode
      case 'identity'
        % nothing necessary
      case 'kdepcomplex'
        % by doing regular PARAFAC ALS, but loop over K
        D = zeros(smodey(2),ncomp,ncomp);
        for ik = 1:smodey(2)
          
          % Calculate Z
          currA = comp{1};
          currB = comp{2}(ik,:);
          currC = comp{3};
          currS = comp{4};
          Z = kr(currC,currB,currA .* exp(1i*2*pi*F(ik).*currS));
          
          % make currdat
          permorder   = [4 setdiff(1:4,4)];
          reshapesize = [smodey(4) smodey(1)*1*smodey(3)];
          Zdat = reshape(permute(Ydat(:,ik,:,:),permorder),reshapesize);
          
          % calculate D
          D(ik,:,:) = Zdat * conj(Z) * pinv(Z'*Z).';
          %D(ik,:,:) = (Zdat * conj(Z)) / (Z'*Z).';
        end
        % save D
        comp{5} = D;
        % normalize D
        comp = normalizecomp(comp,'normD',[],smodey,[],[],0,dispprefix,Dmode);
    end
    %%%%%%%%%%%%%
  end
  
  
  
  
  
  if holdparam(2)==0
    %%%%%%%%%%%%%
    % Update B
    % by doing regular PARAFAC ALS, but loop over K
    B = zeros(smode(2),ncomp);
    for ik = 1:smode(2)
      
      % Calculate Z
      currA = comp{1} .* exp(1i*2*pi*F(ik).*comp{4});
      currC = comp{3};
      switch Dmode
        case 'identity'
          currD = comp{5};
        case 'kdepcomplex'
          currD = squeeze(comp{5}(ik,:,:));
      end
      Z = kr(currD,currC,currA);
      
      % make currdat
      permorder   = [2 setdiff(1:4,2)];
      reshapesize = [1 prod(smodey([1 3 4]))];
      Zdat = reshape(permute(Ydat(:,ik,:,:),permorder),reshapesize);
      
      % calculate currB
      %B(ik,:) = fastnnls(real(Z'*Z),real(Zdat * conj(Z)).');
      B(ik,:) = real(Zdat * conj(Z)) * pinv(real(Z'*Z)).';
    end
    % save B
    comp{2} = B;
    % normalize signB and normB
    comp = normalizecomp(comp,'signB',[],smodey,[],[],0,dispprefix,Dmode);
    comp = normalizecomp(comp,'normB',[],smodey,[],[],0,dispprefix,Dmode);
    %%%%%%%%%%%%%
  end
  
  
  
  
  if holdparam(3)==0
    %%%%%%%%%%%%%
    % Update C
    % by doing regular PARAFAC ALS
    
    % move norm over to A
    comp = normalizecomp(comp,'normC',[],smodey,[],[],0,dispprefix,Dmode);
    
    % Calculate Z
    currA = comp{1};
    currB = comp{2};
    Find  = reshape(repmat(1:smodey(2),[smodey(1) 1]),[1 prod(smodey([1 2]))]);
    Sind  = repmat(1:smodey(1),[1 smodey(2)]);
    currS = exp(1i*2*pi*repmat(F(Find),[1 ncomp]).*comp{4}(Sind,:));
    switch Dmode
      case 'identity'
        currD = comp{5};
        Z = kr(currD,kr(currB,currA) .* currS);
      case 'kdepcomplex'
        currD = reshape(permute(comp{5},[3 1 2]),[ncomp prod(smodey([2 4]))]).';
        Bind  = repmat(reshape(repmat(1:smodey(2),[smodey(1) 1]),[prod(smodey([1 2])) 1]),[smodey(4) 1]);
        Dind  = Bind +  (reshape(repmat(1:smodey(4),[prod(smodey([1 2])) 1]),[prod(smodey([1 2 4])) 1])-1) .* smodey(2);
        Z = repmat(kr(currB,currA) .* currS,[smodey(4) 1]) .* currD(Dind,:);
    end
    
    % make currdat
    permorder   = [3 setdiff(1:4,3)];
    reshapesize = [smodey(3) prod(smodey([1 2 4]))];
    Zdat = reshape(permute(Ydat,permorder),reshapesize);
    
    % calculate currC
    %C = fastnnls(real(Z'*Z),real(Zdat * conj(Z)).');
    C = real(Zdat * conj(Z)) * pinv(real(Z'*Z)).';
    
    % clear Zdat and Z, as they have same size as Ydat
    clear Zdat Z
    
    % save C
    comp{3} = C;
    % normalize signC
    comp = normalizecomp(comp,'signC',[],smodey,[],[],0,dispprefix,Dmode);
    % put back norm in C
    comp = normalizecomp(comp,'normA',[],smodey,[],[],0,dispprefix,Dmode);
    %%%%%%%%%%%%%
  end
  
  
  
  
  
  
  %%%%%%%%%%%%%
  % normalize lambda's and all signs of B and C depending on Dmode
  [comp,Pkl] = normalizecomp(comp,'sigma',sigmacirc,smodey,F,Pkl,0,dispprefix,Dmode);
  [comp,Pkl] = normalizecomp(comp,'allsignBC',[],smodey,[],Pkl,0,dispprefix,Dmode);
  %%%%%%%%%%%%%
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute fit using a QR-orthogonalization shortcut
  % This is equivalent to the code below, but doesn't form the full model array. This trick depends on two things.
  % 1) sum(abs(dat-model).^2) == ssqdat-ssqmodel, and 2) using R of a QR decomposition of each loading matrix instead
  % 2) is easily verifiable: trace((A*B)*(A*B)) ==  trace(R1'R1 + tR2'R2). Number 1 I do not understand yet.
  % perform QR-decomposition of all matrices, except B because of sigma (the R of Pkl*Dk is constant over L)
  % right now, the accuracy of the trick seems to decay when getting closer to the end
  if false %all(smodey>=ncomp) && ~strcmp(Dmode,'identity') % doesn't work for Dmode = 'identity', don't understand why
    % permute Dk for faster acces
    ASkr = cell(1,smode(2));
    for ik = 1:smode(2)
      [dum ASkr{ik}] = qr(comp{1}.*exp(1i*2*pi*F(ik)*comp{4}),0);
    end
    [dum Cr] = qr(comp{3},0);
    PklDkr = cell(smode(2),ncomp);
    Dk = permute(comp{5},[2 3 1]);
    for ik = 1:smode(2)
      currDk = Dk(:,:,ik);
      for il = 1:ncomp
        [dum PklDkr{ik,il}] = qr((currDk.' * Pkl{ik,il}.').',0);
      end
    end
    % compute model
    rmodelx = zeros([ncomp smode(2) ncomp ncomp]);
    for ik = 1:smode(2)
      currASr = ASkr{ik};
      for il = 1:ncomp
        currm   = size(PklDkr{ik,il},1);
        currPklDkr = PklDkr{ik,il};
        rmodelx(:,ik,il,1:currm) = currASr * diag(Cr(il,:)) * diag(comp{2}(ik,:)) * currPklDkr.';
      end
    end
    ssqmodel = sum(abs(rmodelx(:)).^2);
    ssqres   = ssqdat - ssqmodel;
    expvar   = 100 - ((ssqres / ssqdat) * 100);
  else
    ssqres   = calcssqressparse(comp,Pkl,F,smode,datforQ,Dmode); % subfunction for calculating ssqres
    ssqmodel = ssqdat - ssqres;
    expvar   = 100 - ((ssqres / ssqdat) * 100);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Degeneracy failsafe
  % calculate tucker's congruency coefficient, indirectly
  tuckcongr = zeros(ncomp,ncomp);
  switch Dmode
    case 'identity'
      Cn = comp{3};
      Cn = Cn ./ repmat(sqrt(sum(abs(Cn).^2,1)),[smode(3) 1]);
      for ik = 1:smode(2)
        AS    = comp{1} .* exp(1i*2*pi*F(ik)*comp{4});
        AStAS = AS' * AS;
        BktBk = comp{2}(ik,:)' * comp{2}(ik,:);
        for il = 1:smode(3)
          CntCn = (Cn(il,:)' * Cn(il,:));
          PklDk       = conj(Pkl{ik,il}) * comp{5};
          PklDktPklDk = PklDk' * PklDk;
          tuckcongr = tuckcongr + AStAS .* BktBk .* CntCn .* PklDktPklDk;
        end
      end
    case 'kdepcomplex'
      Dk = permute(comp{5},[2 3 1]);
      Cn = comp{3};
      Cn = Cn ./ repmat(sqrt(sum(abs(Cn).^2,1)),[smode(3) 1]);
      for ik = 1:smode(2)
        AS    = comp{1} .* exp(1i*2*pi*F(ik)*comp{4});
        AStAS = AS' * AS;
        BktBk = comp{2}(ik,:)' * comp{2}(ik,:);
        currDk = Dk(:,:,ik);
        for il = 1:smode(3)
          CntCn = (Cn(il,:)' * Cn(il,:));
          PklDk       = conj(Pkl{ik,il}) * currDk;
          PklDktPklDk = PklDk' * PklDk;
          tuckcongr = tuckcongr + AStAS .* BktBk .* CntCn .* PklDktPklDk;
        end
      end
  end
  % remove upper triangle and diagonal, make vector, and take abs
  tuckcongr = tuckcongr(logical(tril(ones(ncomp,ncomp),-1)));
  tuckcongr = abs(tuckcongr);
  if ncomp==1 % if ncomp is one, tuckers congruence cannot be calculated
    tuckcongr = NaN;
  end
  % display warning if too high
  if max(tuckcongr) >= degencrit
    % start counting if counting hasn't started
    if isnan(degenitercount)
      disp([dispprefix 'warning: some components are correlated above critical value, model might be degenerate, stopping after 30 iterations if model hasn''t escaped'])
      degenitercount = 0;
    end
    degenitercount = degenitercount + 1;
  else
    degenitercount = NaN;
  end
  if degenitercount == 30
    disp([dispprefix 'model didn''t escape from possible degeneracy, algorithm stopped'])
    disp([dispprefix 'explained variance by model: ' num2str(expvar,'%-2.1f') '%'])
    break
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  % Display results of current iteration
  if holdparam(4) == 0
  disp([dispprefix 'iteration ' num2str(iter) ' - expvar: ' num2str(expvar,'%-2.1f')   '%  ssqres: ' num2str(ssqres)  '  ssqmodel: ' num2str(ssqmodel),...
    ' - NR-steps: ' num2str(mean(cpnewtloops),'%-2.1f') ' ('  num2str(std(cpnewtloops),'%-2.1f') ') | ', ...
       'LS-steps: ' num2str(mean(cplsloops),'%-2.1f')   ' ('  num2str(std(cplsloops),'%-2.1f')   ') | sigma-nu/uf/ni: ' num2str(sum(~sigupdsuc(:)),'%-4.0f')  '/' num2str(sigupdfail,'%-4.0f')  '/' num2str(sigupdsame,'%-4.0f') ' | max tuckcongr: ' num2str(max(tuckcongr),'%-2.3f') ])
  else
      disp([dispprefix 'iteration ' num2str(iter) ' - expvar: ' num2str(expvar,'%-2.1f')   '%  ssqres: ' num2str(ssqres)  '  ssqmodel: ' num2str(ssqmodel) ' | max tuckcongr: ' num2str(max(tuckcongr),'%-2.3f') ])
  end
  if prevssqres<ssqres && iter>1
    disp(['iteration ' num2str(iter) ' - warning: moved away from solution, ssqres increased by ' num2str(ssqres-prevssqres) ' and ' num2str(((ssqres-prevssqres)/prevssqres)*100) '%'])
    warning('Moved away from solution, something went horribly wrong. Restoring previous best estimate and stopping algorithm')
    comp = prevcomp;
    ssqres = prevssqres;
    Pkl  = calcPkl(comp,F,smode,datforQ,Dmode);
    break
  end
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Algorithm stops
  if ~(abs((ssqres - prevssqres) / prevssqres) > convcrit)
    disp([dispprefix 'convergence criterion of ' num2str(convcrit) ' reached in ' num2str(iter) ' iterations'])
    disp([dispprefix 'explained variance by model: ' num2str(expvar,'%-2.1f') '%'])
    % calculate final expvar without using qr-decomposition and throw error when deviating
    ssqresnonqr = calcssqressparse(comp,Pkl,F,smode,datforQ,Dmode); % subfunction for calculating ssqres
    if (abs(ssqresnonqr - ssqres) / ssqresnonqr) > (eps*1e6)
      warning(['problem with QR-based acceleration: ssqres diff =  ' num2str(ssqres-ssqresnonqr) ' and ' num2str(((ssqres-ssqresnonqr)/ssqresnonqr)*100) '%'])
    end
  elseif (iter == niter)
    disp([dispprefix 'maximum number of iterations = ' num2str(iter) ' reached'])
    disp([dispprefix 'explained variance by model: ' num2str(expvar,'%-2.1f') '%'])
    % following stop criteria only useful during simulations
  elseif ~(abs(ssqres - prevssqres) > precision)
    disp([dispprefix 'difference in residual sums of squares reached precision limit of ' num2str(precision) ' , algorithm stopped'])
    disp([dispprefix 'explained variance by model: ' num2str(expvar,'%-2.1f') '%'])
  elseif ~(abs(ssqres ./ ssqdat) > precision)
    disp([dispprefix 'fraction of residual sums of squares to sums of squares of the data reached precision limit of ' num2str(precision) ' , algorithm stopped'])
    disp([dispprefix 'explained variance by model: ' num2str(expvar,'%-2.1f') '%'])
  elseif ssqres<precision
    disp([dispprefix 'residual sums of squares reached precision, algorithm stopped'])
    disp([dispprefix 'explained variance by model: ' num2str(expvar,'%-2.1f') '%'])
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
end % end main while loop of algorithm (updating component matrices)
% stop conditions fail check
if iter==0
  error('stop conditions were not satisfied at iter = 0')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   SPACE-time ALS        END  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   POSTPROCESSING START  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize sigma's one last time
comp = normalizecomp(comp,'sigma',sigmacirc,smodey,F,Pkl,1,dispprefix,Dmode);

% Sort components based on ssq
compssqres = sum(abs(comp{3}).^2,1);
[dum, sortind] = sort(compssqres,'descend');
% sort all components
comp{1} = comp{1}(:,sortind);
comp{2} = comp{2}(:,sortind);
comp{3} = comp{3}(:,sortind);
comp{4} = comp{4}(:,sortind);
switch Dmode
  case 'identity'
    % do not sort the eye
  case 'kdepcomplex'
    comp{5} = comp{5}(:,:,sortind);
end
for ik = 1:smode(2)
  for il = 1:smode(3)
    Pkl{ik,il} = Pkl{ik,il}(:,sortind);
  end
end
disp([dispprefix 'components have been ordered according to magnitude'])


% Calculate Tucker3 core if desired (used for core consistency diagnostic)
if nargout > 6
  disp([dispprefix 'calculating tucker3 core array'])
  % build ZtZ
  ZtZ = zeros(ncomp^4); % 4 = number of modes
  for ik = 1:smode(2)
    AtA = (comp{1} .* exp(1i*2*pi*F(ik)*comp{4}))' * (comp{1} .* exp(1i*2*pi*F(ik)*comp{4}));
    BtB = comp{2}(ik,:)' * comp{2}(ik,:);
    for il = 1:smode(3)
      CtC = comp{3}(il,:)' * comp{3}(il,:);
      switch Dmode
        case 'identity'
          PD = conj(Pkl{ik,il}) * comp{5};
        case 'kdepcomplex'
          PD = conj(Pkl{ik,il}) * squeeze(comp{5}(ik,:,:));
      end
      PDtPD = PD' * PD;
      ZtZ = ZtZ + kron(PDtPD,kron(CtC,kron(BtB,AtA)));
    end
  end
  % build XtZ
  XtZ = zeros(ncomp^4,1); % 4 = number of modes
  ccind = 0;
  for imcomp = 1:ncomp
    for ilcomp = 1:ncomp
      for ikcomp = 1:ncomp
        for ijcomp = 1:ncomp
          ccind = ccind + 1;
          % calculate model for current component combination for each ik and il
          for ik = 1:smode(2)
            AS = comp{1}(:,ijcomp) .* exp(1i*2*pi*F(ik)*comp{4}(:,ijcomp));
            for il = 1:smode(3)
              switch Dmode
                case 'identity'
                  currPklDk = comp{5}.' * Pkl{ik,il}';
                case 'kdepcomplex'
                  currPklDk = squeeze(comp{5}(ik,:,:)).' * Pkl{ik,il}';
              end
              modelx = AS * diag(comp{3}(il,ilcomp)) * diag(comp{2}(ik,ikcomp)) * currPklDk(imcomp,:);
              % get data for current ik and il
              currdatQ = datforQ{ik,il};
              currdatQ = double(currdatQ); % ensure double precision
              % calc XtZ
              XtZ(ccind) = XtZ(ccind) + currdatQ(:)' * modelx(:);
            end
          end
        end
      end
    end
  end
  % calculate tucker3 core array
  t3core = pinv(real(ZtZ)) * real(XtZ);
end


% Magnitude post-processing of set of C
% set norm of C to 1 per component and save scaling coeffient
normC = sqrt(sum(abs(comp{3}).^2,1));
comp{3} = comp{3} ./ repmat(normC,[smode(3) 1]);
scaling = normC;
disp([dispprefix 'magnitude scaling coefficient of C removed'])


% Calculate Tucker's Congruency coefficient
% do it directly instead of indirectly as in parafac because of model structure
% calculate tucker's congruency coefficient, indirectly
tuckcongr = zeros(ncomp,ncomp);
for ik = 1:smode(2)
  AS    = comp{1} .* exp(1i*2*pi*F(ik)*comp{4});
  AStAS = AS' * AS;
  BktBk = comp{2}(ik,:)' * comp{2}(ik,:);
  switch Dmode
    case 'identity'
      currDk = comp{5};
    case 'kdepcomplex'
      currDk = squeeze(comp{5}(ik,:,:));
  end
  for il = 1:smode(3)
    CntCn = (comp{3}(il,:)' * comp{3}(il,:));
    PklDk       = conj(Pkl{ik,il}) * currDk;
    PklDktPklDk = PklDk' * PklDk;
    tuckcongr = tuckcongr + AStAS .* BktBk .* CntCn .* PklDktPklDk;
  end
end
% remove upper triangle and diagonal, make vector, and take abs
tuckcongr = tuckcongr(logical(tril(ones(ncomp,ncomp),-1)));
tuckcongr = abs(tuckcongr);
if ncomp==1 % if ncomp is one, tuckers congruence cannot be calculated
  tuckcongr = NaN;
end
% display warning if too high
if max(tuckcongr) >= degencrit
  disp([dispprefix 'Warning: some components are correlated above critical value, model might be degenerate'])
end

% apply precision work around to the scaling coefficients and ssqres
scaling = scaling .*  10.^datexpfac;
ssqres  = ssqres  .* (10.^datexpfac).^2; 

% report execution time
exectime = toc(stopwatch);
disp([dispprefix 'execution took ' num2str(exectime,'%.1f') ' seconds'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   POSTPROCESSING END    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      SUBFUNCTIONS START HERE     %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction normalizecomp       %%%%%%%%%%%%
function [comp,Pkl] = normalizecomp(comp,normmode,sigmacirc,smodey,F,Pkl,dispflg,dispprefix,Dmode)
% set stuff for subfunction internal use
ncomp = size(comp{1},2);


%%%%%%%%%%%%%
if strcmp(normmode,'sigma')
  % Normalize sigma's for unicity, and inversly normalize D or Pkl
  % normalize sigma such that the magnitude weighted mean per component is sigmacirc/2
  A = comp{1};
  S = comp{4};
  % set radian conversion factor
  radconv = ((2*pi)./sigmacirc);
  % convert to complex domain
  Scomp = exp(1i*S.*radconv);
  % get mangitude weighted mean and add sigmacirc/2
  magwmeanScomp = mean(Scomp.*A,1) ./ exp(1i*(sigmacirc/2)*radconv);
  % normalize sigma
  Scomp = Scomp ./ repmat(magwmeanScomp,[smodey(1) 1]);
  % move sigma back to the real domain
  S      = angle(Scomp);
  S(S<0) = S(S<0) + (2*pi);
  S      = S ./ radconv;
  % inversily normalize D in the real domain
  magwmeanScomp = angle(magwmeanScomp);
  magwmeanScomp(magwmeanScomp<0) = magwmeanScomp(magwmeanScomp<0) + (2*pi);
  magwmeanScomp = magwmeanScomp ./ radconv;
  switch Dmode
    case 'identity'
      for ik = 1:smodey(2)
        for il = 1:smodey(3)
          currPkl = Pkl{ik,il};
          currm   = size(currPkl,1);
          currPkl = currPkl .* conj(repmat(exp(1i*F(ik)*2*pi*magwmeanScomp),[currm 1])); % conj is necessary because of complex conjugate transpose in model formulation
          Pkl{ik,il} = currPkl;
        end
      end
    case 'kdepcomplex'
      D = comp{5};
      D = D .* exp(1i*permute(repmat(F*2*pi*magwmeanScomp,[1 1 smodey(4)]),[1 3 2]));
      comp{5} = D;
  end
  
  % save normalized coefficients
  comp{4} = S;
  if dispflg
    disp([dispprefix 'sigma has been shifted such that magnitude weighted average is equal to circularity-point/2'])
  end
end
%%%%%%%%%%%%%


%%%%%%%%%%%%%
if strcmp(normmode,'signA')
  % sign normalize A onto D
  A = comp{1};
  D = comp{5};
  % get mean sign, remove from A and apply to D
  meansignA = sign(mean(A,1));
  A = A .* repmat(meansignA,[smodey(1) 1]);
  switch Dmode
    case 'identity'
      D = D .* repmat(meansignA,[ncomp 1]);
    case 'kdepcomplex'
      D = D .* repmat(reshape(meansignA,[1 1 ncomp]), [smodey(2) ncomp 1]);
  end
  % save A and D
  comp{1} = A;
  comp{5} = D;
end
%%%%%%%%%%%%%


%%%%%%%%%%%%%
if strcmp(normmode,'signB')
  % sign normalize B onto D
  B = comp{2};
  D = comp{5};
  switch Dmode
    case 'identity'
      meansignB = sign(mean(B,1));
      B = B .* repmat(meansignB,[smodey(2) 1]);
      D = D .* repmat(meansignB,[ncomp 1]);
    case 'kdepcomplex'
      % get sign, remove from B and flip D 180 degrees (i.e. .* -1)
      signB = sign(B);
      B = B ./ signB;
      D = D .* repmat(reshape(signB,[smodey(2) 1 ncomp]),[1 ncomp 1]);
  end
  % save B and D
  comp{2} = B;
  comp{5} = D;
end
%%%%%%%%%%%%%


%%%%%%%%%%%%%
if strcmp(normmode,'signC')
  % sign normalize C onto D
  C = comp{3};
  D = comp{5};
  % get mean sign, remove from C and apply to D
  meansignC = sign(mean(C,1));
  C = C .* repmat(meansignC,[smodey(3) 1]);
  switch Dmode
    case 'identity'
      D = D .* repmat(meansignC, [ncomp 1]);
    case 'kdepcomplex'
      D = D .* repmat(reshape(meansignC,[1 1 ncomp]), [smodey(2) ncomp 1]);
  end
  % save C and D
  comp{3} = C;
  comp{5} = D;
end
%%%%%%%%%%%%%


%%%%%%%%%%%%%
if strcmp(normmode,'allsignBC')
  % sign normalize B and C onto Pkl, or build a Qkl
  B = comp{2};
  C = comp{3};
  D = permute(comp{5},[2 3 1]);
  % if Dmode is identity, get sign, remove from B and C and flip Pkl 180 degrees (i.e. .* -1)
  signB = sign(B);
  signC = sign(C);
  signB(signB==0) = 1;
  signC(signC==0) = 1;
  % apply normalization dependent on Dmode
  switch Dmode
    case 'identity'
      B = B ./ signB;
      C = C ./ signC;
      for ik = 1:smodey(2)
        for il = 1:smodey(3)
          if any(signB(ik,:)==-1) || any(signC(il,:)==-1)
            currPkl = Pkl{ik,il};
            newPkl = currPkl * diag(signB(ik,:)) * diag(signC(il,:));
            % save newPkl
            Pkl{ik,il} = newPkl;
          end
        end
      end
    case 'kdepcomplex'
      % not applicable
  end
  % save B and C (Pkl is directly modified)
  comp{2} = B;
  comp{3} = C;
end
%%%%%%%%%%%%%


%%%%%%%%%%%%%
if strcmp(normmode,'normA')
  % normalize A onto C
  A = comp{1};
  C = comp{3};
  % get norm, remove from A and apply to C
  normA = sqrt(sum(abs(A).^2,1));
  A = A ./ repmat(normA,[smodey(1) 1]);
  C = C .* repmat(normA,[smodey(3) 1]);
  % save A and C
  comp{1} = A;
  comp{3} = C;
end
%%%%%%%%%%%%%


%%%%%%%%%%%%%
if strcmp(normmode,'normB')
  % normalize B onto C
  B = comp{2};
  C = comp{3};
  % get norm, remove from B and apply to C
  normB = sqrt(sum(abs(B).^2,1));
  B = B ./ repmat(normB,[smodey(2) 1]);
  C = C .* repmat(normB,[smodey(3) 1]);
  % save B and C
  comp{2} = B;
  comp{3} = C;
end
%%%%%%%%%%%%%


%%%%%%%%%%%%%
if strcmp(normmode,'normC')
  % normalize C onto A
  C = comp{3};
  A = comp{1};
  % get norm, remove from C and apply to A
  normC = sqrt(sum(abs(C).^2,1));
  C = C ./ repmat(normC,[smodey(3) 1]);
  A = A .* repmat(normC,[smodey(1) 1]);
  % save C and A
  comp{3} = C;
  comp{1} = A;
end
%%%%%%%%%%%%%


%%%%%%%%%%%%%
if strcmp(normmode,'normD')
  switch Dmode
    case 'identity'
      error('normalizing D should not be performed when Dmode = ''identity''')
    case 'kdepcomplex'
      % normalize D onto B onto C
      D = comp{5};
      B = comp{2};
      C = comp{3};
      % get norm, remove from D and apply to B, and remove from D and apply to B
      for ik = 1:smodey(2)
        currD = squeeze(D(ik,:,:));
        currDnorm = sqrt(sum(abs(currD).^2,1));
        currDnorm(currDnorm<eps*1e6) = 1; % when one B is zero, D cannot be found and has no norm
        D(ik,:,:) = currD ./ repmat(currDnorm,[ncomp 1]);
        B(ik,:)   = B(ik,:) .* currDnorm;
      end
      normB = sqrt(sum(abs(B).^2,1));
      B = B ./ repmat(normB,[smodey(2) 1]);
      C = C .* repmat(normB,[smodey(3) 1]);
      % save D, B and C
      comp{5} = D;
      comp{2} = B;
      comp{3} = C;
  end
end
%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction optimizesigma       %%%%%%%%%%%%
function [sigma,newtcount,lscount] = optimizesigma(sigmaind,sigma,alpha,beta,theta,gamma,eta,F,sigmacirc,cpairind,precision,maxoptimiter,Dmode)
%
% Find optimum for sigma
%
% Code below tries to do on newton-raphson, and does a line search when it can't
%
%
% Some relevant benchmarks:
%    5 component problem, ratio: 1:3.1:7.3
%     Loss function evaluation: 0.000167 seconds
%         Gradient calculation: 0.000520 seconds
% Gradient+Hessian calculation: 0.001214 seconds
%   10 component problem, ratio: 1:3.3:8.4
%     Loss function evaluation: 0.000222 seconds
%         Gradient calculation: 0.000726 seconds
% Gradient+Hessian calculation: 0.001855 seconds
%
%
% set n's
ncomp = numel(sigma);
sigmastart = sigma;

% initialize
maxstep       = (1/2)/F(end);
optimiter     = 0;
% start counting optim bits
newtcount     = 0;
lscount       = 0;
nomorenewtflg = 0;
nomorelsflg   = 0;
% initialize gradient and hessian
[siggradient,sighessian] = sigmagradhess(sigma,beta,theta,gamma,eta,F,cpairind,sigmaind,Dmode);
newtstep  = (sighessian\siggradient).';
startloss = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
while any(abs(newtstep) > precision) && (optimiter<maxoptimiter)
  optimiter = optimiter + 1;
  
  % if hessian is positive definite, continue with newton raphson
  if all(eig(sighessian)>eps*1e6) && ~nomorenewtflg
    newtcount = newtcount + 1;
    % make sure step is not bigger than maxstep
    stepsizenewt  = 1;
    while any(abs(stepsizenewt*newtstep)>maxstep)
      stepsizenewt = stepsizenewt / 2;
    end
    % lower stepsizenewt when lossfunval increases
    currloss = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
    sigmanewt = sigma;
    sigmanewt(sigmaind) = sigma(sigmaind) - stepsizenewt.*newtstep;
    lossnewt = sigmalossfun(sigmanewt,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
    newtstepsziter = 0;
    if lossnewt == currloss
      % difference in sigma's is so low, that the differences they cause in parts of the lossfun dissappear
      % in the total lossfun. This can easily happen if sigma differences (1e-6) create tiny differences (1e-12)
      % in the input to cos if the input is close to pi/0: when alpha is ~1e6 the final differences are smaller
      % than eps, and the difference between lossfunctions becomes zero. A break is thus warrented.
      break
    end
    while (lossnewt > currloss) && (lossnewt>eps*1e6) % if optimization fails, an extra check here on stepsizenewt.*newtstep > precision will likely solve the problem!
      newtstepsziter = newtstepsziter + 1;
      % lower stepsizenewt
      stepsizenewt = stepsizenewt / 2;
      sigmanewt = sigma;
      sigmanewt(sigmaind) = sigma(sigmaind) - stepsizenewt.*newtstep;
      lossnewt = sigmalossfun(sigmanewt,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
      if newtstepsziter == 100;
        disp('warning: decreasing newtstepsize took too long')
      end
    end
    [siggradientnewt,sighessiannewt] = sigmagradhess(sigmanewt,beta,theta,gamma,eta,F,cpairind,sigmaind,Dmode);
    % implement circularity
    onlocircbound = sigmanewt(sigmaind)<(eps*1e6);
    onhicircbound = sigmanewt(sigmaind)>(sigmacirc-(eps*1e6));
    if any(onlocircbound) || any(onhicircbound)
      lossprecirc = lossnewt;
      % correct in complex domain
      radconv = ((2*pi)./sigmacirc);
      Scomp = exp(1i*sigmanewt.*radconv);
      sigmanewt      = angle(Scomp);
      sigmanewt(sigmanewt<0) = sigmanewt(sigmanewt<0) + (2*pi);
      sigmanewt      = sigmanewt ./ radconv;
      [siggradientnewt,sighessiannewt] = sigmagradhess(sigmanewt,beta,theta,gamma,eta,F,cpairind,sigmaind,Dmode);
      losspostcirc = sigmalossfun(sigmanewt,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
      if (abs(lossprecirc-losspostcirc)/lossprecirc)>precision
        warning(['sigma circularity implementation failed in newton-raphson: relative difference is ' num2str(abs(lossprecirc-losspostcirc)/lossprecirc)])
      end
    end
    % check whether hessian is no longer positive definite
    if all(eig(sighessiannewt)>eps*1e6)
      if any(abs(sigmanewt(sigmaind)-sigma(sigmaind))>precision)
        sigma = sigmanewt;
        siggradient = siggradientnewt;
        sighessian  = sighessiannewt;
        newtstep = (sighessian\siggradient).';
      else
        % difference in sigma has reached precision level, stopping and keeping previous estimate. Note, this assumes that because the hessian is positive
        % definite, that newton-raphson is always the best way to optimize. A break is then justified, as line search is not needed after this step.
        break
      end
    else
      nomorenewtflg = 1;
    end
  elseif ~nomorelsflg % gradient descent with line search, perform only once
    % first, find a point b away from a such that the sign of the derivative of the step size function at b is the opposite of that at a
    a = 0;
    b = (maxstep/2)./max(abs(siggradient)); % start with a point maximally 1/4th of a cycle of the higest F away from a
    % evaluate derivative of the step size function at a and b
    sigmasearch = sigma;
    sigmasearch(sigmaind) = sigma(sigmaind) - a.*siggradient.';
    aderiv = sigmagradhess(sigmasearch,beta,theta,gamma,eta,F,cpairind,sigmaind,Dmode).' * -siggradient;
    sigmasearch(sigmaind) = sigma(sigmaind) - b.*siggradient.';
    bderiv = sigmagradhess(sigmasearch,beta,theta,gamma,eta,F,cpairind,sigmaind,Dmode).' * -siggradient;
    bracketiter = 0;
    while sign(aderiv) == sign(bderiv)
      bracketiter = bracketiter + 1;
      b = b*1.2;
      sigmasearch = sigma;
      sigmasearch(sigmaind) = sigma(sigmaind) - b.*siggradient.';
      bderiv = sigmagradhess(sigmasearch,beta,theta,gamma,eta,F,cpairind,sigmaind,Dmode).' * -siggradient;
      if bracketiter == 100
        disp('warning: initial bracketing before line search took too long')
      end
      % sanity checks
      if any(abs(b*siggradient)>maxstep*5) % maxstep*5 should be two and half full maximum cycles (changed to *5 because sometimes there is a saddlepoint in between two local minima)
        disp('warning: initial bracketing before line search possibly jumped out of intended minimum')
      end
    end % while ~(lossb<lossa)
    % bracket found, start bisection
    %
    % Below there are two algorithms for finding the step size with the lowest loss function.
    % The first finds this point using the sign of the derivative (gradient) of the loss function given the step size,
    % the other finds this point by just evaluating the loss function at each step size.
    % The first can be very slow because it calculates the gradient at every iteration, which is computationally intense
    % when Dmode = kdepcomplex. However, using the first can result in finding the global minimum (of the main loop) in much
    % less iterations. Moreover, the second has a lot of difficulties relating to machine precision in the loss function,
    % and the step sizes themselves. It can probably be improved quite a bit. Things to keep in mind there are:
    % 1) the precision of the cosines/sines in the loss function: small lambda's can have no discernable effect on the loss
    % 2) tiny stepsizes can have a large effect on the loss if the gradients are huge. But the stepsizes can easily approach eps...
    % The first (using sign of deriv) doesn't have a stop condition relating to eps of stepsize, while the second does. I don't
    % remember why this was necessary in the second, but not in the first. This likely relates to the above two points.
    %
    % perform bisection on the derivative of the step size function
    %maxbsiter = ceil(log2((b-a)./precision)+1); % precision guaranteed to be reached in this many iterations
    bsiter = 0;
    while (any(abs((b.*siggradient)-(a.*siggradient)) > precision)) && (bsiter < (maxoptimiter))
      lscount = lscount + 1;
      bsiter  = bsiter + 1;
      % set middle point x
      x = (a+b)/2;
      % evaluate derivative of step size function at x
      sigmasearch = sigma;
      sigmasearch(sigmaind) = sigma(sigmaind) - x.*siggradient.';
      xderiv = sigmagradhess(sigmasearch,beta,theta,gamma,eta,F,cpairind,sigmaind,Dmode).' * -siggradient;
      if (sign(xderiv) == sign(aderiv)) && (sign(xderiv) == sign(bderiv))
        error('something went wrong during line search')
      end
      % select the two points with opposite sign
      if sign(xderiv) == sign(aderiv) % next set will be x,b
        a = x;
        aderiv = xderiv;
      elseif sign(xderiv) == sign(bderiv) % next set will be x,a
        b = x;
        bderiv = xderiv;
      end
    end % line search
    % select lowest deriv of the two
    if abs(aderiv)<abs(bderiv)
      x = a;
    elseif abs(bderiv)<abs(aderiv)
      x = b;
    end
    %     % perform bisection on the step size of the loss function
    %     sigmasearch = sigma;
    %     sigmasearch(sigmaind) = sigma(sigmaind) - a.*siggradient.';
    %     lossa = sigmalossfun(sigmasearch,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
    %     sigmasearch = sigma;
    %     sigmasearch(sigmaind) = sigma(sigmaind) - b.*siggradient.';
    %     lossb = sigmalossfun(sigmasearch,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
    %     % if a and b are extremely close together, start with a
    %     x = a;
    %     bsiter = 0;
    %     while (any(abs((b.*siggradient)-(a.*siggradient)) > precision)) && (bsiter < (maxoptimiter)) && (abs((lossa-lossb)./lossa) > eps*1e6) && (abs(a-b) > eps*1e6)
    %       lscount = lscount + 1;
    %       bsiter  = bsiter + 1;
    %       % set middle point x
    %       x = (a+b)/2;
    %       % evaluate loss function function at x
    %       sigmasearch = sigma;
    %       sigmasearch(sigmaind) = sigma(sigmaind) - x.*siggradient.';
    %       lossx = sigmalossfun(sigmasearch,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
    %       if (lossx > lossa) && (lossx > lossb)
    %         error('something went wrong during line search')
    %       end
    %       if lossb>lossa % next is a and x, x is should always be lowest
    %         b = x;
    %         lossb = lossx;
    %       else % next is b and x, x should always be lowest
    %         a = x;
    %         lossa = lossx;
    %       end
    %     end % line search
    % select lowest loss of the two
    %x = x; % x is always lowest
    % update sigma and break main while loop
    sigma(sigmaind) = sigma(sigmaind) - x.*siggradient.';
    [siggradient sighessian] = sigmagradhess(sigma,beta,theta,gamma,eta,F,cpairind,sigmaind,Dmode);
    newtstep = (sighessian\siggradient).';
    % implement circularity
    onlocircbound = sigma(sigmaind)<(eps*1e6);
    onhicircbound = sigma(sigmaind)>(sigmacirc-(eps*1e6));
    if any(onlocircbound) || any(onhicircbound)
      lossprecirc = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
      % correct in complex domain
      radconv = ((2*pi)./sigmacirc);
      Scomp = exp(1i*sigma.*radconv);
      sigma      = angle(Scomp);
      sigma(sigma<0) = sigma(sigma<0) + (2*pi);
      sigma      = sigma ./ radconv;
      [siggradient sighessian] = sigmagradhess(sigma,beta,theta,gamma,eta,F,cpairind,sigmaind,Dmode);
      newtstep = (sighessian\siggradient).';
      losspostcirc = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
      if (abs(lossprecirc-losspostcirc)/lossprecirc)>precision
        warning(['sigma circularity implementation failed in gradient descent: difference is ' num2str(abs(lossprecirc-losspostcirc)/lossprecirc)])
      end
    end
    % if hessian is not positive definite, quit main loop
    if ~all(eig(sighessian)>eps*1e6)
      break
    else
      % if hessian is positve, try to do another newton-raphson, but prevent doing another line search
      nomorelsflg = 1;
    end
  end % if all(eig(sighessian)>precision)
  if nomorelsflg && nomorenewtflg % when both have been done to the best of their abilities, quit
    break
  end
end % main while loop
if (optimiter==maxoptimiter)
  disp('warning: sigma optimization ended because maximum number iterations was reached')
end
% failsafe on lossfunval
endloss = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,cpairind,1:ncomp,Dmode);
if (((endloss-startloss)./endloss)>eps*1e6) && (endloss>eps*1e6) && (startloss>eps*1e6)
  warning('optimization failed, using initial estimate')
  sigma = sigmastart;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction calcsigmagrad       %%%%%%%%%%%%
function [siggradient,sighessian] = sigmagradhess(sigma,beta,theta,gamma,eta,F,spairind,sigmaind,Dmode)
%
% Calculate gradient and hessian matrices of given sigma
% Sigma should always be 1xNComp.
%
% sigmaconst is an 1xN vector indicating which sigma's are considered as active
%
%

% select which sigma's are active
nsigma = numel(sigma);
actsigma = sigmaind;
nactsigma = numel(actsigma);
nactspair = (nactsigma.*(nactsigma-1))./2;

% switch between Dmode
switch Dmode
  
  case 'identity'
    % switch between ncomp == 1 and ncomp > 1
    if nsigma ~= 1
      
      %%%%%%%%%%%%%%%
      %%% Gradient (and some ingredients for Hessian)
      % get pairwise sigma differences and expand gamma and eta for:
      % first order partial derive
      % pairs of identical sigma's of second order partial derivative
      % part2 depends on component interaction and is always zero (because gamma is always zero)
      
      % calculate first order partial derivatives
      % part1 of first order partial derivative
      part1 = sum(beta(:,actsigma) .* bsxfun(@times,(2*pi*F),-sin(2*pi*F*sigma(actsigma) + theta(:,actsigma)))).';
      % part2 of first order partial derivative
      % part2 depends on component interaction and is always zero (because gamma is always zero)
      % store first order partial derivative
      siggradient = part1;
      %%%%%%%%%%%%%%%
      
      
      %%%%%%%%%%%%%%%
      %%% Hessian
      if nargout == 2
        % get pairwise sigma differences and expand gamma and eta for:
        % pairs of non-identical sigma's of second order partial derivative
        % part2 depends on component interaction and is always zero (because gamma is always zero)
        % allocate
        sighessian  = zeros(nactsigma);
        
        % calculate lower diagonal of hessian matrix
        % part1 of second order partial derivative for pairs of non-identical sigma's
        % part1 falls away
        % part2 of second order partial derivative for pairs of non-identical sigma's
        % part2 depends on component interaction and is always zero (because gamma is always zero)
        
        % calculate diagonal of hessian matrix
        % part1 of second order partial derivative for pairs of identical sigma's
        part1 = sum(beta(:,actsigma) .* bsxfun(@times,(2*pi*F).^2,-cos(2*pi*F*sigma(actsigma) + theta(:,actsigma)))).';
        % part2 of second order partial derivative for pairs of identical sigma's
        % part2 depends on component interaction and is always zero (because gamma is always zero)
        sighessian(diag(true(1,nactsigma))) = part1;
      end
      %%%%%%%%%%%%%%%
      
    elseif nsigma == 1
      %%%%%%%%%%%%%%%
      %%% First order derivative (still called gradient here, so optimizesigma can be used when ncomp == 1)
      % calculate first order derivative
      siggradient = sum(beta .* (2*pi*F)  .* -sin((2*pi*F) .* sigma + theta));
      %%%%%%%%%%%%%%%
      
      
      %%%%%%%%%%%%%%%
      %%% Second order derivative (still called hessian here, so optimizesigma can be used when ncomp == 1)
      if nargout == 2
        % calculate second order derivative
        sighessian = sum(beta .* (2*pi*F).^2 .* -cos((2*pi*F) .* sigma + theta)).';
      end
      %%%%%%%%%%%%%%%
    end
    
    
  case 'kdepcomplex'
    % switch between ncomp == 1 and ncomp > 1
    if nsigma ~= 1
      
      %%%%%%%%%%%%%%%
      %%% Gradient (and some ingredients for Hessian)
      % get pairwise sigma differences and expand gamma and eta for:
      % first order partial derive
      % pairs of identical sigma's of second order partial derivative
      spairindexp = NaN(nactsigma,nsigma-1);
      pairsign    = NaN(nactsigma,nsigma-1);
      for isigma = 1:nactsigma
        iactsigma = actsigma(isigma);
        currpairind = find(sum(spairind==iactsigma,2));
        spairindexp(isigma,:) = currpairind;
        pairsign(isigma,:)    = -1 + (spairind(currpairind,2)~=iactsigma).'*2;
      end
      sigmapairdiff1 = reshape(sigma(spairind(spairindexp,1))-sigma(spairind(spairindexp,2)),[nactsigma nsigma-1]).';
      gammaexp1sign  = bsxfun(@times,reshape(gamma(:,spairindexp.'),[numel(F) nsigma-1 nactsigma]),reshape(pairsign.',[1 nsigma-1 nactsigma]));
      gammaexp1      = reshape(gamma(:,spairindexp.'),[numel(F) nsigma-1 nactsigma]);
      etaexp1        = reshape(eta(:,spairindexp.'),[numel(F) nsigma-1 nactsigma]);
      
      % calculate first order partial derivatives
      % part1 of first order partial derivative
      part1 = sum(beta(:,actsigma) .* bsxfun(@times,(2*pi*F),-sin(2*pi*F*sigma(actsigma) + theta(:,actsigma)))).';
      % part2 of first order partial derivative
      part2 = sum(sum(bsxfun(@times,gammaexp1sign,bsxfun(@times,(2*pi*F),-sin(bsxfun(@plus,bsxfun(@times,2*pi*F,reshape(sigmapairdiff1,[1 nsigma-1 nactsigma])),etaexp1)))),1),2);
      % store first order partial derivative
      siggradient = part1 + part2(:);
      %%%%%%%%%%%%%%%
      
      
      %%%%%%%%%%%%%%%
      %%% Hessian
      if nargout == 2
        % get pairwise sigma differences and expand gamma and eta for:
        % pairs of non-identical sigma's of second order partial derivative
        if nactsigma>1
          spairindexp = NaN(nactspair,1);
          cnt = 0;
          for isigma1 = 1:nactsigma-1
            for isigma2 = isigma1+1:nactsigma
              cnt = cnt+1;
              iactsigma1 = actsigma(isigma1);
              iactsigma2 = actsigma(isigma2);
              spairindexp(cnt) =  find((sum(spairind==iactsigma1,2) + sum(spairind==iactsigma2,2))==2);
            end
          end
          sigmapairdiff2 = reshape(sigma(spairind(spairindexp,1))-sigma(spairind(spairindexp,2)),[nactspair 1]).';
          gammaexp2      = reshape(gamma(:,spairindexp.'),[numel(F) 1 nactspair]);
          etaexp2        = reshape(eta(:,spairindexp.'),[numel(F) 1 nactspair]);
        end
        % allocate
        sighessian  = NaN(nactsigma);
        
        % calculate lower diagonal of hessian matrix
        % part1 of second order partial derivative for pairs of non-identical sigma's
        % part1 falls away
        % part2 of second order partial derivative for pairs of non-identical sigma's
        if nactsigma>1
          part2 = sum(bsxfun(@times,gammaexp2,bsxfun(@times,-(2*pi*F).^2,-cos(bsxfun(@plus,bsxfun(@times,2*pi*F,reshape(sigmapairdiff2,[1 1 nactspair])),etaexp2)))),1);
          sighessian(tril(true(nactsigma),-1)) = part2;
          tmp = sighessian.';
          sighessian(triu(true(nactsigma))) = tmp(triu(true(nactsigma)));
        end
        
        % calculate diagonal of hessian matrix
        % part1 of second order partial derivative for pairs of identical sigma's
        part1 = sum(beta(:,actsigma) .* bsxfun(@times,(2*pi*F).^2,-cos(2*pi*F*sigma(actsigma) + theta(:,actsigma)))).';
        % part2 of second order partial derivative for pairs of identical sigma's
        part2 = sum(sum(bsxfun(@times,gammaexp1,bsxfun(@times,(2*pi*F).^2,-cos(bsxfun(@plus,bsxfun(@times,2*pi*F,reshape(sigmapairdiff1,[1 nsigma-1 nactsigma])),etaexp1)))),1),2);
        sighessian(diag(true(1,nactsigma))) = part1 + part2(:);
      end
      %%%%%%%%%%%%%%%
      
    elseif nsigma == 1
      %%%%%%%%%%%%%%%
      %%% First order derivative (still called gradient here, so optimizesigma can be used when ncomp == 1)
      % calculate first order derivative
      siggradient = sum(beta .* (2*pi*F)  .* -sin((2*pi*F) .* sigma + theta));
      %%%%%%%%%%%%%%%
      
      
      %%%%%%%%%%%%%%%
      %%% Second order derivative (still called hessian here, so optimizesigma can be used when ncomp == 1)
      if nargout == 2
        % calculate second order derivative
        sighessian = sum(beta .* (2*pi*F).^2 .* -cos((2*pi*F) .* sigma + theta)).';
      end
      %%%%%%%%%%%%%%%
    end
    
end % switch Dmode

%%%%%%%%%
% TEST CODE (about twice as slow for a 20x25x30 problem with 5 components)
% % allocate
% oldsiggradient = NaN(nactsigma,1);
% oldsighessian  = NaN(nactsigma);
%
% % loop over sigma to calculate first order partial derivatives
% for isigma = 1:nactsigma
%   iactsigma    = actsigma(isigma);
%   currsigma    = sigma(iactsigma);
%   ispair       = find(sum(spairind==iactsigma,2));
%   currspairind = spairind(ispair,:);
%   pairsign     = -1 + (currspairind(:,2)~=iactsigma).'*2;
%
%   % part1 of first order partial derivative
%   part1 = sum(beta(:,iactsigma) .* (2*pi*F) .* -sin(2*pi*F.*currsigma + theta(:,iactsigma)));
%
%   % part2 of first order partial derivative
%   part2 = sum(sum(gamma(:,ispair) .* bsxfun(@times,(2*pi*F),bsxfun(@times,pairsign,-sin(bsxfun(@times,2*pi*F,(sigma(currspairind(:,1))-sigma(currspairind(:,2)))) + eta(:,ispair)))))); % -sign only present when it's the second of 2 pairs
%
%   % store first order partial derivative
%   oldsiggradient(isigma) = part1 + part2;
% end
%
% % loop over sigma pairs to calculate lower diagonal of hessian matrix
% if nactsigma > 1
%   for isigma1 = 1:nactsigma-1
%     for isigma2 = isigma1+1:nactsigma
%       iactsigma1 = actsigma(isigma1);
%       iactsigma2 = actsigma(isigma2);
%       ispair       = find((sum(spairind==iactsigma1,2) + sum(spairind==iactsigma2,2))==2);
%       currspairind = spairind(ispair,:);
%
%       % part1 of second order partial derivative
%       % part1 falls away
%       % part2 of second order partial derivative
%       part2 = sum(sum(gamma(:,ispair) .* -(2*pi*F).^2 .* -cos(bsxfun(@times,2*pi*F,(sigma(currspairind(:,1))-sigma(currspairind(:,2)))) + eta(:,ispair)))); % -sign always present
%
%       % store second order partial derivative in hessian matrix
%       oldsighessian(isigma1,isigma2) = part2;
%     end
%   end
% end
% % loop over pairs of identical sigma's  to calculate diagonal of hessian matrix
% for isigma = 1:nactsigma
%   iactsigma    = actsigma(isigma);
%   currsigma    = sigma(iactsigma);
%   ispair       = find(sum(spairind==iactsigma,2));
%   currspairind = spairind(ispair,:);
%
%   % part1 of second order partial derivative
%   part1 = sum(beta(:,iactsigma) .* (2*pi*F).^2 .* -cos(2*pi*F.*currsigma + theta(:,iactsigma)));
%
%   % part2 of second order partial derivative
%   part2 = sum(sum(gamma(:,ispair) .* bsxfun(@times,(2*pi*F).^2, -cos(bsxfun(@times,2*pi*F,(sigma(currspairind(:,1))-sigma(currspairind(:,2)))) + eta(:,ispair))))); % -sign never present
%
%   % store second order partial derivative in hessian matrix
%   oldsighessian(isigma,isigma) = part1 + part2;
% end
% % make sighessian diagonally symmetric
% tmp = oldsighessian';
% oldsighessian(tril(true(nactsigma))) = tmp(tril(true(nactsigma)));
% sighessian-oldsighessian
% siggradient-oldsiggradient


% % use finite differences to check above
% if exist('alpha','var')
%   delta = [1e-2 1e-4 1e-6];
%   % gradient first
%   fx = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,spairind,1:nsigma,Dmode);
%   gradientfd = NaN(nactsigma,3);
%   for idelt = 1:3
%     for isigma = 1:nactsigma
%       i    = actsigma(isigma);
%       sigip = sigma;
%       sigip(i) = sigip(i) + delta(idelt);
%       sigim = sigma;
%       sigim(i) = sigim(i) - delta(idelt);
%       fxsigip = sigmalossfun(sigip,alpha,beta,theta,gamma,eta,F,spairind,1:nsigma,Dmode);
%       fxsigim = sigmalossfun(sigim,alpha,beta,theta,gamma,eta,F,spairind,1:nsigma,Dmode);
%       gradientfd(isigma,idelt) = (fxsigip-fxsigim)/(2*delta(idelt));
%     end
%   end
%   cat(2,siggradient,gradientfd)
%   hessianfd = NaN(3,nactsigma*nactsigma);
%   for idelt = 1:3
%     hesstmp = NaN(nactsigma);
%     % loop over sigma pairs to calculate lower diagonal of hessian matrix
%     for isigma1 = 1:nactsigma-1
%       for isigma2 = isigma1+1:nactsigma
%         i = actsigma(isigma1);
%         j = actsigma(isigma2); % partial derivative of i of the partial derivative of j
%         sigipjm = sigma;
%         sigipjm(i) = sigipjm(i) + delta(idelt);
%         sigipjm(j) = sigipjm(j) - delta(idelt);
%         sigimjp = sigma;
%         sigimjp(i) = sigimjp(i) - delta(idelt);
%         sigimjp(j) = sigimjp(j) + delta(idelt);
%         sigipjp = sigma;
%         sigipjp([i j]) = sigipjp([i j]) + delta(idelt);
%         sigimjm = sigma;
%         sigimjm([i j]) = sigimjm([i j]) - delta(idelt);
%         fxsigipjm = sigmalossfun(sigipjm,alpha,beta,theta,gamma,eta,F,spairind,1:nsigma,Dmode);
%         fxsigimjp = sigmalossfun(sigimjp,alpha,beta,theta,gamma,eta,F,spairind,1:nsigma,Dmode);
%         fxsigipjp = sigmalossfun(sigipjp,alpha,beta,theta,gamma,eta,F,spairind,1:nsigma,Dmode);
%         fxsigimjm = sigmalossfun(sigimjm,alpha,beta,theta,gamma,eta,F,spairind,1:nsigma,Dmode);
%         hesstmp(isigma1,isigma2) = (fxsigipjp - fxsigipjm - fxsigimjp + fxsigimjm)/(4*delta(idelt).*delta(idelt));
%       end
%     end
%     % loop over pairs of identical sigma's  to calculate diagonal of hessian matrix
%     for isigma = 1:nactsigma
%       i = actsigma(isigma); % partial derivative of i of the partial derivative of i
%       sigip = sigma;
%       sigip(i) = sigip(i) + delta(idelt);
%       sigim = sigma;
%       sigim(i) = sigim(i) - delta(idelt);
%       fxsigip = sigmalossfun(sigip,alpha,beta,theta,gamma,eta,F,spairind,1:nsigma,Dmode);
%       fxsigim = sigmalossfun(sigim,alpha,beta,theta,gamma,eta,F,spairind,1:nsigma,Dmode);
%       hesstmp(isigma,isigma) = (fxsigip-(2*fx)+fxsigim)/(delta(idelt).*delta(idelt));
%     end
%     tmp = hesstmp';
%     hesstmp(tril(true(nactsigma))) = tmp(tril(true(nactsigma)));
%     hessianfd(idelt,:) = hesstmp(:);
%   end
%   cat(2,sighessian(:),hessianfd.')
% end
%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction sigmalossfun        %%%%%%%%%%%%
function [lossfunval] = sigmalossfun(sigma,alpha,beta,theta,gamma,eta,F,spairind,sigmaind,Dmode)
%
% Calculate loss function value given sigma
% Sigma should always be 1xNComp.
%
% sigmaind indicates which sigma's to exclusively use for calcuating lossfun value
% Except alpha, which will always be used. It's irrelevant, but will result in a
% lossfunval of 0<>eps when a perfect noiseless situation is used as a test case.
%

% select which sigma's are active
nsigma = numel(sigma);

% switch between Dmode
switch Dmode
  
  case 'identity'
    % switch between ncomp == 1 and ncomp > 1
    if nsigma ~= 1
      % part2 depends on component interaction and is always zero (because gamma is always zero)
      
      % calculate first part of loss function
      part1 = alpha + sum(sum(beta(:,sigmaind) .* cos(2*pi*F*sigma(sigmaind) + theta(:,sigmaind))));
      % calculate second part of loss function
      % part2 depends on component interaction and is always zero (because gamma is always zero)
      
      % output loss function value
      lossfunval = part1;
      
    elseif nsigma == 1
      % calculate loss function value
      lossfunval = alpha + sum(beta .* cos((2*pi*F) .* sigma + theta));
    end
    
  case 'kdepcomplex'
    % switch between ncomp == 1 and ncomp > 1
    if nsigma ~= 1
      if ~(nsigma==numel(sigmaind))
        % rmeove ignored parts from gamma, eta, spairind
        sigmarem = 1:nsigma;
        sigmarem(sigmaind) = [];
        nsigmarem = numel(sigmarem);
        remind = cell(1,nsigmarem);
        for isigmarem = 1:nsigmarem
          currsigmarem = sigmarem(isigmarem);
          remind{isigmarem} = find(sum(spairind==currsigmarem,2));
        end
        remind = vertcat(remind{:});
        gamma(:,remind) = [];
        eta(:,remind) = [];
        spairind(remind,:) = [];
      end
      
      % calculate first part of loss function
      part1 = alpha + sum(sum(beta(:,sigmaind) .* cos(2*pi*F*sigma(sigmaind) + theta(:,sigmaind))));
      % calculate second part of loss function
      part2 = sum(sum(gamma .* cos(2*pi*F*(sigma(spairind(:,1))-sigma(spairind(:,2))) + eta)));
      
      % output loss function value
      lossfunval = part1 + part2;
      
    elseif nsigma == 1
      % calculate loss function value
      lossfunval = alpha + sum(beta .* cos((2*pi*F) .* sigma + theta));
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction calcmagmodely          %%%%%%%%%%%%
function [model] = calcmagmodely(comp,smodey,compind,Dmode,ij,Bind,Cind,Dind)

% sanity check
if numel(compind)~=1
  error('calcmodely only suited for single component models')
end

% calcute based on Dmode
switch Dmode
  
  case 'identity'
    % Calculate model using kr-products
    model = zeros(smodey(2).*smodey(3).*smodey(4),1);
    
    % compute requested model
    for icomp = 1:numel(compind)
      currcomp = compind(icomp);
      
      % expand and element-wise multiply
      compmodel = kr(comp{5}(:,currcomp),comp{3}(:,currcomp),comp{2}(:,currcomp),comp{1}(ij,currcomp));
      
      % save model
      model = model + compmodel;
    end
    % reshape
    model = reshape(model, smodey([2 3 4]));
    
  case 'kdepcomplex'
    % Calculate model using vectorized expansion
    % set Aind
    Aind = repmat(ij,[1 prod(smodey([2 3 4]))]).';
    model = zeros(smodey(2).*smodey(3).*smodey(4),1);
    
    % compute requested model
    for icomp = 1:numel(compind)
      currcomp = compind(icomp);
      
      % expand and element-wise multiply
      % special treatment for D's
      currD = comp{5}(:,:,currcomp);
      currD = currD(:);
      compmodel = comp{1}(Aind,currcomp) .* comp{2}(Bind,currcomp) .* comp{3}(Cind,currcomp) .* abs(currD(Dind));
      
      % save model
      model = model + compmodel;
    end
    % reshape
    model = reshape(model, smodey([2 3 4]));
    
end


% OLD
% % Calculate model using vectorized massive singleton expansion
% model = zeros(smodey(1).*smodey(2).*smodey(3).*smodey(4),1);
%
% if numel(compind)~=1
%   error('calcmodely only suited for single component models')
% end
%
% % compute requested model
% for icomp = 1:numel(compind)
%   currcomp = compind(icomp);
%
%   % expand and element-wise multiply
%   switch Dmode
%     case 'identity'
%       compmodel = kr(comp{5}(:,currcomp),comp{3}(:,currcomp),comp{2}(:,currcomp),comp{1}(:,currcomp));
%     case 'kdepcomplex'
%       % special treatment for D's
%       currD = comp{5}(:,:,currcomp);
%       currD = currD(:);
%       compmodel = comp{1}(Aind,currcomp) .* comp{2}(Bind,currcomp) .* comp{3}(Cind,currcomp) .* abs(currD(Dind));
%   end
%
%   % save model
%   model = model + compmodel;
% end
% % reshape
% model = reshape(model, smodey);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction calcPkl             %%%%%%%%%%%%
function [Pkl] = calcPkl(comp,F,smode,datforQ,Dmode)
%
% Calculate Pkl by Q
%
switch Dmode
  case 'identity'
    % nothing necessary
  case 'kdepcomplex'
    % permute Dk
    comp{5} = permute(comp{5},[2 3 1]);
end
Pkl = cell(smode([2 3]));
for ik = 1:smode(2)
  ASk = comp{1} .* exp(1i*2*pi*F(ik)*comp{4});
  for il = 1:smode(3)
    % select currdatforQ
    currdatforQ = datforQ{ik,il};
    currdatforQ = double(currdatforQ); % ensure double precision
    
    % calculate partial model Z for Q
    switch Dmode
      case 'identity'
        ZforQ = ASk * diag(comp{2}(ik,:)) * diag(comp{3}(il,:)) * comp{5};
      case 'kdepcomplex'
        ZforQ = ASk * diag(comp{2}(ik,:)) * diag(comp{3}(il,:)) * comp{5}(:,:,ik).';
    end
    
    % calculate Qkl
    Qkl = currdatforQ' * ZforQ;
    % calculate Pkl
    [U,S,V] = svd(Qkl,'econ');
    currP = U*V';
    
    % save Pkl
    Pkl{ik,il} = currP;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction calcYdat            %%%%%%%%%%%%
function Ydat = calcYdat(datforQ,Pkl,smodey)
%
% Calculate Ydat
%
% calculate Y
Ydat = zeros(smodey);
for ik = 1:smodey(2)
  for il = 1:smodey(3)
    % select currdatforQ
    currdatforQ = datforQ{ik,il};
    currdatforQ = double(currdatforQ); % ensure double precision
    
    % calculate Ykl
    Ydat(:,ik,il,:) = currdatforQ * Pkl{ik,il}; % Pkl is always double precision
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction calcssqressparse          %%%%%%%%%%%%
function ssqres = calcssqressparse(comp,Pkl,F,smode,datforQ,Dmode)
%
% Calculate sums of squared resisuals by looping over ik and il, i.e. 'memory sparsely'.
% When using data-arrays of e.g. >1GB this not only saves memory, but
% also reduces computation time.
%

switch Dmode
  case 'identity'
    ssqres = 0;
    for ik = 1:smode(2)
      AS = comp{1} .* exp(1i*2*pi*F(ik)*comp{4});
      for il = 1:smode(3)
        
        % calculate model
        currPklDk = comp{5}.' * Pkl{ik,il}';
        modelx = AS * diag(comp{3}(il,:)) * diag(comp{2}(ik,:)) * currPklDk;
        
        % select data
        currdatforQ = datforQ{ik,il};
        currdatforQ = double(currdatforQ); % ensure double precision
        ssqres = ssqres + sum(sum(abs(currdatforQ-modelx).^2,1),2);
      end
    end
  case 'kdepcomplex'
    ssqres = 0;
    % permute Dk
    Dk = permute(comp{5},[2 3 1]);
    for ik = 1:smode(2)
      AS = comp{1} .* exp(1i*2*pi*F(ik)*comp{4});
      currDk = Dk(:,:,ik).';
      for il = 1:smode(3)
        
        % calculate model
        currPklDk = currDk * Pkl{ik,il}';
        modelx = AS * diag(comp{3}(il,:)) * diag(comp{2}(ik,:)) * currPklDk;
        
        % select data
        currdatforQ = datforQ{ik,il};
        currdatforQ = double(currdatforQ); % ensure double precision
        ssqres = ssqres + sum(sum(abs(currdatforQ-modelx).^2,1),2);
      end
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction Linear search ALS   %%%%%%%%%%%%%%%%%
function [newcomp,newssqres] = linsearch(datforQ,comp,prevcomp,F,smode,iter,prevssqres,ssqres,dispprefix,Dmode,holdparam,sigmacirc)
%
% This subfunction searches linearly for expected loading vectors
% to 'skip' several iterations
%
% ssqres is calculated precisely, without orthogonalizations
%
%
% ATM it's quite conservative
%
% This function is heavily inspired by the N-way toolbox (http://www.models.kvl.dk/source/nwaytoolbox)


% set which params to keep non-negative
nonnegparams = 1; % only A
nparam = 5;

% set which parameter sets to accelerate
paramaccind = ones(1,5);
paramaccind(holdparam) = 0;

% set delta and other things
delta = iter ^ 1/3;
switch Dmode
  case 'identity'
    paramaccind(5) = 0;
  case 'kdepcomplex'
    % nothing necessary
end
ncomp = size(comp{1},2);
ssqaccel = [];
paramaccind = find(paramaccind);

% set difference of current loadings with previous loadings
dcomp = cell(1,nparam);
for iparam = 1:nparam
  if iparam~=4
    dcomp{iparam} = comp{iparam} - prevcomp{iparam};
  else
    % take care of circularity
    curr = comp{iparam};
    prev = prevcomp{iparam};
    % get to complex domain
    radconv = ((2*pi)./sigmacirc);
    curr = exp(1i*curr.*radconv);
    prev = exp(1i*prev.*radconv);
    % get diff
    diffcomp = curr.*conj(prev);
    diffreal = angle(diffcomp);
    diffreal(diffreal<0) = diffreal(diffreal<0) + (2*pi);
    diffreal = diffreal ./ radconv;
    % save
    dcomp{iparam} = diffreal;
  end
end

% calculate expected ssqres of delta iterations later
comppred = comp;
for iparam = paramaccind
  if iparam~=4
    comppred{iparam} = comp{iparam} + delta*dcomp{iparam};
  else
    % apply delta in complex domain to take care of circularity
    % extract and convert
    radconv = ((2*pi)./sigmacirc);
    Scomp = exp(1i*comp{iparam}.*radconv);
    deltadcomp = exp(1i*dcomp{iparam}.*delta.*radconv);
    % apply delta
    Scomp = Scomp .* deltadcomp;
    % move to real
    S      = angle(Scomp);
    S(S<0) = S(S<0) + (2*pi);
    S      = S ./ radconv;
    % save
    comppred{iparam} = S;
  end
end
Pkl         = calcPkl(comppred,F,smode,datforQ,Dmode);
ssqresdelta = calcssqressparse(comppred,Pkl,F,smode,datforQ,Dmode); % subfunction for calculating ssqres

% set isnegflg
isnegflg = 0;
for inonnegparam = 1:numel(nonnegparams)
  isnegflg = isnegflg || any(sign(comppred{nonnegparams(inonnegparam)}(:))==-1);
end

% set ssqres and ssq collection variable for picking best delta at the end
ssqaccel{1} = [0 ssqres 0];
ssqaccel{2} = [delta ssqresdelta isnegflg];

% decrease delta when expected ssqres is larger than current ssqres
while (ssqresdelta > ssqres) && (delta >= 1)
  delta = delta * 0.25;
  if (delta >= 1)
    % calculate expected ssqres of delta iterations later
    comppred = comp;
    for iparam = paramaccind
      if iparam~=4
        comppred{iparam} = comp{iparam} + delta*dcomp{iparam};
      else
        % apply delta in complex domain to take care of circularity
        % extract and convert
        radconv = ((2*pi)./sigmacirc);
        Scomp = exp(1i*comp{iparam}.*radconv);
        deltadcomp = exp(1i*dcomp{iparam}.*delta.*radconv);
        % apply delta
        Scomp = Scomp .* deltadcomp;
        % move to real
        S      = angle(Scomp);
        S(S<0) = S(S<0) + (2*pi);
        S      = S ./ radconv;
        % save
        comppred{iparam} = S;
      end
    end
    Pkl   = calcPkl(comppred,F,smode,datforQ,Dmode);
    ssqresdelta = calcssqressparse(comppred,Pkl,F,smode,datforQ,Dmode); % subfunction for calculating ssqres
    isnegflg = 0;
    for inonnegparam = 1:numel(nonnegparams)
      isnegflg = isnegflg || any(sign(comppred{nonnegparams(inonnegparam)}(:))==-1);
    end
    ssqaccel{end+1} = [delta ssqresdelta isnegflg];
    % debugging  disp(['delta ' num2str(delta)])
    % debugging  disp(['ssqres   ' num2str(ssqresdelta)])
  else
    delta = delta * 2; % undo previous delta decrease for next while loop
    break
  end
end

% increase delta until ssqres no longer get's better, but only if the extrapolated ssqres is at least 2 times the previous ssqres increase
ssqresdeltaold = ssqres;
ssqresdeltanew = ssqresdelta;
prevssqresinc = prevssqres - ssqres;
while (ssqresdeltanew < ssqresdeltaold) && (((ssqresdeltaold - ssqresdeltanew) / prevssqresinc) > 2)
  delta = delta *1.25;
  % calculate expected ssqres of delta iterations later
  comppred = comp;
  for iparam = paramaccind
    if iparam~=4
      comppred{iparam} = comp{iparam} + delta*dcomp{iparam};
    else
      % apply delta in complex domain to take care of circularity
      % extract and convert
      radconv = ((2*pi)./sigmacirc);
      Scomp = exp(1i*comp{iparam}.*radconv);
      deltadcomp = exp(1i*dcomp{iparam}.*delta.*radconv);
      % apply delta
      Scomp = Scomp .* deltadcomp;
      % move to real
      S      = angle(Scomp);
      S(S<0) = S(S<0) + (2*pi);
      S      = S ./ radconv;
      % save
      comppred{iparam} = S;
    end
  end
  Pkl   = calcPkl(comppred,F,smode,datforQ,Dmode);
  ssqresdeltaold = ssqresdeltanew;
  ssqresdeltanew = calcssqressparse(comppred,Pkl,F,smode,datforQ,Dmode); % subfunction for calculating ssqres
  % debugging  disp(['inc loop done  ' num2str(delta) '  ' num2str(ssqresdeltanew)])
  if (ssqresdeltanew < ssqresdeltaold)
    isnegflg = 0;
    for inonnegparam = 1:numel(nonnegparams)
      isnegflg = isnegflg || any(sign(comppred{nonnegparams(inonnegparam)}(:))==-1);
    end
    ssqaccel{end+1} = [delta ssqresdeltanew isnegflg];
    % debugging  disp(['deltainc ' num2str(delta)])
    % debugging  disp(['ssqresinc   ' num2str(ssqresdeltanew)])
  end
end

% set number of ssqcalcs
ssqcalcnum = numel(ssqaccel);

% remove accels that led to neg values in nonnegparams
remind = [];
actconstr = 0;
for iaccel = 1:ssqcalcnum
  if ssqaccel{iaccel}(3)==1
    remind = [remind iaccel];
    actconstr = 1;
  end
end
ssqaccel(remind) = [];

% pick best delta by picking lowest ssqres
ssqaccel = sortrows(vertcat(ssqaccel{:}),2);
delta  = ssqaccel(1,1);
delta  = double(delta); % ensure double precision
deltassqres = ssqaccel(1,2);


% create new component loading matrices
newcomp = comp;
for iparam = paramaccind
  if iparam~=4
    newcomp{iparam} = comp{iparam} + delta*dcomp{iparam};
  else
    % apply delta in complex domain to take care of circularity
    % extract and convert
    radconv = ((2*pi)./sigmacirc);
    Scomp = exp(1i*comp{iparam}.*radconv);
    deltadcomp = exp(1i*dcomp{iparam}.*delta.*radconv);
    % apply delta
    Scomp = Scomp .* deltadcomp;
    % move to real
    S      = angle(Scomp);
    S(S<0) = S(S<0) + (2*pi);
    S      = S ./ radconv;
    % save
    newcomp{iparam} = S;
  end
end
newssqres = deltassqres;
if ~actconstr
  disp([dispprefix 'acceleration performed with delta = ' num2str(delta,'%-8.2f') ' and ssqres = '  num2str(deltassqres) ' ('  num2str(ssqcalcnum-1) ' ssq calcs)'])
else
  disp([dispprefix 'acceleration performed with delta = ' num2str(delta,'%-8.2f') ' and ssqres = '  num2str(deltassqres) ' ('  num2str(ssqcalcnum-1) ' ssq calcs) - active constraint'])
end



















