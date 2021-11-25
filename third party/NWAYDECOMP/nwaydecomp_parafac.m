function [comp,ssqres,expvar,scaling,tuckcongr,t3core] = nwaydecomp_parafac(dat,ncomp,varargin)

% NWAYDECOMP_PARAFAC is a low-level function for nd_nwaydecomposition and is used
% to perform nway-decomposition using the PARAFAC model. The model is modified such
% that (specified) parameter vectors can be real-valued with complex-valued input.
% It is described and used in the publication below (please cite when used):
%
%    van der Meij, R., Kahana, M. J., and Maris, E. (2012). Phase-amplitude coupling in
%        human ECoG is spatially distributed and phase diverse. Journal of Neuroscience,
%        32(1), 111-123.
%
%
% Use as
%   [comp,ssqres,expvar,scaling,tuckcongr,t3core] = nwaydecomp_parafac(dat,ncomp,...)
%
% Input:
%   dat   = array containing data to be decomposed
%   ncomp = number indicating number of components
%
% Output:
%   comp       = cell array containing component loadings per mode (all loading vectors have frobenius norm = 1)
%   ssqres     = sums of squares of the residuals
%   expvar     = percentage explained variance of the data by the model
%   scaling    = scaling coefficients belonging to the first mode (magnitude scaling) and the first complex mode (phase shift)
%                all component loading vectors have norm = 1 and, if complex, have magnitude weighted mean phase of 0 in the output
%                the scaling coefficients that have been removed are put in scaling, which is a 1x1 cell-array containing an ncomp*1 vector
%                if no complex modes are present and a 1x2 cell-array containing an extra ncomp*1 vector if they are
%   tuckcongr  = tuckers congruence coefficents between components, high values mean some components are highly correlated, which is a sign of
%                a degenerate model
%       t3core = Tucker3 core array, if the n-way parafac model is valid, this array should resemble an identity-array
%
%
%
%
% Additional options should be specified in key-value pairs and can be
%   'compmodes'   = vector with length equal to ndims(dat) with 0 or 1
%                   indicating whether component parameters should be complex or not
%   'niter'       = maximum number of iterations (default = 2500)
%   'convcrit'    = convergence criterion (default = 1e-6)
%   'startval'    = previously computed start-values
%   'dispprefix'  = prefix added to all disp-calls, handy when function is used in many loops after each other
%   'holdmodes'   = vector of 0s and 1s indicating whether certain modes are not updated in each ALS-iteration
%   'randomseed'  = scalar, seed to use for Matlab's random number generator
%
%
%
%
%
% This function is inspired by the N-way toolbox (http://www.models.kvl.dk/source/nwaytoolbox)
% and the Triple SPICE project (http://www.ece.umn.edu/users/nikos/public_html/3SPICE/code.html)
% A fantastic resource on PARAFAC/2 is written by Rasmus Bro http://www.models.kvl.dk/users/rasmus/
% (Multi-way Analysis book or monograph)
%
%
%
%   TO DO: merge some of the subfunctions of the PARAFAC models into externals
%
%
%
%

%
% Copyright (C) 2010-present, Roemer van der Meij, roemervandermeij AT gmail DOT com
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

% start execution timer
stopwatch = tic;

% Get the optional input arguments
keyvalcheck(varargin, 'optional', {'compmodes','niter','convcrit','startval','dispprefix','holdmodes','randomseed'});
compmodes   = keyval('compmodes', varargin);
niter       = keyval('niter', varargin);       if isempty(niter),        niter        = 2500;                    end
convcrit    = keyval('convcrit', varargin);    if isempty(convcrit),     convcrit     = 1e-6;                    end
randomseed  = keyval('randomseed', varargin);  if isempty(randomseed),   randomseed   = round(sum(clock.*1e6));  end
startval    = keyval('startval', varargin);
dispprefix  = keyval('dispprefix', varargin);  if isempty(dispprefix),   dispprefix   = [];                      end
holdmodes   = keyval('holdmodes', varargin);   if isempty(holdmodes),    holdmodes    = [];                      end

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

% set random seed (using clock as default)
rng(randomseed)

% set size and number of modes and compflg (and default compmodes)
nmode = ndims(dat);
smode = size(dat);
compflg = ~isreal(dat);
if isempty(compmodes),
  compmodes = ones(1,ndims(dat)) .* compflg;
end
modeind = 1:nmode; % for small speedup when array is small

% display input data
dimstring = num2str(smode(1));
for imode = 2:nmode
  dimstring = [dimstring 'x' num2str(smode(imode))];
end
disp([dispprefix 'data is complex array with dimensions ' dimstring])
disp([dispprefix 'a PARAFAC-model with ' num2str(ncomp) ' components will be estimated '])
disp([dispprefix 'maximum number of iterations = ' num2str(niter)])
disp([dispprefix 'convergence criterion = ' num2str(convcrit)])
disp([dispprefix 'random seed = ' num2str(randomseed)])

% throw errors for real/complex input with complex/real component matrices
if ~compflg && sum(compmodes) ~= 0
  error('complex input needed for some component matrices to be complex')
elseif compflg && sum(compmodes) == 0
  error('real-valued input needed for all component matrices to be real-valued')
end

% throw error when compmodes is not the same size as data
if nmode ~= length(compmodes)
  error('length of compmodes needs to be equal to ndims(dat)')
end

% ensure working with double precision
if ~isa(dat,'double')
  dat = double(dat);
end

% throw errors related to NaNs
if any(isnan(dat(:)))
  error('dat cannot contain NaNs')
end

% default holdmodes
if isempty(holdmodes)
  holdmodes = zeros(1,nmode);
end

% unfold dat's for computation during ALS
disp([dispprefix 'unfolding data into ' num2str(nmode) ' separate permutations'])
datorig = dat;
dat     = cell(1,nmode);
modeindex = 1:nmode;
for imode = 1:nmode
  permorder   = [imode setdiff(modeindex,imode)];
  reshapesize = [smode(imode) prod(smode(setdiff(modeindex,imode)))];
  dat{imode}  = reshape(permute(datorig,permorder),reshapesize);
end
datmain = dat{1}; % this is identical to dat{1}, put here for code transparency
clear datorig % no longer needed


% produce random start values from -1 to 1, if no start-values were provided
if isempty(startval)
  comp = cell(1,nmode);
  for imode = 1:nmode
    compsize    = [smode(imode) ncomp];
    if compmodes(imode)==1
      comp{imode} = complex((rand(compsize)*2)-1,(rand(compsize)*2)-1);
    elseif compmodes(imode)==0
      comp{imode} = rand(compsize);
    end
  end
else
  comp = startval;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   PARAFAC ALS   START  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The below code calculates an ALS estimate of the PARAFAC model based on the following structure:
% ((x) = Khatri-rao-bro product)
%
% Model (4-way):
% X = A*(D(x)C(x)B)' + E
%
% Estimate A by minimizing X - AZ'  (frobenius norm squared, x = unfolded)
% where Z = D(x)C(x)B
% ---> A = X*Z * inv(Z'Z)
%
%
%
%
%


% Set some  important variables
ssqdat     = sum(sum(abs(datmain).^2));
ssqres     = ssqdat;
prevssqres = 2 * ssqres;
iter       = 0;


% start main while loop of algorithm (updating component matrices)
disp([dispprefix 'starting ALS algorithm using QR-based fit estimation'])
while (abs((ssqres - prevssqres) / prevssqres) > convcrit) && (iter < niter) &&  (abs(ssqres) / ssqdat) > eps*1e3
  
  
  % Count iter
  iter = iter + 1;
  
  % Perform linear search every 4 iterations if relative ssqres increase is smaller than 10% (from the perspective of the previous iteration)
  if rem((iter-1),4) == 0 && ((prevssqres / ssqres) <= 1.05) && (((iter-1)^1/3) >= 2) && (iter < (niter-2))
    [comp] = linsearch(datmain,comp,prevcomp,iter-1,prevssqres,ssqres,dispprefix); % subfunction for performing linear search
  end
  
  % Set previous stuff (important for linear search in next iteration)
  prevcomp = comp;
  prevssqres  = ssqres;
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Update all component matrices per mode (this is the ALS routine)
  for imode = 1:nmode
    if holdmodes(imode)==0
      
      % set remaining mode, those to calculate the krb's over
      remmodes = modeind(modeind~=imode); % remaining modes
      
      
      % Calculate Z by a series of nested khatri-rao-bro products (using kr)
      Z = kr(comp(remmodes(end:-1:1)));
      
      
      % Old code, kept here for future use
      %       % Calculate Z by building tempZ over nested khatri-rao-bro products
      %       % Z is calculated by a nested series of khatri-rao-bro products
      %       % E.g. the Z for the first mode with 4 modes: Z = krbcomp(comp{4},krbcomp(comp{3},comp{2}));
      %       nremmodes = length(remmodes); % number of remaining modes
      %       nkrbcomps = nremmodes-1; % number of khatri-rao-bro products to calculate
      %       tempZ     = comp{remmodes(1)}; % second element of first krbcomp is always first remaining mode
      %       for ikrbcomp = 1:nkrbcomps % loop over nested calculations of khatri-rao-bro productsda
      %         tempZ = krbcomp(comp{remmodes(ikrbcomp+1)},tempZ);
      %       end
      %       Z = tempZ;
      %       % Update the component matrix for the current mode
      %       comp{imode} = dat{imode} * conj(Z) * inv(Z'*Z).'; % .' is equal to conj in this case
      
      
      
      % calculate ZctZ directly, which is faster than calculating Z first, especially when decomposing large arrays
      ZctZ = comp{remmodes(end)}'*comp{remmodes(end)};
      for iremmode = 1:numel(remmodes)-1
        ZctZ = ZctZ .* (comp{remmodes(iremmode)}'*comp{remmodes(iremmode)});
      end
      
      % Update the component matrix for the current mode
      % if complex, and currmode shouldn' be, take real of output that would otherwise be computed using catted real and imag parts, it's equivalent
      if compmodes(imode)==0 && compflg % computation below is identical (taking real of real data), but split up for code transparency
        %comp{imode} = real(dat{imode} * conj(Z)) * inv(real(ZctZ)).'; % real(Z'*X) = [Zre Zim]' * [Xre Xim];b  (.' is equal to conj in this case)
        comp{imode} = real(dat{imode} * conj(Z)) / real(ZctZ);
      else
        %comp{imode} = (dat{imode} * conj(Z)) * conj(inv(ZctZ));
        comp{imode} = (dat{imode} * conj(Z)) / conj(ZctZ);
      end
      
    end
  end % end of looping over modes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % normalize components, for scaling, phase and permutation indeterminancy
  comp = normalizecomp(comp,compmodes,0,dispprefix);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute fit with QR-decomposition of loadings
  % (much faster and ssqres approaches actual ssqres)
  rcomp = cell(1,nmode);
  % perform QR-decomposition of all component-matrices
  for imode = 1:nmode
    [dum rcomp{imode} ] = qr(comp{imode},0);
  end
  model     = calcmodel(rcomp); % subfunction for calculating model
  ssqmodel  = sum(sum(abs(model).^2));
  ssqres    = ssqdat - ssqmodel;
  expvar    = (ssqmodel / ssqdat) * 100;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  % Display results of current iteration
  disp([dispprefix 'iteration ' num2str(iter) ' - expvar: ' num2str(expvar,'%-2.1f')   '%  ssqres: ' num2str(ssqres)  '  ssqmodel: ' num2str(ssqmodel)])
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Algorithm stops
  % calculate final explained variance, not based on QR-decomposition
  if ~(abs((ssqres - prevssqres) / prevssqres) > convcrit) || (iter == niter)
    model     = calcmodel(comp); % subfunction for calculating model
    ssqres    = sum(sum(abs(datmain - model).^2));
    qrexpvar  = expvar;
    expvar    = 100 - ((ssqres / ssqdat) * 100);
    if abs(qrexpvar - expvar) > 1e-3
      error(['problem with QR-based acceleration: expvar = ' num2str(expvar,'%-2.4f') '%  QR-based expvar = ' num2str(qrexpvar,'%-2.4f') '%'])
    end
  end
  if ~(abs((ssqres - prevssqres) / prevssqres) > convcrit)
    disp([dispprefix 'convergence criterion of ' num2str(convcrit) ' reached in ' num2str(iter) ' iterations'])
    disp([dispprefix 'explained variance by model: ' num2str(expvar,'%-2.1f') '%'])
  elseif (iter == niter)
    disp([dispprefix 'maximum number of iterations = ' num2str(iter) ' reached'])
    disp([dispprefix 'explained variance by model: ' num2str(expvar,'%-2.1f') '%'])
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
end % end main while loop of algorithm (updating component matrices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   PARAFAC ALS   END    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   POSTPROCESSING START  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize components, for scaling, phase and permutation indeterminancy
comp = normalizecomp(comp,compmodes,1,dispprefix);

% Calculate Core Consistency Diagnostic if desired
if nargout > 5
  disp([dispprefix 'calculating tucker3 core array'])
  % Calculate ZtZ 
  % building it over nested kronecker products of ctc's (in 3-way case, ZtZ = kron(c3tc3,kron(c2tc2,c1tc1)) )
  tempZtZ    = comp{1}'*comp{1}; % second element of first kronprod is always first ctc
  for imode = 2:nmode % loop over nested calculations of kron products
    tempZtZ = kron(comp{imode}'*comp{imode},tempZtZ);
  end
  ZtZ = tempZtZ;
  % Calculate XtZ
  XtZ = zeros(ncomp^nmode,1);
  % determine component combinations
  modecompind = zeros(ncomp.^nmode,nmode);
  for imode = 1:nmode
    step1 = reshape(repmat((1:ncomp)',[1 max([1 ncomp.^(imode-1)])])',[ncomp*max([1 ncomp.^(imode-1)]) 1]);
    step2 = repmat(step1,[ncomp.^(nmode-imode) 1]);
    modecompind(:,imode) = step2;
  end
  % build XtZ per component combination
  for iccomb = 1:(ncomp.^nmode)
    tempckron = comp{1}(:,modecompind(iccomb,1)); % second element of first kronprod is always start at first mode
    for imode = 2:nmode % loop over nested calculations of kron products
      tempckron = kron(comp{imode}(:,modecompind(iccomb,imode)),tempckron);
    end
     XtZ(iccomb) = dat{1}(:)' * tempckron;
  end
  t3core = pinv(real(ZtZ)) * real(XtZ); % real(Z'*X) = [Zre Zim]' * [Xre Xim]; 
end

% Compute ssq per component and sort
ssqcomp = diag(comp{1}'*comp{1});
[ssqcomp,sortorder] = sort(ssqcomp,'descend');
for imode = 1:nmode % sort per mode
  comp{imode} = comp{imode}(:,sortorder);
end
disp([dispprefix 'components have been ordered according to magnitude'])


% Continue postprocessing
% set norm of mode 1 per loading vector to 1 and save scaling coeffient
for icomp = 1:ncomp
  % set mode1
  mode1     = comp{1}(:,icomp);
  mode1norm = norm(mode1,'fro');
  mode1     = mode1 ./ mode1norm;
  comp{1}(:,icomp)  = mode1;
  scaling{1}(icomp) = mode1norm;
end
disp([dispprefix 'first mode magnitude scaling coefficient removed'])
% set average magnitude weighted phase to 0 for the first complex loading vector per component and save scaling coefficient
compmodeindex = find(compmodes==1);
if sum(compmodes)~=0 % only do this if there is a complex mode
  for icomp = 1:ncomp
    % set mode1
    mode1      = comp{compmodeindex(1)}(:,icomp); % compmodeindex determined above
    meanangle  = angle(mean(mode1)); % angle of the mean weighted by magnitude
    phaseshift = exp(-1i*meanangle);
    mode1      = mode1 .* phaseshift;
    comp{compmodeindex(1)}(:,icomp) = mode1;
    scaling{2}(icomp) = phaseshift;
  end
  disp([dispprefix 'first complex mode has been phase shifted so average magnitude-weighted-phase = 0'])
end
disp([dispprefix 'post-processing finished'])


% Calculate Tucker's Congruency coefficient
tuckcongr = ones(ncomp,ncomp);
for imode = 1:nmode
  innprod = comp{imode}' * comp{imode};
  tuckcongr = tuckcongr .* innprod;
end
% remove upper triangle and diagonal, make vector, and take abs
tuckcongr = tril(tuckcongr,-1);
tuckcongr(tuckcongr==0) = [];
tuckcongr = abs(tuckcongr);
if ncomp==1 % if ncomp is one, tuckers congruence cannot be calculated
  tuckcongr = NaN;
end
% display warning if too high
if max(tuckcongr) >= 0.85
  disp([dispprefix 'Warning: some components are highly correlated, model might be degenerate'])
end

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
%%%%%%%%%%%%%%%%%%%%  Subfunction normalizecomp   %%%%%%%%%%%%
function comp = normalizecomp(comp,compmodes,dispflg,dispprefix)
% set stuff for subfunction internal use
ncomp = size(comp{1},2);
nmode = numel(comp);

% Magnitude postprocessing of components (make frobenius norm per loading vector norm = 1)
% make sure all explained variance is contained in first mode
for icomp = 1:ncomp
  % set mode1, which get's all the variance
  mode1 = comp{1}(:,icomp);
  % loop over the other modes
  for imode = 2:nmode
    currmode     = comp{imode}(:,icomp); % current mode of current compoment
    currmodenorm = norm(currmode,'fro');
    mode1        = mode1    .* currmodenorm;
    currmode     = currmode ./ currmodenorm;
    comp{imode}(:,icomp) = currmode;
  end
  % set mode1 back into original format
  comp{1}(:,icomp) = mode1;
end
if dispflg
  disp([dispprefix 'components have been magnitude normalized in all but the first mode'])
end

% Phase postprocessing of complex components (make average magnitude weigthed phase per complex loading vector 0)
% make sure the scaling of phases is contained in first complex mode
compmodeindex = find(compmodes==1);
if sum(compmodes)>=2 % only do this if there is more than one complex mode
  for icomp = 1:ncomp
    % set mode1 (the first complex mode), which get's the phase scaling
    mode1 = comp{compmodeindex(1)}(:,icomp);
    % loop over the other modes
    for imode = compmodeindex(2:end)
      currmode   = comp{imode}(:,icomp); % current mode of current compoment
      meanangle  = angle(mean(currmode)); % angle of the mean weighted by magnitude
      phaseshift = exp(-1i*meanangle);
      mode1      = mode1    ./ phaseshift;
      currmode   = currmode .* phaseshift;
      comp{imode}(:,icomp) = currmode;
    end
    % set mode1 back into original format
    comp{compmodeindex(1)}(:,icomp) = mode1;
  end
  if dispflg
    disp([dispprefix 'components have been phase shifted so average magnitude-weighted-phase = 0 with respect to the first complex mode'])
  end
end


% Apply sign convention, all real-valued modes are set to be 'mainly' positive, on a per component basis
% If a loading vector in a mode has a negative sum(value) then the entire loading vector is scaled by -1,
% the respective loading vector in the first mode is then also scaled by -1 (inv(-1) to be precise)
% set mode1 which get's all the scaling
mode1         = comp{1};
signsmode1    = ones(1,ncomp);
% loop over the other modes
for imode = 2:nmode
  if isreal(comp{imode}) % don't scale the complex modes, it's irrelevant
    currmode      = comp{imode};
    signscurrmode = sign(sum(currmode,1));
    signsmode1    = signsmode1 .* signscurrmode; % update signs to be used at mode1
    comp{imode}   = currmode * diag(signscurrmode);
  end
end
% Apply all scaling to mode1
comp{1} = mode1 * diag(signsmode1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction Calculate model     %%%%%%%%%%%%
function [model] = calcmodel(comp)

% Calculate model by building it over a series of nested khatri-rao-products
% E.g. the model for 4 modes: model = comp{1} * krbcomp(comp{4},krbcomp(comp{3},comp{2})).';
model = comp{1} * kr(comp(end:-1:2)).';

% Old code, kept here for future use
% % set several variables
% nmode    = length(comp); % number of modes
% nkrbcomps = nmode-2; % number of khatri-rao-bro products to calculate
% tempmodel = comp{2}; % second element of first krbcomp is always second mode
%
% % start series of khatri-rao-bro products
% for ikrbcomp = 1:nkrbcomps %
%   tempmodel = krbcomp(comp{ikrbcomp+2},tempmodel); %
% end
%
% % last step in calculating model
% model = comp{1} * tempmodel.';





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Subfunction Linear search ALS   %%%%%%%%%%%%%%%%%
function [newcomp] = linsearch(dat,comp,prevcomp,iter,prevssqres,ssqres,dispprefix)
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

% set delta and other things
delta = iter ^ 1/3;
nmode = length(comp);


% set difference of current loadings with previous loadings
dcomp = cell(1,nmode);
for imode = 1:nmode
  dcomp{imode} = comp{imode} - prevcomp{imode};
end


% calculate expected ssqres of delta iterations later
comppred = cell(1,nmode);
for imode = 1:nmode
  comppred{imode} = comp{imode} + delta*dcomp{imode};
end
model = calcmodel(comppred); % subfunction for calculating model
ssqresdelta = sum(sum(abs(dat - model).^2));


% set ssqres and ssq collection variable for picking best delta at the end
ssqaccel{1} = [0 ssqres];
ssqaccel{2} = [delta ssqresdelta];

% decrease delta when expected ssqres is larger than current ssqres
while (ssqresdelta > ssqres) && (delta >= 1)
  delta = delta * 0.5;
  if (delta >= 1)
    % calculate expected ssqres of delta iterations later
    comppred = cell(1,nmode);
    for imode = 1:nmode
      comppred{imode} = comp{imode} + delta*dcomp{imode};
    end
    model = calcmodel(comppred); % subfunction for calculating model
    ssqresdelta = sum(sum(abs(dat - model).^2));
    ssqaccel{end+1} = [delta ssqresdelta];
    % debugging  disp(['delta ' num2str(delta)])
    % debugging  disp(['ssqres   ' num2str(ssqresdelta)])
  else
    delta = delta * 2; % undo previous delta decrease for next while loop
    break
  end
end

% increase delta until ssqres no longer get's better, but only if the extrapolated ssqres is at least 4 times the previous ssqres increase
ssqresdeltaold = ssqres;
ssqresdeltanew = ssqresdelta;
prevssqresinc = prevssqres - ssqres;
while (ssqresdeltanew < ssqresdeltaold) && (((ssqresdeltaold - ssqresdeltanew) / prevssqresinc) > 2)
  delta = delta *1.5;
  % calculate expected ssqres of delta iterations later
  comppred = cell(1,nmode);
  for imode = 1:nmode
    comppred{imode} = comp{imode} + delta*dcomp{imode};
  end
  model = calcmodel(comppred); % subfunction for calculating model
  ssqresdeltaold = ssqresdeltanew;
  ssqresdeltanew = sum(sum(abs(dat - model).^2));
  % debugging  disp(['inc loop done  ' num2str(delta) '  ' num2str(ssqresdeltanew)])
  if (ssqresdeltanew < ssqresdeltaold)
    ssqaccel{end+1} = [delta ssqresdeltanew];
    % debugging  disp(['deltainc ' num2str(delta)])
    % debugging  disp(['ssqresinc   ' num2str(ssqresdeltanew)])
  end
end


% pick best delta by picking lowest ssqres
ssqaccel = sortrows(vertcat(ssqaccel{:}),2);
delta  = ssqaccel(1,1);
deltassqres = ssqaccel(1,2);


% create new component loading matrices
newcomp = cell(1,nmode);
for imode = 1:nmode
  newcomp{imode} = comp{imode} + delta*dcomp{imode};
end
disp([dispprefix 'acceleration performed with delta = ' num2str(delta,'%-8.2f') ' and ssqres = '  num2str(deltassqres) ' ('  num2str(length(ssqaccel(:,1))-1) ' ssq calcs)'])
















