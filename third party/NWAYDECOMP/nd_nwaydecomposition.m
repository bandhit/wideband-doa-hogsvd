function [nwaycomp] = nd_nwaydecomposition(cfg,data)

% ND_NWAYDECOMPOSITION Decomposes an N-way array of electrophysiological data into components.
%
% Please first read the README, acompanying this toolbox.
%
% A numerical N-way array needs to be present, and indicated in the cfg.
%
% Supported N-way models are: SPACE-FSP, SPACE-time, PARAFAC (modified), PARAFAC2 (modified)
% This function can be used to decompose numerical arrays, which are complex-valued or real-valued, into
% components. How a component is modeled and what it represents, as well as the required structure of the
% numerical array, depends on the model. This functions does two important things:
% 1) it takes care of randomly initializing each model
% 2) it can be used to estimate the number of components to extract
%
% Ad1:
% Each of the supported models needs to be randomly initialized multiple times in order to avoid local minima
% of the loss functions of each model, and to avoid degenerate solutions. Once it converges to the same
% solution from multiple random starting points, it can be assumed the global minimum is reached.
%
% Ad2:
% The number of components to extract from the numerical array needs to be determined empirically (similar to ICA).
% This is the case for all of the supported models. This function can be used for this purpose using 4 different
% strategies:
% 1) split-reliablility of the array. This strategy increasse the number of components until a reliability criterion
%    is no longer reached. This criterion is based on a statistic that assesses the similarity between components of the full data
%    and of componets of splits of the data (e.g. sets of trials). The separate splits of the data need to be given in a
%    separate field in the data, next to the full N-way array. (The splits can have dimensions that are different from the full array)
% 2) core-consistency diagnostic. This approach uses a statistic which can be viewed as a measures of noise being modelled
%    and, as such, is as an indication of whether the model with a certain number of components is still appropriate
%    (mostly appropriate for PARAFAC). This ranges from 0 to 1 (perfect)
% 3) minium increase in explained variance. A simple procedure that increases the number of components until the new component
%    no longer results in a certain increase of % explained variance.
% 4) degeneracy. This procedure keeps on increasing the number of components until some become degenerate. This uses a statistic
%    denoted as the Tucker's congruency coefficient, which ranges from 0 to 1 (fully degenerate).
%
% The settings of the above strategies are specified in the cfg.
%
%
% Models SPACE-time and SPACE-FSP are fully described in the following publication.
% Please cite when either of them are used:
%    van der Meij R, Jacobs J, Maris E (2015). Uncovering phase-coupled oscillatory networks in
% 	      electrophysiological data. Human Brain Mapping
% In the case of applying SPACE to extracranial recordings, please also cite the second reference paper:
%    van der Meij R, van Ede F, Maris E (2016). Rhythmic Components in Extracranial Brain
%         Signals Reveal Multifaceted Task Modulation of Overlapping Neuronal Activity. PLOS One
%
% The PARAFAC and PARAFAC2 models are modified such that parameter matrices can be real-valued when the input array is complex-valued.
% For additional info on the models and the split-half procedure, see the following publication, please cite when either of
% the models are used:
%    van der Meij, R., Kahana, M. J., and Maris, E. (2012). Phase-amplitude coupling in
%        human ECoG is spatially distributed and phase diverse. Journal of Neuroscience,
%        32(1), 111-123.
%
% For additional information of the models, please also see their low-level functions. And, of course,
% the README acompanying this toolbox.
%
%
%
%
% Use as
%   [nwaycomp] = nd_nwaydecomposition(cfg, data)
%
% The input data can be any type of FieldTrip-style data structure, as long as the field containing the data
% to be decomposed contains a numerical array. For the SPACE models, the input needs to be Fourier coefficients.
% These Fourier coefficients can be provided in 2 ways. They can either come from (1) custom code together with dimord
% of 'chan_freq_epoch_tap' (with the Fourier coefficients following this dimensionality; nepoch can be 1),
% or (2) they can can come as output from ft_freqanalysis (cfg.output = 'fourier' or 'powandcsd'; DPSS tapering is supported).
% In the case of the later, the 'rpt' (trial) dimension will be used as 'epochs' in SPACE terminology. The time dimension will be
% used as 'tapers' in SPACE terminology. If multiple tapers are present per time-point, these will be handled accordingly.
% Additionally, if you used method = 'mtmconvol', and frequency-dependent window-lengths, it is highly recommended to supply
% cfg.fsample, containing the sampling rate of the data in Hz.
%
%
%   cfg.model                = 'parafac', 'parafac2', 'parafac2cp', 'spacetime', 'spacefsp'
%   cfg.datparam             = string, containing field name of data to be decomposed (must be a numerical array)
%   cfg.randstart            = 'no' or number indicating amount of random starts for final decomposition (default = 50)
%   cfg.numiter              = number of iterations to perform (default = 2500)
%   cfg.convcrit             = number, convergence criterion (default = 1e-8)
%   cfg.degencrit            = number, degeneracy criterion (default = 0.7)
%   cfg.ncomp                = number of components to extract
%   cfg.ncompest             = 'no', 'splitrel', 'corcondiag', 'minexpvarinc', or 'degeneracy' (default = 'no'), see below (FIXME: complexicity of minexpvarinc and degeneracy should same as others)
%   cfg.ncompestrandstart    = 'no' or number indicating amount of random starts for estimating component number (default = cfg.randstart)
%   cfg.ncompeststart        = number, starting number of components to try to extract (default = 1) (used in splitrel/corcondiag)
%   cfg.ncompestend          = number, maximum number of components to try to extract (default = 50) (used in splitrel/corcondiag)
%   cfg.ncompeststep         = number, forward stepsize in ncomp estimation (default = 1) (backward is always 1; used in splitrel/corcondiag)
%   cfg.ncompestsrdatparam   = (for 'splitrel'): string containing field-name of partitioned data. Data should be kept in 1xN cell-array, each split in one cell
%                              when using SPACE, one can also specify 'oddeven' as cfg.ncompestsrdatparam. In this case the data will be partioned using odd/even trials/epochs
%   cfg.ncompestsrcritval    = (for 'splitrel'): scalar, or 1Xnparam vector, critical value to use for selecting number of components using splif-reliability (default = 0.7)
%                              when using SPACE, the default for the trial profile and between-component coherency (in case Dmode = identity) is set to 0
%   cfg.ncompestvarinc       = (for 'minexpvarinc'): minimal required increase in explained variance when increasing number of compononents by cfg.ncompeststep
%   cfg.ncompestcorconval    = (for 'corcondiag'): minimum value of the core consistency diagnostic for increasing the number of components, between 0 and 1 (default is 0.7)
%
%      Algorithm specific options:
%        PARAFAC/2(CP)
%   cfg.complexdims          = vector of 0's and 1's with length equal to number of dimensions in data, indicating dimensions to keep complex
%                              (default = 0/1 for each dimension if real/complex)
%        PARAFAC2(CP)
%   cfg.specialdims          = vector with length equal to ndims(dat) with 1, 2 and 3
%                              indicating special modes: 1: outer dim of inner-product                     (i.e. the utility dims)   (ndim-2 dims must be this)
%                                                        2: inner dim of inner-product                     (i.e. the incomplete dim) (only one allowed)
%                                                        3: dim over which inner-products will be computed (i.e. the estimating dim) (only one allowed)
%        SPACEFSP/SPACETIME
%   cfg.Dmode                = string, 'identity', 'kdepcomplex', type of D to estimate/use. Default = 'identity', 'kdepcomplex' is not advised.
%
%
%      -Using distributed computing to run random starts in parallel-
%   cfg.distcomp.system          = string, distributed computing system to use, 'torque' or 'matlabpct' ('torque' requires the qsub FieldTrip module on path, 'matlabpct' implementation is via parfor)
%   cfg.distcomp.timreq          = number, (torque only) maximum time requirement in seconds of a random start (default = 60*60*24*1 (1 days))
%   cfg.distcomp.memreq          = number, (torque only) maximum memory requirement in bytes of a random start (default is computed)
%   cfg.distcomp.inputsaveprefix = string, (torque only) path/filename prefix for temporarily saving input data with a random name (default, saving is determined by the queue system)
%   cfg.distcomp.matlabcmd       = string, (torque only) command to execute matlab (e.g. '/usr/local/MATLAB/R2012b/bin/matlab') (default = 'matlab')
%   cfg.distcomp.torquequeue     = string, (torque only) name of Torque queue to submit to (default = 'batch')
%   cfg.distcomp.torquestack     = number, (torque only) number of random initializations to stack in one Torque job (default is 1, i.e. no stacking)
%   cfg.distcomp.mpctpoolsize    = number, (matlabpct only) number of workers to use (default is determined by matlab)
%   cfg.distcomp.mpctcluster     = Cluster object, (matlabpct only) Cluster object specifying PCT Cluster profile/parameters, see matlab help PARCLUSTER
%
%
%       CFG.NCOMPEST - Methods for determining the number of components to extract
%             splitrel: Determines the number of *reliable* components to extract. Extract components from the full data and N splits of the data (e.g. odd/even trials), and judge
%                         similiarity (coef. 0<->1) between components from the main data and each of the splits. If similarity of each parameter surpasses its criterion for each of
%                         the splits, increase the number of components. Otherwise decrease till the criterion is satisfied.
%                         The criterion is set per parameter using cfg.ncompestsrcritval. Set the criterion to zero to ignore parameters.
%                         Applicable and advised for all models. See cfg.ncompestsrcritval & cfg.ncompestsrdatparam & cfg.ncompeststart/end/step
%                         To obtain a split-reliability estimate for a fixed number of componenents, set cfg.ncompeststart/end to the same number, and set cfg.ncompestsrcritval to NaN for the parameters that
%                         should determine the spilt-reliability coefficient (0 for the rest).
%                         (see any of the three reference paper, and Bro 1998, Multi-way Analysis in the Food Industry.)
%           corcondiag: Extract components, and compute the Core Consistency Diagnostic (coef. 0<->1). This coefficient indicates whether the components reflect the N-way
%                         structure of the data, or reflects noise. If the coefficients is lower than the criterion, decrease the number of components, otherwise, increase the number.
%                         Not advised for anything but PARAFAC. See cfg.ncompestcorconval & cfg.ncompeststart/end/step
%                         (see Bro 1998, Multi-way Analysis in the Food Industry.)
%         minexpvarinc: Increase the number of components until the increase of explained variance w.r.t. the previous number of components is smaller than cfg.ncompestvarinc
%                         See cfg.ncompeststart/end. Stepsize is always 1.
%           degeneracy: Increase the number of components until the components become degenerate in at least one of the parameters.
%                         Not advised for anything but PARAFAC. See cfg.degencrit & cfg.ncompeststart/end. Stepsize is always 1.
%
%
%         Output nwaycomp:
%             label: cell-array containing labels for each channel
%            dimord: string, dimord of input data
%              comp: cell-array, each component consists of 1 cell, each of these consist of 1 cell per dimension/parameter
%            expvar: number, variance explained by the model
%         tuckcongr: vector, Tucker's congruence coefficient per component
%           scaling: vector orcell-array, dep. on model, containing magnitude/phase scaling coefficients
%               cfg: input cfg and its history
%
%    Possible additional output fields:
%            t3core: a Tucker3 model core-array (vectorized) (not possible for all models)
%        randomstat: structure containing statistics of random estimation of final decomposition (if randomly started)
%      splitrelstat: structure containing statistics for split-half component number estimation procedure
%    corcondiagstat: structure containing statistics for corcondiag component number estimation procedure
%
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.

%  Undocumented options:
% (experimental)
% cfg.distcomp.system          = 'p2p'
% cfg.distcomp.p2presubdel     = scalar, resubmission delay for p2p in seconds (default = 60*60*24*3 (3 days))  (for p2p)
% cfg.distcomp.qsuboptions     = string, (torque only) additional options command-line options for qsub specified as a string
% cfg.ncompestsrcritjudge      = 'meanoversplitscons/lib' or 'minofsplitscons/lib', lib = picking SRC form best combination of randinits, cons = picking SRC of randinits with highest expvar
% cfg.checkpointpath           = path to use for checkpointing, save entire workspace minus data to path with unique ID, useful for dist. comp. jobs getting cancelled etc
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
%ft_preamble debug % ft_preamble_debug currently leads to qsubfeval saving a copy of the input data
ft_preamble loadvar data

%%% backwards compatability per August 2016 for oldsplithalf
% first, check for renamed splitrel options
if isfield(cfg,'ncompest') && strcmp(cfg.ncompest,'splithalf')
  warning(['You are using the old style split-reliability. Consider using the new method (see documentation), which computes '...
    ' a reliability coefficient between components of the full data, and components of each split (which can be more than one split).'...
    ' When doing so, the output fields will be slightly different. '])
  oldsplithalf = true;
else
  oldsplithalf = false;
end
cfg = ft_checkconfig(cfg, 'renamedval',  {'ncompest', 'splithalf', 'splitrel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'ncompestshdatparam', 'ncompestsrdatparam'});
cfg = ft_checkconfig(cfg, 'renamed',     {'ncompestshcritval', 'ncompestsrcritval'});
% check for only two splits
if oldsplithalf && ~any(strcmp(cfg.ncompestsrdatparam,{'oddeven','oddevenavg'})) && numel(data.(cfg.ncompestsrdatparam))~=2
  error('Old style split-reliability can only be used with two splits')
end
%%% backwards compatability per August 2016 for oldsplithalf


% Set defaults
cfg.model               = ft_getopt(cfg, 'model',                  []);
cfg.datparam            = ft_getopt(cfg, 'datparam',               []);
cfg.complexdims         = ft_getopt(cfg, 'complexdims',            []);
cfg.randstart           = ft_getopt(cfg, 'randstart',              50);
cfg.numiter             = ft_getopt(cfg, 'numiter',                2500);
cfg.convcrit            = ft_getopt(cfg, 'convcrit',               1e-8);
cfg.degencrit           = ft_getopt(cfg, 'degencrit',              0.7);
cfg.ncomp               = ft_getopt(cfg, 'ncomp',                  []);
cfg.ncompest            = ft_getopt(cfg, 'ncompest',               'no');
cfg.ncompestrandstart   = ft_getopt(cfg, 'ncompestrandstart',      cfg.randstart);
cfg.ncompeststart       = ft_getopt(cfg, 'ncompeststart',          1);
cfg.ncompestend         = ft_getopt(cfg, 'ncompestend',            50);
cfg.ncompeststep        = ft_getopt(cfg, 'ncompeststep',           1);
cfg.ncompestsrdatparam  = ft_getopt(cfg, 'ncompestsrdatparam',     []);
cfg.ncompestsrcritval   = ft_getopt(cfg, 'ncompestsrcritval',      0.7); % expanded to all paramameters later
cfg.ncompestsrcritjudge = ft_getopt(cfg, 'ncompestsrcritjudge',    'minofsplitslib');
cfg.specialdims         = ft_getopt(cfg, 'specialdims',            []); % parafac2 specific
cfg.ncompestvarinc      = ft_getopt(cfg, 'ncompestvarinc',         []);
cfg.Dmode               = ft_getopt(cfg, 'Dmode',                  'identity'); %  spacefsp/spacetime specific
cfg.ncompestcorconval   = ft_getopt(cfg, 'ncompestcorconval',      0.7);
cfg.t3core              = ft_getopt(cfg, 't3core',                 'no');
cfg.checkpointpath      = ft_getopt(cfg, 'checkpointpath',         []);

% set distributed computing random starting defaults and throw errors
cfg.distcomp                  = ft_getopt(cfg, 'distcomp',                    []);
cfg.distcomp.system           = ft_getopt(cfg.distcomp, 'system',             []);
cfg.distcomp.memreq           = ft_getopt(cfg.distcomp, 'memreq',             []);
cfg.distcomp.timreq           = ft_getopt(cfg.distcomp, 'timreq',             60*60*24*1);
cfg.distcomp.inputsaveprefix  = ft_getopt(cfg.distcomp, 'inputsaveprefix',    []); % i.e. current dir
cfg.distcomp.matlabcmd        = ft_getopt(cfg.distcomp, 'matlabcmd',          'matlab'); % i.e. current dir
cfg.distcomp.torquequeue      = ft_getopt(cfg.distcomp, 'torquequeue',        'batch');
cfg.distcomp.torquestack      = ft_getopt(cfg.distcomp, 'torquestack',        1);
cfg.distcomp.p2presubdel      = ft_getopt(cfg.distcomp, 'p2presubdel',        60*60*24*3);
cfg.distcomp.qsuboptions      = ft_getopt(cfg.distcomp, 'qsuboptions',        []);
cfg.distcomp.mpctpoolsize     = ft_getopt(cfg.distcomp, 'mpctpoolsize',       []);
cfg.distcomp.mpctcluster      = ft_getopt(cfg.distcomp, 'mpctcluster',        []);
if strcmp(cfg.distcomp.system,'p2p') && isempty(cfg.distcomp.p2presubdel)
  error('need to specifiy cfg.distcomp.p2presubdel')
end


% check essential cfg options
if isempty(cfg.model)
  error('please specify cfg.model')
else
  cfg.model = lower(cfg.model); % make sure its lowercase
end
if isempty(cfg.datparam)
  error('you need to specify cfg.datparam')
end

% make a specific check for presence cfg.trials and cfg.channel
if isfield(cfg,'trials') || isfield(cfg,'channel')
  error('cfg.trials and cfg.channel are not supported')
end

% Check the data provided and get dimensions of data for for error checking below (parse filename if data.(cfg.datparam) is not data)
if ~isnumeric(data.(cfg.datparam))
  filevars = whos('-file', data.(cfg.datparam));
  % check file structure
  if numel(filevars)>1
    error('data filename can only contain a single variable')
  end
  if ~any(strcmp(filevars.class,{'single','double'}))
    error('array in data filename needs to be a numeric single or double array')
  end
  ndimsdat     = numel(filevars.size);
  datprecision = filevars.class;
  datcomplex   = filevars.complex;
else
  if ~isfloat(data.(cfg.datparam))
    error('field specified by cfg.datparam needs to contain a numeric single or double array')
  end
  ndimsdat     = ndims(data.(cfg.datparam));
  datprecision = class(data.(cfg.datparam));
  datcomplex   = ~isreal(data.(cfg.datparam));
end

% default options for parafac/parafac2/parafac2cp
if strncmp(cfg.model,'parafac',7) && isempty(cfg.complexdims)
  cfg.complexdims = ones(1,ndimsdat) .* datcomplex;
end

% Make sure a dimord is present (in case one uses this outside of FT)
if ~isfield(data,'dimord')
  error(['Input structure needs to contain a dimord field. '...
    'This field identifies the dimensions in the most important data-containing field. See FieldTrip wiki for further details'])
end


% Ncomp and Ncompest related errors
if (isempty(cfg.ncomp) && strcmp(cfg.ncompest,'no'))
  error('you either need to estimate the number of components or set cfg.ncomp')
end
if ~isempty(cfg.ncomp) && ~strcmp(cfg.ncompest,'no')
  warning('cfg.ncomp cannot be used when estimating number of components, cfg.ncomp is ignored')
  cfg.ncomp = [];
end
if strcmp(cfg.ncompest,'splitrel') && isempty(cfg.ncompestsrdatparam)
  error('you need to specify cfg.ncompestsrdatparam')
end
if strcmp(cfg.ncompest,'minexpvarinc') && isempty(cfg.ncompestvarinc)
  error('must set cfg.ncompestvarinc, the minimal required increase in explained variance')
end
if strncmp(cfg.model,'parafac2',8) && strcmp(cfg.ncompest,'corcondiag')
  error('at the moment corcondiag cannot be used when cfg.model is parafac2/cp')
end
if numel(cfg.ncompeststart) ~= 1 || numel(cfg.ncompestend) ~= 1 || numel(cfg.ncompeststep) ~= 1
  error('improper cfg.ncompestXXX')
end

% check splitrel and ncompestsrcritval
if strcmp(cfg.ncompest,'splitrel')
  switch cfg.model
    case {'parafac','parafac2','parafac2cp','spacetime','spacefsp'}
      % implemented
    otherwise
      error('model not supported')
  end
end
if strcmp(cfg.ncompest,'splitrel')
  if (numel(cfg.ncompestsrcritval)==1) % expand using defaults for SPACE
    if strncmp(cfg.model,'parafac',7)
      cfg.ncompestsrcritval = repmat(cfg.ncompestsrcritval,[1 ndimsdat]);
    elseif strcmp(cfg.model,'spacetime') || strcmp(cfg.model,'spacefsp')
      cfg.ncompestsrcritval = repmat(cfg.ncompestsrcritval,[1 5]);
      cfg.ncompestsrcritval(3) = 0; % set the trial profile to zero by default
      if strcmp(cfg.Dmode,'identity')
        cfg.ncompestsrcritval(5) = 0; % set D to zero by default
      end
    end
  else
    if (strncmp(cfg.model,'parafac',7) && (numel(cfg.ncompestsrcritval)~=ndimsdat))... % FIXME: likley should contain check for PARAFAC2
        || (strcmp(cfg.model,'spacetime') && (numel(cfg.ncompestsrcritval)~=5))...
        || (strcmp(cfg.model,'spacefsp') && (numel(cfg.ncompestsrcritval)~=5))
      error('improper size of cfg.ncompestsrcritval')
    end
  end
  if all(cfg.ncompestsrcritval==0)
    error('cfg.ncompestsrcritval cannot be zero for all parameters')
  end
end

% distcomp check
if ~isempty(cfg.distcomp.system)
  switch cfg.distcomp.system
    case 'torque'
      if ~ft_hastoolbox('qsub')
        error('FieldTrip qsub module needs to be on the path in order to used Torque distributed computing')
      end
      
    case 'matlabpct'
      if ~ft_hastoolbox('distcomp')
        error('MATLAB Parallel Computing Toolbox needs to be present in order to use it for distributed computing')
      end
      
    case 'p2p'
      % do nothing
      
    otherwise
      error('specified distributed computing system not supported')
  end
end

% disp progress
disp(['N-way decomposition using ' upper(cfg.model) ' will be performed on ''' cfg.datparam ''''])
if ~strcmp(cfg.ncompest,'no')
  disp(['number of components to extract will be determined using method ''' cfg.ncompest ''''])
  disp([num2str(cfg.ncompestrandstart) ' random initializations will be used during component number estimation' ])
  disp([num2str(cfg.randstart) ' random initializations will be used for the final decomposition' ])
else
  disp([num2str(cfg.ncomp) ' components will be extracted'])
  disp([num2str(cfg.randstart) ' random initializations will be used ' ])
end
disp(['decompositions will stop once the convergence criterion of ' num2str(cfg.convcrit) ' has been reached, or after ' num2str(cfg.numiter) ' iterations' ])
disp(['degeneracy criterion that will be used is a Tuckers congruence coefficient of ' num2str(cfg.degencrit)])


% Model-specific errors
% parafac2/parafac2cp
if strncmp(cfg.model,'parafac2',8) && ndimsdat ~= numel(cfg.specialdims)
  error('length of cfg.specialdims should be equal to number of dimensions in data')
end
if strncmp(cfg.model,'parafac2',8) && ~isfield(data,'ssqdatnoncp')
  error('sneaky trick required')
end
% parafac/parafac2/parafac2cp
if strncmp(cfg.model,'parafac',7) && ndimsdat ~= numel(cfg.complexdims)
  error('length of cfg.complexdims should be equal to number of dimensions in data')
end
if strcmp(cfg.ncompest,'splitrel')
  if strncmp(cfg.model,'parafac',7) && (numel(cfg.complexdims) ~= numel(cfg.ncompestsrcritval))
    error('length of cfg.ncompestsrcritval should be equal to the number of dimensions in the data')
  end
end


% Handle input for SPACE models
% the below code also applies a trick to dramatically reduce memory in some case
if any(strcmp(cfg.model,{'spacefsp','spacetime'}))
  if strcmp(ft_datatype(data),'freq') && ~strcmp(data.dimord,'chan_freq_epoch_tap') % ft_datatype is a bit eager, add check for dimord
    % Handle output from ft_freqanalysis
    if ~any(strcmp(cfg.datparam,{'fourierspctrm','crsspctrm'}))
      error('cfg.datparam should be either ''fourierspctrm'' or ''crsspctrm'' when input is output from ft_freqanalysis')
    end
    % disp progress
    disp(['input for SPACE is the output of ft_freqanalysis, reorganizing ' cfg.datparam ' to have dimensionality ''chan_freq_epoch_tap'''])
    if ~isnumeric(data.(cfg.datparam))
      error('specifying data field as filename is only possible with manually constructed chan_freq_epoch_tap')
    end
    % throw error based on specestmethod, this has to do with current nsample normalization of Fourier coefficients (not important for mtmfft)
    specestmethod = data.cfg.method;
    if ~any(strcmp(specestmethod,{'mtmfft','mtmconvol'}))
      error(['distortion-free scaling of frequency-dependent time-windows lengths over frequency using ft_freqanalysis with method = ' specestmethod ' is not guaranteed'])
    end
    if strcmp(specestmethod,'mtmconvol') && ~isfield(cfg,'fsample')
      warning(['Your input resulted from ft_freqanalysis with method = mtmconvol. If you''re also using frequency-dependent window-lengths ' ...
        'it is highly recommended to supply cfg.fsample, containing the sampling rate of the data in Hz, to correct for ' ...
        'frequency-dependent distortions of power'])
    else
      % disp progress
      disp('applying correction for double scaling in ft_freqanalysis with method = mtmconvol')
    end
    % failsafe error for if this ever becomes supported in ft_freqananalysis (which it shouldn't)
    if strcmp(cfg.datparam,'fourierspctrm') && (strcmp(data.cfg.keeptrials,'no') || strcmp(data.cfg.keeptrials,'no'))
      error('fourierspctrm must have been computed using keeptrials and keeptapers = yes')
    end
    % failsafe error for if this ever becomes supported in ft_freqananalysis (which it shouldn't)
    if strcmp(cfg.datparam,'crsspctrm') && strcmp(data.cfg.keeptrials,'yes') && strcmp(data.cfg.keeptrials,'no')
      error('crsspctrm must have been computed using keeptrials = yes and keeptapers = no') % this can likely be detected below if it gets implemented
    end
    % failsafe error for if this ever becomes supported in ft_freqananalysis (which it might)
    if any(any(diff(data.cumtapcnt,1,2)))
      error(['Variable number of tapers over frequency is not supported using output from ft_freqanalysis.'...
        'In order to do this using a custom ''chan_freq_epoch_tap'' see the tutorial on rhythmic components'...
        'and the code below this error message.'])
    end
    % make sure channel cmb representation of data is full in case of crssspctrm
    if strcmp(cfg.datparam,'crssspctrm') % treat crsspctrm as an exception to fourierspctrm in the below
      data = ft_checkdata(data,'cmbrepresentation','full');
    end
    % if no trials are present, add trials as singular dimension to ensure SPACE input is 4-way
    if ~strncmp(data.dimord,'rpt',3)
      data.dimord = ['rpt_' data.dimord];
      data.(cfg.datparam) = permute(data.(cfg.datparam),[ndimsdat+1 1:ndimsdat]);
      % disp progress
      disp('input data no longer has trial dimension, adding singleton dimension')
    end
    % set
    ntrial = size(data.cumtapcnt,1);
    nfreq  = size(data.cumtapcnt,2);
    nchan  = numel(data.label);
    ntaper = nchan;
    dat    = complex(NaN(nchan,nfreq,ntrial,ntaper,datprecision),NaN(nchan,nfreq,ntrial,ntaper,datprecision));
    tapcnt = data.cumtapcnt(:,1); % explicitly only use ntap of first freq due to above
    % construct new dat
    for itrial = 1:ntrial
      trialtapind = sum(tapcnt(1:itrial-1))+1:sum(tapcnt(1:itrial-1))+tapcnt(itrial); % ow
      for ifreq = 1:nfreq
        % first, select data and create csd
        if strcmp(cfg.datparam,'fourierspctrm')
          currfour = permute(data.(cfg.datparam)(trialtapind,:,ifreq,:),[2 1 3 4]); % will work regardless of time dim presence/absence
          currfour = currfour(:,:); % this unfolds the dimensions other than chan, will work regardless of time dim presence/absence
          % get rid of NaNs (should always be the same over channels)
          currfour(:,isnan(currfour(1,:))) = [];
          % ensure double precision for calculations below
          currfour = double(currfour);
          %%%%%%
          % UNDO double scaling in ft_specest_mtmconvol
          % There is currently a double scaling applied in ft_specest_mtmconvol. This will likely not hurt other analyses,
          % but causes a distortion over frequencies that SPACE is sensitive for.
          % This scaling is only dependent on frequency for fourierspctrm
          % This scaling requires cfg.fsample to be given (... :( )
          if strcmp(data.cfg.method,'mtmconvol') && isfield(cfg,'fsample')
            % reconstruct tinwinsample
            t_ftimwin     = ft_findcfg(data.cfg,'t_ftimwin');
            timwinnsample = round(t_ftimwin(ifreq) .* cfg.fsample);
            % undo additional the scaling by sqrt(2./ timwinnsample)
            currfour      = currfour ./ sqrt(2 ./ timwinnsample);
          end
          %%%%%%
          % obtain count of time-points and tapers
          currntimetap = size(currfour,2);
          % compute the csd
          csd = currfour*currfour';
          % correct for number of time-points and tapers (if fourierspctrm, individual tapers and time-points (if mtmconvol) are not aggegrated, enforced above)
          % this step is CRUCIAL if we want to interpret the loadings of the trial profile
          csd = csd ./ currntimetap;
        else
          currcsd = permute(data.(cfg.datparam)(itrial,:,:,ifreq,:),[2 3 5 1 4]); % will work regardless of time dim presence/absence
          % get rid of NaNs (should always be the same over channel-pairs)
          currcsd(:,:,isnan(squeeze(currcsd(1,1,:)))) = [];
          % ensure double precision for calculations below
          currcsd = double(currcsd);
          %%%%%%
          % UNDO double scaling in ft_specest_mtmconvol
          % There is currently a double scaling applied in ft_specest_mtmconvol. This will likely not hurt other analyses,
          % but causes a distortion over frequencies that SPACE is sensitive for.
          % This scaling is only dependent on frequency for fourierspctrm
          % This scaling requires cfg.fsample to be given (... :( )
          if strcmp(data.cfg.method,'mtmconvol') && isfield(cfg,'fsample')
            % reconstruct tinwinsample
            t_ftimwin     = ft_findcfg(out.cfg,'t_ftimwin');
            timwinnsample = round(t_ftimwin(ifreq) .* cfg.fsample);
            % undo additional the scaling by (sqrt(2./ timwinnsample)).^2 = (2./ timwinnsample)
            % the scaling is performed on the level of individual taper-specific Fourier coefficients, and currcsd contains the summed cross-products of these
            currcsd       = currcsd ./ (2 ./ timwinnsample);
          end
          %%%%%%
          % obtain count of time-points (tapers are never kept if crsspctrm, enforced above)
          currntime = size(currcsd,3);
          % obtain csd
          csd = sum(currcsd,3);
          % correct for number of time-points (if crsspctrm, time-points (if mtmconvol) are not aggegrated, but tapers are averaged, enforced above)
          % this step is CRUCIAL if we want to interpret the loadings of the trial profile
          csd = csd ./  currntime; % this is a rather circumstantial way to do this, but is such to keep it similar to above
        end
        %%%%%%%%%
        % Reduce SPACE memory load and computation time by replacing each chan_taper matrix by the
        % Eigenvectors of its chan_chan cross-products weighted by sqrt(Eigenvalues).
        % This is possible because (1) SPACE only uses the cross-products of the chan_taper matrices
        % (i.e. the frequency- and trial-specific CSD) and (2) the Eigendecomposition of a symmetric
        % matrix A is A = VLV'.
        % As such, VL^.5 has the same cross-products as the original chan_taper matrix.
        [V L] = eig(csd);
        L     = diag(L);
        tol   = max(size(csd))*eps(max(L)); % compute tol using matlabs default
        zeroL = L<tol;
        eigweigth = V(:,~zeroL)*diag(sqrt(L(~zeroL)));
        % positive semi-definite failsafe
        if any(L<-tol) || any(~isreal(L))
          error('csd not positive semidefinite')
        end
        %%%%%%%%%
        % save in dat
        currm = size(eigweigth,2);
        dat(:,ifreq,itrial,1:currm) = eigweigth;
      end
    end
    % trim dat (actual ntaper can be lower than ntaper, which depends on the data, hence the need for trimming)
    notnan = logical(squeeze(sum(sum(~isnan(squeeze(dat(1,:,:,:))),2),1)));
    dat = dat(:,:,:,notnan);
    % clear old dat
    data.(cfg.datparam) = [];
  elseif strcmp(data.dimord,'chan_freq_epoch_tap')
    % Handle output from custom code with dimensions 'chan_freq_epoch_tap'
    % disp progress
    disp('input for SPACE is non-FieldTrip')
    if ndimsdat~=4 || ~datcomplex
      error('input for SPACE needs to be complex-valued and have 4 dimensions organized as ''chan_freq_epoch_tap''')
    end
    if ~isnumeric(data.(cfg.datparam))
      warning('not attempting memory and computation time optimization because data is specified as filename')
    else
      % only apply trick if ntaper exceeds nchan
      if size(data.(cfg.datparam),4)>size(data.(cfg.datparam),1)
        % disp progress
        disp('number of tapers exceeds number of channels, applying Eigendecomposition to reduce memory and computation time')
        nepoch = size(data.(cfg.datparam),3);
        nfreq  = size(data.(cfg.datparam),2);
        nchan  = size(data.(cfg.datparam),1);
        ntaper = nchan;
        dat    = complex(NaN(nchan,nfreq,nepoch,ntaper,datprecision),NaN(nchan,nfreq,nepoch,ntaper,datprecision));
        % construct new dat
        for iepoch = 1:nepoch
          for ifreq = 1:nfreq
            % first, select data and create csd
            currfour = permute(data.(cfg.datparam)(:,ifreq,iepoch,:),[1 4 2 3]);
            % get rid of NaNs (should always be the same over channels)
            currfour(:,isnan(currfour(1,:))) = [];
            % ensure double precision for calculations below
            currfour = double(currfour);
            % compute the csd
            csd = currfour*currfour';
            %%%%%%%%%
            % Reduce SPACE memory load and computation time by replacing each chan_taper matrix by the
            % Eigenvectors of its chan_chan cross-products weighted by sqrt(Eigenvalues).
            % This is possible because (1) SPACE only uses the cross-products of the chan_taper matrices
            % (i.e. the frequency- and trial-specific CSD) and (2) the Eigendecomposition of a symmetric
            % matrix A is A = VLV'.
            % As such, VL^.5 has the same cross-products as the original chan_taper matrix.
            [V L] = eig(csd);
            L     = diag(L);
            tol   = max(size(csd))*eps(max(L)); % compute tol using matlabs default
            zeroL = L<tol;
            eigweigth = V(:,~zeroL)*diag(sqrt(L(~zeroL)));
            % positive semi-definite failsafe
            if any(L<-tol) || any(~isreal(L))
              error('csd not positive semidefinite')
            end
            %%%%%%%%%
            % save in dat
            currm = size(eigweigth,2);
            dat(:,ifreq,iepoch,1:currm) = eigweigth;
          end
        end
        % trim dat (actual ntaper can be lower than ntaper, which depends on the data, hence the need for trimming)
        notnan = logical(squeeze(sum(sum(~isnan(squeeze(dat(1,:,:,:))),2),1)));
        dat = dat(:,:,:,notnan);
        % clear old dat
        data.(cfg.datparam) = [];
      end
    end
  else
    error('Input data structure is not supported for SPACE-FSP or SPACE-time. Please see function help for supported input.')
  end
end


% Set several easy to work with variables
if ~exist('dat','var')
  dat         = data.(cfg.datparam);
end
model          = cfg.model;
niter          = cfg.numiter;
convcrit       = cfg.convcrit;
degencrit      = cfg.degencrit;
ncomp          = cfg.ncomp;
nrand          = cfg.randstart;
nrandestcomp   = cfg.ncompestrandstart;
estnum         = [cfg.ncompeststart cfg.ncompestend cfg.ncompeststep];
estsrcritval   = cfg.ncompestsrcritval;
estsrcritjudge = cfg.ncompestsrcritjudge;
distcomp       = cfg.distcomp;
expvarinc      = cfg.ncompestvarinc;
corconval      = cfg.ncompestcorconval;
% set model specific ones and create modelopt
switch cfg.model
  case 'parafac'
    modelopt = {'compmodes',cfg.complexdims};
  case 'parafac2'
    modelopt = {'compmodes',cfg.complexdims,'specmodes',cfg.specialdims};
  case 'parafac2cp'
    if strcmp(cfg.ncompest,'splitrel')
      % do hacks for parafac2cp
      modelopt = {'compmodes',cfg.complexdims,'ssqdatnoncp',data.ssqdatnoncp,'specmodes',cfg.specialdims,'ssqdatnoncppart1',data.ssqdatnoncppart1,'ssqdatnoncppart2',data.ssqdatnoncppart2};
    else
      modelopt = {'compmodes',cfg.complexdims,'ssqdatnoncp',data.ssqdatnoncp,'specmodes',cfg.specialdims};
    end
  case 'spacetime'
    modelopt = {'freq',data.freq,'Dmode',cfg.Dmode};
  case 'spacefsp'
    modelopt = {'Dmode',cfg.Dmode};
  otherwise
    error('model not supported')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Component estimation
switch cfg.ncompest
  
  case 'splitrel'
    
    if any(strcmp(cfg.ncompestsrdatparam,{'oddeven','oddevenavg'}))
      if ~strcmp(model,'spacefsp') && ~strcmp(model,'spacetime')
        error('cfg.ncompestsrdatparam = ''oddeven/oddevenavg'' is only supported for SPACE-time and SPACE-FSP. Please provide partioned data in a 1x2 cell-array and specify its field name in cfg.ncompestsrdatparam')
      end
      if size(dat,3)==1
        error('splitrel procedure with automatic segmentation is only suitable for when the provided number trials/epochs is bigger than 1')
      end
      %  disp progress
      disp('creating split-half datasets using odd/even trial numbers')
      % extract partitions
      datsplit = cell(1,2);
      datsplit{1} = dat(:,:,1:2:size(dat,3),:);
      datsplit{2} = dat(:,:,2:2:size(dat,3),:);
    else
      % extract partitions
      nsplit = numel(data.(cfg.ncompestsrdatparam));
      datsplit = cell(1,nsplit);
      for isplit = 1:nsplit
        datsplit{isplit} = data.(cfg.ncompestsrdatparam){isplit};
      end
    end
    
    % perform splitrel component number estimate
    [ncomp, splitrelstat] = splitrel(model, dat, datsplit, nrandestcomp, estnum, estsrcritval, estsrcritjudge, niter, convcrit, degencrit, distcomp, oldsplithalf, modelopt{:}, 'checkpointpath', cfg.checkpointpath); % subfunction
    
    % extract startval if nrand is the same
    if ~oldsplithalf && (nrandestcomp==nrand)
      if ~isempty(splitrelstat.randomstatfullsucc)
        startval   = splitrelstat.randomstatfullsucc.startvalglobmin;
        randomstat = splitrelstat.randomstatfullsucc;
      else
        startval   = splitrelstat.randomstatfullfail.startvalglobmin;
        randomstat = splitrelstat.randomstatfullfail;
      end
    end
    
  case 'degeneracy'
    % Warn about component estimation
    warning('this is a very liberal method for component estimation')
    
    % estimate ncomp
    [ncomp] = degeneracy(model, dat, nrandestcomp, estnum, niter, convcrit, degencrit, distcomp, modelopt{:}, 'checkpointpath', cfg.checkpointpath); % subfunction
    
  case 'minexpvarinc'
    % estimate ncomp
    [ncomp] = minexpvarinc(model, dat, nrandestcomp, estnum, niter, convcrit, degencrit, distcomp, expvarinc, modelopt{:}, 'checkpointpath', cfg.checkpointpath); % subfunction
    
  case 'corcondiag'
    % estimate ncomp
    [ncomp, corcondiagstat] = corcondiag(model, dat, nrandestcomp, estnum, niter, convcrit, degencrit, distcomp, corconval, modelopt{:}, 'checkpointpath', cfg.checkpointpath); % subfunction
    
    % extract startval if nrand is the same
    if (nrandestcomp==nrand)
      if ~isempty(corcondiagstat.randomstatsucc)
        startval   = corcondiagstat.randomstatsucc.startvalglobmin;
        randomstat = corcondiagstat.randomstatsucc;
      else
        startval   = corcondiagstat.randomstatfail.startvalglobmin;
        randomstat = corcondiagstat.randomstatfail;
      end
    end
    
  case 'no'
    % do nothing, ncomp is set
    
  otherwise
    error('method for estimating number of components not supported')
    
end % switch ncompest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final decomposition
disp('performing FINAL decomposition')
% call low-level function and get start values if requested
% set up general ooptions
opt = {'niter', niter, 'convcrit', convcrit};
switch cfg.model
  case 'parafac'
    if isnumeric(nrand)
      if ~exist('startval','var') && ~exist('randomstat','var')
        [startval, randomstat] = randomstart(model, dat, ncomp, nrand, niter, convcrit, degencrit, distcomp, [], modelopt{:}); % subfunction
      end
      if strcmp(cfg.t3core,'yes')
        [comp, dum, expvar, scaling, tuckcongr, t3core] = feval(['nwaydecomp_' model], dat, ncomp, 'startval', startval, opt{:}, modelopt{:});
      else
        [comp, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, 'startval', startval, opt{:}, modelopt{:});
      end
    else
      if strcmp(cfg.t3core,'yes')
        [comp, dum, expvar, scaling, tuckcongr, t3core] = feval(['nwaydecomp_' model], dat, ncomp, opt{:}, modelopt{:});
      else
        [comp, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, opt{:}, modelopt{:});
      end
    end
  case 'parafac2'
    if isnumeric(nrand)
      if ~exist('startval','var') && ~exist('randomstat','var')
        [startval, randomstat] = randomstart(model, dat, ncomp, nrand, niter, convcrit, degencrit, distcomp, [], modelopt{:}); % subfunction
      end
      [comp, P, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, cfg.specialdims, 'startval', startval, opt{:}, modelopt{:});
    else
      [comp, P, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, cfg.specialdims, opt{:}, modelopt{:});
    end
  case 'parafac2cp'
    if isnumeric(nrand)
      if ~exist('startval','var') && ~exist('randomstat','var')
        [startval, randomstat] = randomstart(model, dat, ncomp, nrand, niter, convcrit, degencrit, distcomp, [], modelopt{:}); % subfunction
      end
      [comp, P, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, cfg.specialdims, 'startval', startval, opt{:}, modelopt{1:4});
    else
      [comp, P, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, cfg.specialdims, 'niter', niter, opt{:}, modelopt{1:4});
    end
  case 'spacetime'
    if isnumeric(nrand)
      if ~exist('startval','var') && ~exist('randomstat','var')
        [startval, randomstat] = randomstart(model, dat, ncomp, nrand, niter, convcrit, degencrit, distcomp, [], modelopt{:}); % subfunction
      end
      if strcmp(cfg.t3core,'yes')
        [comp, dum, dum, expvar, scaling, tuckcongr, t3core] = feval(['nwaydecomp_' model], dat, ncomp, data.freq, 'Dmode', cfg.Dmode, 'startval', startval, opt{:});
      else
        [comp, dum, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, data.freq, 'Dmode', cfg.Dmode, 'startval', startval, opt{:});
      end
    else
      if strcmp(cfg.t3core,'yes')
        [comp, dum, dum, expvar, scaling, tuckcongr, t3core] = feval(['nwaydecomp_' model], dat, ncomp, data.freq, 'Dmode', cfg.Dmode, opt{:});
      else
        [comp, dum, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, data.freq, 'Dmode', cfg.Dmode, opt{:});
      end
    end
  case 'spacefsp'
    if isnumeric(nrand)
      if ~exist('startval','var') && ~exist('randomstat','var')
        [startval, randomstat] = randomstart(model, dat, ncomp, nrand, niter, convcrit, degencrit, distcomp, [], modelopt{:}); % subfunction
      end
      if strcmp(cfg.t3core,'yes')
        [comp, dum, dum, expvar, scaling, tuckcongr, t3core] = feval(['nwaydecomp_' model], dat, ncomp, 'Dmode', cfg.Dmode, 'startval', startval, opt{:});
      else
        [comp, dum, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, 'Dmode', cfg.Dmode, 'startval', startval, opt{:});
      end
    else
      if strcmp(cfg.t3core,'yes')
        [comp, dum, dum, expvar, scaling, tuckcongr, t3core] = feval(['nwaydecomp_' model], dat, ncomp, 'Dmode', cfg.Dmode, opt{:});
      else
        [comp, dum, dum, expvar, scaling, tuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, 'Dmode', cfg.Dmode, opt{:});
      end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transform comp to compoment-specific cell-arrays
outputcomp = cell(1,ncomp);
switch model
  case {'parafac','parafac2','parafac2cp'}
    % comp always consist of a loading matrix per dimension
    for icomp = 1:ncomp
      for idim = 1:numel(comp)
        outputcomp{icomp}{idim} = comp{idim}(:,icomp);
      end
    end
  case 'spacetime'
    for icomp = 1:ncomp
      % A,B,C,S
      for iparam = 1:4
        outputcomp{icomp}{iparam} = comp{iparam}(:,icomp);
      end
      % D
      if strcmp(cfg.Dmode,'kdepcomplex')
        outputcomp{icomp}{5} = comp{5}(:,:,icomp);
      else % identity
        outputcomp{icomp}{5} = comp{5}(:,icomp);
      end
    end
  case 'spacefsp'
    for icomp = 1:ncomp
      % A,B,C
      for iparam = 1:3
        outputcomp{icomp}{iparam} = comp{iparam}(:,icomp);
      end
      % L
      outputcomp{icomp}{4} = comp{4}(:,:,icomp);
      % D
      if strcmp(cfg.Dmode,'kdepcomplex')
        outputcomp{icomp}{5} = comp{5}(:,:,icomp);
      else % identity
        outputcomp{icomp}{5} = comp{5}(:,icomp);
      end
    end
end

% set dimord
switch model
  case {'parafac','parafac2','parafac2cp'}
    dimord = data.dimord;
  case 'spacetime'
    dimord = 'A_B_C_S_D';
  case 'spacefsp'
    dimord = 'A_B_C_L_D';
end

% Construct output structure
nwaycomp.label      = data.label;
nwaycomp.dimord     = dimord;
nwaycomp.comp       = outputcomp;
if any(strcmp(model,{'parafac2','parafac2cp'}))
  nwaycomp.P        = P;
end
nwaycomp.expvar     = expvar;
nwaycomp.tuckcongr  = tuckcongr;
nwaycomp.scaling    = scaling;
if isnumeric(nrand)
  nwaycomp.randomstat = randomstat;
end
if strcmp(cfg.ncompest,'splitrel')
  if ~oldsplithalf
    nwaycomp.splitrelstat = splitrelstat;
  else
    nwaycomp.splithalfstat = splitrelstat;
  end
end
if strcmp(cfg.ncompest,'corcondiag')
  nwaycomp.corcondiagstat = corcondiagstat;
end
if exist('t3core','var')
  nwaycomp.t3core = t3core;
end


% add certain fields if present in input (some might be mutually exclusive, or undesired, adding all currenlty for completeness)
% general
fieldnames = {'grad','elec','trialinfo'};
nwaycomp = copyfields(data, nwaycomp, fieldnames);
% ft_freqanalysis/connectivityanalysis/timelockanalysis
fieldnames = {'freq','time','dof','labelcmb'}; % this explicitly does not contain 'cumtapcnt','cumsumcnt', as these are controlled for
nwaycomp = copyfields(data, nwaycomp, fieldnames);
% ft_componentanalysis
fieldnames = {'topo','topolabel','unmixing'};
nwaycomp = copyfields(data, nwaycomp, fieldnames);
% ft_sourceanalysis
fieldnames = {'pos','inside','outside','leadfield','dim','tri','transform'};
nwaycomp = copyfields(data, nwaycomp, fieldnames);


% includes certain fields for backwards compatability
fieldnames = {'ampfreq','phasfreq','freq_old','label_old'};
nwaycomp = copyfields(data, nwaycomp, fieldnames);


% do the general cleanup and bookkeeping at the end of the function
%ft_postamble debug % ft_preamble_debug currently leads to qsubfeval saving a copy of the input data
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous data
ft_postamble history nwaycomp
ft_postamble savevar nwaycomp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Subfunction for cfg.ncompest = 'degeneracy'            %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ncomp] = degeneracy(model, dat, nrand, estnum, niter, convcrit, degencrit, distcomp, varargin)

% get model specific options from keyval
% not necessary, not explicitly used
checkpointpath = keyval('checkpointpath', varargin);
if ~isempty(checkpointpath)
  docheckpoint = true;
  % create random identifier
  checkpointid = num2str(sum(clock.*1e6));
else
  docheckpoint = false;
end

% Display start
disp('degen-only: performing component number estimation by only checking degeneracy at each solution')


% Estimate number of components by incremently increasing number and calculating split-half correlations
ncomp  = 0;
for incomp = 1:estnum(2)
  disp(['degen-only: number of components = ' num2str(incomp) ' of max ' num2str(estnum(2))]);
  disp('degen-only: performing decomposition');
  
  
  % Get decompositions for current incomp
  % get start values for decomposition for current incomp
  [dum, randomstat] = randomstart(model, dat, incomp, nrand, niter, convcrit, degencrit, distcomp, ['degen-only ncomp = ' num2str(incomp) ' - '], varargin{:}); % subfunction
  
  % see if there are any non-degenerate start values and set flag if otherwise
  if length(randomstat.degeninit)==nrand
    % try again with another round
    [dum, randomstat] = randomstart(model, dat, incomp, nrand, niter, convcrit, degencrit, distcomp, ['degen-only ncomp = ' num2str(incomp) ' - '], varargin{:}); % subfunction
    if length(randomstat.degeninit)==nrand
      degenflg = true;
    else
      degenflg = false;
    end
  else
    degenflg = false;
  end
  
  
  % Act on degeneracy of solutions
  % first do a check for degeneracy (as split half makes no sense if solutions are degenerate)
  if degenflg
    disp('degen-only: random initializations only returned likely degenerate solutions')
    disp(['degen-only: final number of components = ' num2str(ncomp)]);
    break
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % save entire workspace except for the data in checkpoint if requested
  if docheckpoint
    cpfn    = [checkpointpath '_' checkpointid '_' 'incomp' num2str(incomp) '.mat'];
    varlist = who;
    varlist(strncmp(varlist,'dat',3)) = [];
    save(cpfn,varlist{:},'-v7.3')
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % If all checks are passed, update ncomp
  ncomp = incomp;
end % incomp



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Subfunction for cfg.ncompest = 'minexpvarinc'          %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ncomp] = minexpvarinc(model, dat, nrand, estnum, niter, convcrit, degencrit, distcomp, expvarinc, varargin)

% get model specific options from keyval
% not necessary, not explicitly used
checkpointpath = keyval('checkpointpath', varargin);
if ~isempty(checkpointpath)
  docheckpoint = true;
  % create random identifier
  checkpointid = num2str(sum(clock.*1e6));
else
  docheckpoint = false;
end

% Display start
disp('minexpvarinc: performing minimum increase in expvar component number estimation')


% Estimate number of components by incremently increasing number and calculating split-half correlations
% save current expvar
lastexpvar = 0;
ncomp  = 0;
for incomp = 1:estnum(2)
  disp(['minexpvarinc: number of components = ' num2str(incomp) ' of max ' num2str(estnum(2))]);
  disp('minexpvarinc: performing decomposition');
  
  
  % Get decompositions for current incomp
  % get start values for decomposition for current incomp
  [dum, randomstat] = randomstart(model, dat, incomp, nrand, niter, convcrit, degencrit, distcomp, ['minexpvarinc ncomp = ' num2str(incomp) ' - '], varargin{:}); % subfunction
  currexpvar = randomstat.expvar(1);
  
  % see if there are any non-degenerate start values and set flag if otherwise
  if length(randomstat.degeninit)==nrand
    degenflg = true;
  else
    degenflg = false;
  end
  
  
  % Act on degeneracy of solutions
  % first do a check for degeneracy (as split half makes no sense if solutions are degenerate)
  if degenflg
    disp('minexpvarinc: random initializations only returned likely degenerate solutions')
    disp(['minexpvarinc: final number of components = ' num2str(ncomp)]);
    break
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % save entire workspace except for the data in checkpoint if requested
  if docheckpoint
    cpfn    = [checkpointpath '_' checkpointid '_' 'incomp' num2str(incomp) '.mat'];
    varlist = who;
    varlist(strncmp(varlist,'dat',3)) = [];
    save(cpfn,varlist{:},'-v7.3')
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Act on explained variance increase
  if incomp ~= 1
    disp(['minexpvarinc: going from ncomp = ' num2str(ncomp) ' to ' num2str(incomp) '  lead to ' num2str((currexpvar - lastexpvar),'%-2.2f') '% increase in explained variance'])
    if (currexpvar - lastexpvar) < expvarinc
      disp('minexpvarinc: increase in explained variance not sufficient')
      disp(['minexpvarinc: final number of components = ' num2str(ncomp)]);
      break
    end
  end
  
  
  % If all checks are passed, update ncomp and lastexpvar
  ncomp = incomp;
  lastexpvar = currexpvar;
end % incomp



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Subfunction for cfg.ncompest = 'corcondiag'           %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ncomp, corcondiagstat] = corcondiag(model, dat, nrand, estnum, niter, convcrit, degencrit, distcomp, corconval, varargin)

% get model specific options from keyval
switch model
  case 'parafac'
    compmodes = keyval('compmodes', varargin);
  case 'parafac2'
    compmodes = keyval('compmodes', varargin);
    specmodes = keyval('specmodes', varargin);
  case 'parafac2cp'
    compmodes   = keyval('compmodes',   varargin);
    specmodes   = keyval('specmodes',   varargin);
    ssqdatnoncp = keyval('ssqdatnoncp', varargin);
  case 'spacetime'
    freq  = keyval('freq', varargin);
    Dmode = keyval('Dmode', varargin);
  case 'spacefsp'
    Dmode = keyval('Dmode', varargin);
  otherwise
    error('model not supported')
end
checkpointpath = keyval('checkpointpath', varargin);
if ~isempty(checkpointpath)
  docheckpoint = true;
  % create random identifier
  checkpointid = num2str(sum(clock.*1e6));
else
  docheckpoint = false;
end

% Display start
disp('corcondiag: performing component number estimation based on core consistency diagnostic')

% Estimate number of components by comparing the t3core array the identity-array of the same size
allcomp       = [];
allt3core     = [];
allcorcondiag = [];
allrandomstat = [];
ncompsucc     = [];
incomp        = estnum(1);
ncompfound    = false;
while ~ncompfound % the logic used here is identical as in splitrel, they should be changed SIMULTANEOUSLY or both subfunctions should be merged with switches
  disp(['corcondiag: number of components = ' num2str(incomp) ' of max ' num2str(estnum(2))]);
  disp('corcondiag: performing decomposition');
  
  
  % Get decompositions for current incomp
  % get quickly computed start values for decomposition for current incomp
  [startval, randomstat] = randomstart(model, dat, incomp, nrand, niter, convcrit, degencrit, distcomp, ['corcondiag ncomp = ' num2str(incomp) ' - '], varargin{:}); % subfunction
  
  
  % see if there are any non-degenerate start values and set flag if otherwise
  if numel(randomstat.degeninit)==nrand
    % try again with another round
    [startval, randomstat] = randomstart(model, dat, incomp, nrand, niter, convcrit, degencrit, distcomp, ['corcondiag ncomp = ' num2str(incomp) ', second try due to degeneracy - '], varargin{:}); % subfunction
    if numel(randomstat.degeninit)==nrand
      degenflg = true;
    else
      degenflg = false; % good to go
    end
  else
    degenflg = false;
  end
  
  
  % get final decompositions for current incomp
  opt = {'niter', niter, 'convcrit', convcrit, 'dispprefix',['corcondiag ncomp = ' num2str(incomp) ': ']};
  switch model
    case 'parafac'
      [estcomp,dum,dum,dum,dum,t3core] = feval(['nwaydecomp_' model], dat, incomp, 'startval', startval, 'compmodes', compmodes, opt{:});
      %     case 'parafac2' FIXME: parafac2/cp needs to output a t3core, and crosscompcongruence needs to be aware of specmodes (see splitrel)
      %       [estcomp,dum,dum,dum,dum,dum,t3core] = feval(['nwaydecomp_' model], dat, incomp, specmodes, 'startval', startval, 'compmodes', compmodes, opt{:});
      %     case 'parafac2cp'
      %       [estcomp,dum,dum,dum,dum,dum,t3core] = feval(['nwaydecomp_' model], dat, incomp, specmodes, 'startval', startval, 'compmodes', compmodes, opt{:}, 'ssqdatnoncp', ssqdatnoncp);
    case 'spacetime'
      [estcomp,dum,dum,dum,dum,dum,t3core] = feval(['nwaydecomp_' model], dat, incomp, freq, 'Dmode', Dmode, 'startval', startval, opt{:});
    case 'spacefsp'
      [estcomp,dum,dum,dum,dum,dum,t3core] = feval(['nwaydecomp_' model], dat, incomp, 'Dmode', Dmode, 'startval', startval, opt{:});
    otherwise
      error('model not yet supported in corcondiag component number estimation')
  end
  
  % create identity array and compute core consistency diagnostic
  if incomp==1
    corcondiag = 1;
  else
    ndimsdat = round((log(numel(t3core))/log(incomp))); % get ndimsdat out of the t3core, necessary when dat is a filename containing the data
    ident = zeros(incomp^ndimsdat,1);
    for iident = 1:incomp
      ind = iident;
      for idim = 2:ndimsdat
        ind = ind +incomp.^(idim-1) .* (iident-1);
      end
      ident(ind) = 1;
    end
    t3core = t3core(:); % vectorize
    corcondiag = 1 - (sum((ident-t3core).^2)/sum(t3core.^2));
  end
  
  % save incomp specifics
  allcomp{incomp}       = estcomp;
  allt3core{incomp}     = t3core;
  allcorcondiag{incomp} = corcondiag;
  allrandomstat{incomp} = randomstat;
  
  
  % set critfailflg
  if corcondiag < corconval
    critfailflg = true;
  else
    critfailflg = false;
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % save entire workspace except for the data in checkpoint if requested
  if docheckpoint
    cpfn    = [checkpointpath '_' checkpointid '_' 'incomp' num2str(incomp) '.mat'];
    varlist = who;
    varlist(strncmp(varlist,'dat',3)) = [];
    save(cpfn,varlist{:},'-v7.3')
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% determine incomp succes
  % Act on critfailflg and degenflg
  disp(['corcondiag: core consistency diagnostic at ncomp = ' num2str(incomp) ' was ' num2str(corcondiag,'%-2.2f')])
  if critfailflg || degenflg
    %%%%% FAIL
    
    % When current incomp fails, decrease incomp.
    % If incomp = 1, stop procedure
    % If not, then either a previous incomp has been performed or not
    % If so, they should have been a success. Check this.
    % If not, set incomp to round(incomp/2)
    % If yes, then it was either incomp-1 or not
    % If incomp-1 was performed, set incomp to incomp-1 and stop procedure
    % If incomp-1 has not been performed, find closest success. Set new incomp to the point haldway
    
    % update fail fields
    if ~degenflg
      disp(['corcondiag: minimum core consistency diagnostic of ' num2str(corconval) ' has not been reached at ncomp = ' num2str(incomp)]);
      t3corefail     = t3core;
      corcondiagfail = corcondiag;
      randomstatfail = randomstat;
      stopreason = 'core consistency diagnostic';
    elseif degenflg
      disp('corcondiag: random initializations only returned likely degenerate solutions')
      t3corefail     = []; % no t3core is computed
      corcondiagfail = []; % no t3core, so no corcondiag
      randomstatfail = randomstat;
      stopreason     = 'degeneracy';
    end
    ncompsucc{incomp} = false;
    
    % check for incomp == 1
    if incomp == 1
      warning('corcondiag: NO COMPONENTS REACHED CORE CONSISTENCY DIAGNOSTIC CRITERION');
      % set 'succes' status and final number of components
      ncompfound     = true;
      ncomp          = 1;
      t3coresucc     = [];
      corcondiagsucc = [];
      randomstatsucc = [];
      stopreason = 'core consistency diagnostic criterion fail at ncomp = 1';
      
    else % find a new incomp
      
      % first, sanity check
      if any(~[ncompsucc{1:incomp-1}])
        % this should not be possible
        error('unexpected error in ncomp estimation')
      end
      
      % find closest previous incomp
      prevsucc = find(~cellfun(@isempty,ncompsucc(1:incomp-1)),1,'last');
      if isempty(prevsucc)
        % no previous successes found, set new incomp to halfway to 0
        incomp = round(incomp/2);
      else
        % previous success found, determine new incomp based on prevsucc
        % if prevsucc is incomp-1, stop procedure
        if prevsucc == (incomp-1)
          ncomp      = incomp-1;
          ncompfound = true;
        else
          % if not, set new incomp to point in between prevsucc and incomp
          incomp = prevsucc + round((incomp-prevsucc)/2); % the increase is always at least one, as the diff between prevsucc and incomp is always >1
        end
      end
    end
  else
    %%%%% SUCCESS
    
    % When current incomp succeeds, increase incomp.
    % Then, incomp+i can either be at the maximum or not.
    % If at the maximum, succes = true.
    % If not at the maximum, increment.
    % Then, there have either been attempts at incomp+i or there have not.
    % If not, increment incomp with stepsize limited by the maximum.
    % If there have been incomp+i's, they can only have failed, and none should be present at incomp+stepsize+i
    % Check for this.
    % Then, find the closest fail. If it is incomp+1, succes = true. If not, increment with half the distances to the closest fail.
    
    % update succes fields and ncompsucc
    t3coresucc     = t3core;
    corcondiagsucc = corcondiag;
    randomstatsucc = randomstat;
    ncompsucc{incomp} = true;
    
    % check whether maximum has been reached, and increment incomp otherwise
    if incomp==estnum(2)
      disp(['corcondiag: succesfully reached a priori determined maximum of ' num2str(estnum(2)) ' components']);
      t3corefail     = [];
      corcondiagfail = [];
      randomstatfail = [];
      stopreason = ['corcondiag stopped: reached a priori maximum of ' num2str(estnum(2)) ' components'];
      % set succes status
      ncomp      = incomp;
      ncompfound = true;
      
    else % keep on incrementing
      % check for solutions at incomp+
      if isempty([ncompsucc{incomp+1:end}])
        % no incomp+ solutions found, increment incomp with step (and check for ncompestend)
        incomp = min(incomp + estnum(3),estnum(2));
        
      else
        % incomp+ fails detected
        % first, sanity check
        if any([ncompsucc{incomp+1:end}])
          % incomp+i should never be successful
          error('unexpected error in ncomp estimation')
        end
        
        % find closest fail
        nextfail = incomp + find(~cellfun(@isempty,ncompsucc(incomp+1:end)),1);
        
        % incomp+i found, check whether next ncomp was failing
        if nextfail == (incomp+1)
          % next ncomp failed, set succes status for current
          ncomp      = incomp;
          ncompfound = true;
        else
          % next incomp was not the fail, pick the point halfway to the next fail
          incomp = incomp + round((nextfail-incomp)./2);
        end
      end
    end
    
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
end % incomp
disp(['corcondiag: final number of components = ' num2str(ncomp)]);

% create corcondiagstat
corcondiagstat.ncomp           = ncomp;
corcondiagstat.corcondiagsucc  = corcondiagsucc;
corcondiagstat.corcondiagfail  = corcondiagfail;
corcondiagstat.t3coresucc      = t3coresucc;
corcondiagstat.t3corefail      = t3corefail;
corcondiagstat.randomstatsucc  = randomstatsucc;
corcondiagstat.randomstatfail  = randomstatfail;
corcondiagstat.stopreason      = stopreason;
corcondiagstat.crosscompcongr  = [];
corcondiagstat.allcomp         = allcomp;
corcondiagstat.allt3core       = allt3core;
corcondiagstat.allcorcondiag   = allcorcondiag;
corcondiagstat.allrandomstat   = allrandomstat;
corcondiagstat.ncompsucc       = ncompsucc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Subfunction for cfg.ncompest = 'splitrel'             %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ncomp,splitrelstat] = splitrel(model, datfull, datsplit, nrand, estnum, estsrccritval, estsrcritjudge, niter, convcrit, degencrit, distcomp, oldsplithalf, varargin)

% get model specific options from keyval
switch model
  case 'parafac'
    compmodes = keyval('compmodes', varargin);
  case 'parafac2'
    compmodes = keyval('compmodes', varargin);
    specmodes = keyval('specmodes', varargin);
  case 'parafac2cp'
    compmodes   = keyval('compmodes',   varargin);
    specmodes   = keyval('specmodes',   varargin);
    ssqdatnoncppart1 = keyval('ssqdatnoncppart1', varargin);
    ssqdatnoncppart2 = keyval('ssqdatnoncppart2', varargin);
  case 'spacetime'
    freq  = keyval('freq', varargin);
    Dmode = keyval('Dmode', varargin);
  case 'spacefsp'
    Dmode = keyval('Dmode', varargin);
  otherwise
    error('model not supported')
end
checkpointpath = keyval('checkpointpath', varargin);
if ~isempty(checkpointpath)
  docheckpoint = true;
  % create random identifier
  checkpointid = num2str(sum(clock.*1e6));
else
  docheckpoint = false;
end

% get N
nsplit = numel(datsplit);

%%% backwards compatability per August 2016 for oldsplithalf
% reverse the old style flag
newsplitrel = ~oldsplithalf;
%%% backwards compatability per August 2016 for oldsplithalf

% Display start
disp(['split-reliability: performing split-reliability component number estimation using ' num2str(nsplit) ' splits'])

% Estimate number of components by incremently increasing number and calculating split-rel correlations
allcompsrc         = [];
allpartcombcompsrc = [];
allrandomstatsplit = [];
allrandomstatfull  = [];
incomp  = estnum(1);
ncompfound  = false;
while ~ncompfound % the logic used here is identical as in corcondiag, they should be changed SIMULTANEOUSLY or both subfunctions should be merged with switches
  
  disp(['split-reliability: number of components = ' num2str(incomp) ' of max ' num2str(estnum(2))]);
  disp('split-reliability: performing decomposition and computing split-half coefficients');
  
  
  % Get decompositions for current incomp
  % for the full data
  if newsplitrel
    [dum, randomstatfull] = randomstart(model, datfull, incomp, nrand, niter, convcrit, degencrit, distcomp, ['split-reliability - full data  ncomp = ' num2str(incomp) ' - '], varargin{:}); % subfunction
  else
    %%% backwards compatability per August 2016 for oldsplithalf
    dum            = [];
    randomstatfull = [];
    %%% backwards compatability per August 2016 for oldsplithalf
  end
  % for the splits
  randomstatsplit = cell(1,nsplit);
  for isplit = 1:nsplit
    [dum, randomstatsplit{isplit}] = randomstart(model, datsplit{isplit}, incomp, nrand, niter, convcrit, degencrit, distcomp, ['split-reliability - part ' num2str(isplit) '/' num2str(nsplit) '  ncomp = ' num2str(incomp) ' - '], varargin{:}); % subfunction
  end
  
  % see if there are any non-degenerate start values for the full data and set flag if otherwise (and try again if it only goes for one partition)
  if newsplitrel
    degenflg = numel(randomstatfull.degeninit)==nrand;
    if degenflg % try again
      [dum, randomstatfull] = randomstart(model, datfull, incomp, nrand, niter, convcrit, degencrit, distcomp, ['split-reliability - full data  ncomp = ' num2str(incomp) ', second try due to degeneracy - '], varargin{:}); % subfunction
      degenflg = numel(randomstatfull.degeninit)==nrand;
    end
  else
    degenflg = false;
  end
  
  % see if there are any non-degenerate start values for the splits and set flag if otherwise (and try again if it only goes for one partition)
  if ~degenflg % only continue if the full data had nondegenerate start values
    degencnt = 0;
    for isplit = 1:nsplit
      degencnt = degencnt + (numel(randomstatsplit{isplit}.degeninit)==nrand);
    end
    if degencnt == nsplit
      degenflg = true;
    else % if only a subset of splits had no nondegenerate starts, try those again
      degencnt = 0;
      for isplit = 1:nsplit
        if randomstatsplit{isplit}.degeninit == nrand
          [dum, randomstatsplit{isplit}] = randomstart(model, datsplit{isplit}, incomp, nrand, niter, convcrit, degencrit, distcomp, ['split-reliability - part ' num2str(isplit) '/' num2str(nsplit) '  ncomp = ' num2str(incomp) ', second try due to degeneracy - '], varargin{:}); % subfunction
        end
        degencnt = degencnt + (numel(randomstatsplit{isplit}.degeninit)==nrand);
      end
      % check degeneracy again
      if degencnt>0
        degenflg = true;
      else
        degenflg = false;
      end
    end
  end
  
  
  %%%%
  % The splitrel logic is now as follows: the splitrel criterion fails, if the splitrel coeffcient is not surpassed for any combination of random starts of the splits with the main
  %%%
  
  % extract all non-degenerate startvalues
  % for the full data
  if newsplitrel
    fullcomp        = randomstatfull.startvalall(setdiff(1:nrand,randomstatfull.degeninit));
    nndegenrandfull = numel(fullcomp);
  else
    %%% backwards compatability per August 2016 for oldsplithalf
    % full data now becomes first split
    fullcomp        = randomstatsplit{1}.startvalall(setdiff(1:nrand,randomstatsplit{1}.degeninit));
    nndegenrandfull = numel(fullcomp);
    %%% backwards compatability per August 2016 for oldsplithalf
  end
  % for splits
  nndegenrandsplit = NaN(1,nsplit);
  splitcomp        = cell(1,nsplit);
  for isplit = 1:nsplit
    splitcomp{isplit}        = randomstatsplit{isplit}.startvalall(setdiff(1:nrand,randomstatsplit{isplit}.degeninit));
    nndegenrandsplit(isplit) = numel(splitcomp{isplit});
  end
  
  % compute splitrel coeffcients for all possible pairs of random starts from the splits with the main
  partcombcompsrc = cell(1,nsplit);
  for isplit = 1:nsplit
    compsrc = NaN(incomp,numel(splitcomp{1}{1})); % NaN in case of degenflg and the below is not executed
    partcombcompsrc{isplit} = cell(nndegenrandfull,nndegenrandsplit(isplit));
    for irandfull = 1:nndegenrandfull
      for irandsplit = 1:nndegenrandsplit(isplit)
        % set current estcomp
        currcomp = cell(1,2);
        currcomp{1} = fullcomp{irandfull};
        currcomp{2} = splitcomp{isplit}{irandsplit};
        % compute component congruence for all possible pairs between splits
        compcongr = zeros(incomp,incomp,length(currcomp{1}));
        for iparam = 1:length(currcomp{1})
          % perform model specific stuff
          switch model
            case {'parafac','parafac2','parafac2cp'}
              paramc1 = currcomp{1}{iparam};
              paramc2 = currcomp{2}{iparam};
              % normalize
              paramc1 = bsxfun(@rdivide,paramc1,sqrt(sum(abs(paramc1).^2,1)));
              paramc2 = bsxfun(@rdivide,paramc2,sqrt(sum(abs(paramc2).^2,1)));
              % put in compsrc
              compcongr(:,:,iparam) = abs(paramc1' * paramc2);
            case {'spacetime','spacefsp'}
              switch iparam
                case {1,2}
                  paramc1 = currcomp{1}{iparam};
                  paramc2 = currcomp{2}{iparam};
                  % normalize
                  paramc1 = bsxfun(@rdivide,paramc1,sqrt(sum(abs(paramc1).^2,1)));
                  paramc2 = bsxfun(@rdivide,paramc2,sqrt(sum(abs(paramc2).^2,1)));
                  % put in compsrc
                  compcongr(:,:,iparam) = abs(paramc1' * paramc2);
                case 3
                  paramc1 = currcomp{1}{iparam};
                  paramc2 = currcomp{2}{iparam};
                  % normalize
                  paramc1 = bsxfun(@rdivide,paramc1,sqrt(sum(abs(paramc1).^2,1)));
                  paramc2 = bsxfun(@rdivide,paramc2,sqrt(sum(abs(paramc2).^2,1)));
                  % put in compsrc
                  if size(paramc1,1) == size(paramc2,1)
                    compcongr(:,:,iparam) = abs(paramc1' * paramc2);
                  else
                    compcongr(:,:,iparam) = 0; % congruence can't be computed, set to maximally incongruent (0)
                  end
                case 4
                  switch model
                    case 'spacetime'
                      % create frequency specific phases weighted by spatial maps and frequency profiles
                      A1 = currcomp{1}{1};
                      A2 = currcomp{2}{1};
                      B1 = currcomp{1}{2};
                      B2 = currcomp{2}{2};
                      S1 = currcomp{1}{4};
                      S2 = currcomp{2}{4};
                      % normalize
                      A1 = bsxfun(@rdivide,A1,sqrt(sum(abs(A1).^2,1)));
                      A2 = bsxfun(@rdivide,A2,sqrt(sum(abs(A2).^2,1)));
                      B1 = bsxfun(@rdivide,B1,sqrt(sum(abs(B1).^2,1)));
                      B2 = bsxfun(@rdivide,B2,sqrt(sum(abs(B2).^2,1)));
                      % construct spatial phase maps
                      Scomp1 = ipermute(exp(1i*2*pi*bsxfun(@times,permute(S1,[3 1 2]),freq')),[3 1 2]);
                      Scomp2 = ipermute(exp(1i*2*pi*bsxfun(@times,permute(S2,[3 1 2]),freq')),[3 1 2]);
                      % scale with A
                      Scomp1 = bsxfun(@times,Scomp1,A1);
                      Scomp2 = bsxfun(@times,Scomp2,A2);
                      % compute splitrelcoef over freqs,
                      srcoverfreq = zeros(incomp,incomp,size(B1,1));
                      for ifreq = 1:size(B1,1)
                        srcoverfreq(:,:,ifreq) = abs(Scomp1(:,:,ifreq)'*Scomp2(:,:,ifreq));
                      end
                      % weight with average B and combine over freq
                      Bweight    = bsxfun(@times,permute(B1',[1 3 2]),permute(B2',[3 1 2]));
                      Bweight    = bsxfun(@rdivide,Bweight,sum(Bweight,3));
                      srcsumfreq = sum(srcoverfreq .* Bweight,3);
                      % put in compsrc
                      compcongr(:,:,iparam) = srcsumfreq;
                    case 'spacefsp'
                      % create frequency specific phases weighted by spatial maps and frequency profiles
                      A1 = currcomp{1}{1};
                      A2 = currcomp{2}{1};
                      B1 = currcomp{1}{2};
                      B2 = currcomp{2}{2};
                      L1 = currcomp{1}{4};
                      L2 = currcomp{2}{4};
                      % normalize
                      A1 = bsxfun(@rdivide,A1,sqrt(sum(abs(A1).^2,1)));
                      A2 = bsxfun(@rdivide,A2,sqrt(sum(abs(A2).^2,1)));
                      B1 = bsxfun(@rdivide,B1,sqrt(sum(abs(B1).^2,1)));
                      B2 = bsxfun(@rdivide,B2,sqrt(sum(abs(B2).^2,1)));
                      % construct spatial phase maps
                      Lcomp1 = exp(1i*2*pi*permute(L1,[1 3 2]));
                      Lcomp2 = exp(1i*2*pi*permute(L2,[1 3 2]));
                      % scale with A
                      Lcomp1 = bsxfun(@times,Lcomp1,A1);
                      Lcomp2 = bsxfun(@times,Lcomp2,A2);
                      % compute splitrelcoef over freqs,
                      srcoverfreq = zeros(incomp,incomp,size(B1,1));
                      for ifreq = 1:size(B1,1)
                        srcoverfreq(:,:,ifreq) = abs(Lcomp1(:,:,ifreq)'*Lcomp2(:,:,ifreq));
                      end
                      % weight with average B and combine over freq
                      Bweight    = bsxfun(@times,permute(B1',[1 3 2]),permute(B2',[3 1 2]));
                      Bweight    = bsxfun(@rdivide,Bweight,sum(Bweight,3));
                      srcsumfreq = sum(srcoverfreq .* Bweight,3);
                      % put in compsrc
                      compcongr(:,:,iparam) = srcsumfreq;
                  end
                case 5
                  switch Dmode
                    case 'identity'
                      % D is fixed with arbitrary order, make its splitrel coefficient irrelevant
                      compcongr(:,:,iparam) = 1;
                    case 'kdepcomplex'
                      B1 = currcomp{1}{2};
                      B2 = currcomp{2}{2};
                      D1 = currcomp{1}{5};
                      D2 = currcomp{2}{5};
                      % weight with B
                      D1 = bsxfun(@times,D1,permute(B1,[1 3 2]));
                      D2 = bsxfun(@times,D2,permute(B2,[1 3 2]));
                      % vectorize
                      paramc1 = reshape(permute(D1,[3 1 2]),[incomp numel(B1)]).';
                      paramc2 = reshape(permute(D2,[3 1 2]),[incomp numel(B2)]).';
                      % normalize
                      paramc1 = bsxfun(@rdivide,paramc1,sqrt(sum(abs(paramc1).^2,1)));
                      paramc2 = bsxfun(@rdivide,paramc2,sqrt(sum(abs(paramc2).^2,1)));
                      % put in compsrc
                      compcongr(:,:,iparam) = abs(paramc1' * paramc2);
                  end
              end
          end
        end
        % get splitrel coefficients by selecting most-similair unique pairings, but only for those that matter for the splitrel coeff
        % 'clean' compcongr of estsrcritval==0
        compcongr(:,:,estsrccritval==0) = NaN;
        compsrc  = zeros(incomp,length(currcomp{1}));
        congrsum = nansum(compcongr,3);
        % match from perspective of first split (main) (i.e. find components of split 2 that match those of split 1)
        % do so by starting from the component-pair with the highest similarity, then the next most similar, etc.
        set1ind   = zeros(1,incomp);
        set2ind   = zeros(1,incomp);
        for icomp = 1:incomp
          [dum, set1ind(icomp)] = max(max(congrsum,[],2));
          [dum, set2ind(icomp)] = max(congrsum(set1ind(icomp),:));
          congrsum(set1ind(icomp),:) = 0;
          congrsum(:,set2ind(icomp)) = 0;
        end
        % sanity check
        if any(diff(sort(set1ind))==0) || any(diff(sort(set2ind))==0)
          error('some components were selected multiple times')
        end
        % sort for convenience
        [set1ind, sortind] = sort(set1ind);
        set2ind = set2ind(sortind);
        for iparam = 1:length(currcomp{1})
          compsrc(:,iparam) = diag(compcongr(set1ind,set2ind,iparam));
        end
        % save compsrc
        partcombcompsrc{isplit}{irandfull,irandsplit} = compsrc;
      end % irandsplit
    end % irandfull
  end % isplit
  
  
  % determine failure or succes of current incomp
  if ~degenflg
    
    % deal with compsrcs
    compsrc = NaN([incomp numel(splitcomp{1}{1}) nsplit]);
    switch estsrcritjudge
      
      case {'meanoversplitscons','minofsplitscons'}
        % CONSERVATIVE pick best the first of each split
        for isplit = 1:nsplit
          compsrc(:,:,isplit) = partcombcompsrc{isplit}{1,1};
        end
      
      case {'meanoversplitslib','minofsplitslib'}
        % LIBERAL pick best possible compsrc to pass on per split, 'best' = randinitcomb with highest minimum of src - criterion
        for isplit = 1:nsplit
          currcompsrc = partcombcompsrc{isplit};
          for irandfull = 1:nndegenrandfull
            for irandsplit = 1:nndegenrandsplit(isplit)
              currcompsrc{irandfull,irandsplit} = currcompsrc{irandfull,irandsplit} - repmat(estsrccritval,[incomp 1]);
              currcompsrc{irandfull,irandsplit} = currcompsrc{irandfull,irandsplit}(:,estsrccritval~=0);
            end
          end
          maxminsrccoeff = cellfun(@min,cellfun(@min,currcompsrc,'uniformoutput',0));
          [dum maxind] = max(maxminsrccoeff(:));
          [rowind,colind] = ind2sub([nndegenrandfull nndegenrandsplit(isplit)],maxind);
          compsrc(:,:,isplit) = partcombcompsrc{isplit}{rowind,colind};
        end
        
      otherwise
        error('specified cfg.ncompestsrcritjudge not supported')
    end
    
    %%% backwards compatability per August 2016 for oldsplithalf
    if ~newsplitrel
      if ~isempty(compsrc)
        compsrc = compsrc(:,:,2);
      end
    end
    %%% backwards compatability per August 2016 for oldsplithalf
    
    % deal with multiple splits
    switch estsrcritjudge
      case {'meanoversplitscons','meanoversplitslib'}
        compsrc = mean(compsrc,3);
      case {'minofsplitscons','minofsplitslib'}
        % do nothing
      otherwise
        error('specified cfg.ncompestsrcritjudge not supported')
    end
    
    % check whether the best possible compsrc's failed
    if any(any(any(compsrc < repmat(estsrccritval,[incomp 1 size(compsrc,3)]))))
      critfailflg = true;
    else
      critfailflg = false;
    end
    
  else
    critfailflg = false; % most accurate given degenflg is false
  end
  
  % save incomp specifics
  allpartcombcompsrc{incomp}  = partcombcompsrc;
  allcompsrc{incomp}          = compsrc;
  allrandomstatfull{incomp}   = randomstatfull;
  allrandomstatsplit{incomp}  = randomstatsplit;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % save entire workspace except for the data in checkpoint if requested
  if docheckpoint
    cpfn    = [checkpointpath '_' checkpointid '_' 'incomp' num2str(incomp) '.mat'];
    varlist = who;
    varlist(strncmp(varlist,'dat',3)) = [];
    save(cpfn,varlist{:},'-v7.3')
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% determine incomp succes
  % Act on critfailflg and degenflg
  disp(['split-reliability: lowest relevant absolute split-reliability coefficient were ' num2str(min(min(compsrc(:,estsrccritval~=0,:),[],1),[],2),'% .2f')]);
  if critfailflg || degenflg
    %%%%% FAIL
    
    % When current incomp fails, decrease incomp.
    % If incomp = 1, stop procedure
    % If not, then either a previous incomp has been performed or not
    % If so, they should have been a success. Check this.
    % If not, set incomp to round(incomp/2)
    % If yes, then it was either incomp-1 or not
    % If incomp-1 was performed, set incomp to incomp-1 and stop procedure
    % If incomp-1 has not been performed, find closest success. Set new incomp to the point haldway
    
    % update fail fields
    if ~degenflg
      disp(['split-reliability: one or more components did not reach split-reliability criterion of ' num2str(estsrccritval)]);
      compsrcfail         = compsrc;
      randomstatfullfail  = randomstatfull;
      randomstatsplitfail = randomstatsplit;
      stopreason          = 'split-reliability criterion';
    elseif degenflg
      disp('split-reliability: random initializations only returned likely degenerate solutions')
      compsrcfail         = compsrc;
      randomstatfullfail  = randomstatfull;
      randomstatsplitfail = randomstatsplit;
      stopreason          = 'degeneracy';
    end
    ncompsucc{incomp} = false;
    
    % check for incomp == 1
    if incomp == 1
      warning('split-reliability: NO COMPONENTS REACHED SPLIT-RELIABILITY CRITERION');
      % set 'succes' status and final number of components
      ncompfound         = true;
      ncomp               = 1;
      compsrcsucc         = [];
      randomstatfullsucc  = [];
      randomstatsplitsucc = [];
      stopreason          = 'split-reliability criterion fail at ncomp = 1';
      
    else % find a new incomp
      
      % first, sanity check
      if any(~[ncompsucc{1:incomp-1}])
        % this should not be possible
        error('unexpected error in ncomp estimation')
      end
      
      % find closest previous incomp
      prevsucc = find(~cellfun(@isempty,ncompsucc(1:incomp-1)),1,'last');
      if isempty(prevsucc)
        % no previous successes found, set new incomp to halfway to 0
        incomp = round(incomp/2);
      else
        % previous success found, determine new incomp based on prevsucc
        % if prevsucc is incomp-1, stop procedure
        if prevsucc == (incomp-1)
          ncomp      = incomp-1;
          ncompfound = true;
        else
          % if not, set new incomp to point in between prevsucc and incomp
          incomp = prevsucc + round((incomp-prevsucc)/2); % the increase is always at least one, as the diff between prevsucc and incomp is always >1
        end
      end
    end
  else
    %%%%% SUCCESS
    
    % When current incomp succeeds, increase incomp.
    % Then, incomp+i can either be at the maximum or not.
    % If at the maximum, succes = true.
    % If not at the maximum, increment.
    % Then, there have either been attempts at incomp+i or there have not.
    % If not, increment incomp with stepsize limited by the maximum.
    % If there have been incomp+i's, they can only have failed, and none should be present at incomp+stepsize+i
    % Check for this.
    % Then, find the closest fail. If it is incomp+1, succes = true. If not, increment with half the distances to the closest fail.
    
    % update succes fields and ncompsucc
    compsrcsucc         = compsrc;
    randomstatfullsucc  = randomstatfull;
    randomstatsplitsucc = randomstatsplit;
    ncompsucc{incomp}   = true;
    
    % check whether maximum has been reached, and increment incomp otherwise
    if incomp==estnum(2)
      disp(['split-reliability: succesfully reached a priori determined maximum of ' num2str(estnum(2)) ' components']);
      compsrcfail         = [];
      randomstatfullfail  = [];
      randomstatsplitfail = [];
      stopreason          = ['reached a priori maximum of ' num2str(estnum(2)) ' components'];
      % set succes status
      ncomp      = incomp;
      ncompfound = true;
      
    else % keep on incrementing
      % check for solutions at incomp+
      if isempty([ncompsucc{incomp+1:end}])
        % no incomp+ solutions found, increment incomp with step (and check for ncompestend)
        incomp = min(incomp + estnum(3),estnum(2));
        
      else
        % incomp+ fails detected
        % first, sanity check
        if any([ncompsucc{incomp+1:end}])
          % incomp+i should never be successful
          error('unexpected error in ncomp estimation')
        end
        
        % find closest fail
        nextfail = incomp + find(~cellfun(@isempty,ncompsucc(incomp+1:end)),1);
        
        % incomp+i found, check whether next ncomp was failing
        if nextfail == (incomp+1)
          % next ncomp failed, set succes status for current
          ncomp      = incomp;
          ncompfound = true;
        else
          % next incomp was not the fail, pick the point halfway to the next fail
          incomp = incomp + round((nextfail-incomp)./2);
        end
      end
    end
    
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
end % incomp
disp(['split-reliability: final number of components = ' num2str(ncomp)]);


% create splitrelstat
if newsplitrel
  splitrelstat.ncomp               = ncomp;
  splitrelstat.splitrelcsucc       = compsrcsucc;
  splitrelstat.splitrelcfail       = compsrcfail;
  splitrelstat.randomstatfullsucc  = randomstatfullsucc;
  splitrelstat.randomstatfullfail  = randomstatfullfail;
  splitrelstat.randomstatsplitsucc = randomstatsplitsucc;
  splitrelstat.randomstatsplitfail = randomstatsplitfail;
  splitrelstat.stopreason          = stopreason;
  splitrelstat.allcompsrc          = allcompsrc;
  splitrelstat.allpartcombcompsrc  = allpartcombcompsrc;
  splitrelstat.allrandomstatfull   = allrandomstatfull;
  splitrelstat.allrandomstatsplit  = allrandomstatsplit;
else
  %%% backwards compatability per August 2016 for oldsplithalf
  splitrelstat.ncomp             = ncomp;
  splitrelstat.splithcsucc       = compsrcsucc;
  splitrelstat.splithcfail       = compsrcfail;
  if isempty(compsrcsucc)
    splitrelstat.randomstat1succ = [];
    splitrelstat.randomstat2succ = [];
  else
    splitrelstat.randomstat1succ = randomstatsplitsucc{1};
    splitrelstat.randomstat2succ = randomstatsplitsucc{2};
  end
  if isempty(compsrcfail)
    splitrelstat.randomstat1fail = [];
    splitrelstat.randomstat2fail = [];
  else
    splitrelstat.randomstat1fail = randomstatsplitfail{1};
    splitrelstat.randomstat2fail = randomstatsplitfail{2};
  end
  splitrelstat.failreason        = stopreason;
  splitrelstat.allcomp           = [];
  splitrelstat.allrandomstat     = allrandomstatsplit;
  splitrelstat.allcompsh         = allcompsrc;
  splitrelstat.allpartcombcompsh = allpartcombcompsrc;
  for icorr = 1:numel(allcompsrc)
    if ~isempty(allcompsrc{icorr})
      splitrelstat.allcompsh{icorr}         = splitrelstat.allcompsh{icorr}(:,:,2); % first is between split1 and split1, 2 is between split1 and split2
      splitrelstat.allpartcombcompsh{icorr} = allpartcombcompsrc{icorr}{2}; % first is between split1 and split1, 2 is between split1 and split2
    end
  end
  %%% backwards compatability per August 2016 for oldsplithalf
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Subfunction for calculating randomstart-values         %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [startval, randomstat] = randomstart(model, dat, ncomp, nrand, niter, convcrit, degencrit, distcomp, dispprefix, varargin)

% get model specific options from keyval
switch model
  case 'parafac'
    compmodes = keyval('compmodes', varargin);
  case 'parafac2'
    compmodes = keyval('compmodes', varargin);
    specmodes = keyval('specmodes', varargin);
  case 'parafac2cp'
    compmodes   = keyval('compmodes',   varargin);
    specmodes   = keyval('specmodes',   varargin);
    ssqdatnoncp = keyval('ssqdatnoncp', varargin);
  case 'spacetime'
    freq  = keyval('freq', varargin);
    Dmode = keyval('Dmode', varargin);
  case 'spacefsp'
    Dmode = keyval('Dmode', varargin);
  otherwise
    error('model not supported')
end


% start random decompositions
niter = round(niter); % just to make sure it's an integer
if ~isempty(distcomp.system)
  switch distcomp.system
    
    case {'torque','p2p'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% QSUB distribution of random initializations
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      disp([dispprefix 'random start: determining optimal starting values from random initialization']);
      disp([dispprefix 'random start: decomposition of ' num2str(nrand) ' random initializations will be distributed using ' distcomp.system]);
      
      % prepare cell arrays for input
      cellncomp         = repmat({ncomp},[nrand 1]);
      cellniterkey      = repmat({'niter'},[nrand 1]);
      cellniterval      = repmat({niter},[nrand 1]);
      cellconvcritkey   = repmat({'convcrit'},[nrand 1]);
      cellconvcritval   = repmat({convcrit},[nrand 1]);
      celldispprefixkey = repmat({'dispprefix'},[nrand 1]);
      celldispprefixval = repmat({[dispprefix 'random start: ']},[nrand 1]);
      % create a cell-array containing copys of the data, or save and pass temporary filename
      % implemented for all models (file is deleted below)
      if isnumeric(dat) && ~isempty(distcomp.inputsaveprefix) && ischar(distcomp.inputsaveprefix)
        % make a random file name
        rng(sum(clock.*1e6))
        randname = tempname;
        filename = [distcomp.inputsaveprefix 'nwaytemp_' randname(end-6:end) '.mat'];
        % put in cell-array input
        celldat  = repmat({filename},[nrand 1]);
        % save data
        save(filename,'dat','-v7.3')
      else
        celldat = repmat({dat},[nrand 1]);
      end
      
      % set up general options
      opt = {cellniterkey, cellniterval, cellconvcritkey, cellconvcritval, celldispprefixkey, celldispprefixval};
      if isempty(distcomp.memreq)
        if ~isnumeric(dat)
          s = whos('-file', dat);
        else
          s = whos('dat');
        end
        memreq = s.bytes * 4; % probably enough for the fourier algorithms, probably not for the parafac ones
      else
        memreq = distcomp.memreq;
      end
      
      % set up distributed computing specific options
      switch distcomp.system
        case 'p2p'
          distcompopt = {'uniformoutput','false','resubmitdelay',distcomp.p2presubdel,'timreq', distcomp.timreq};
          distcompfun = @peercellfun;
        case 'torque'
          % increment timreq (should really be done based on size of data)
          distcomp.timreq = distcomp.timreq * ceil(ncomp/10);
          distcompopt = {'backend','torque','queue',distcomp.torquequeue,'timreq', distcomp.timreq,'matlabcmd',distcomp.matlabcmd,'stack',distcomp.torquestack,'options',['-V ', distcomp.qsuboptions],'sleep',30};
          distcompfun = @qsubcellfun;
        otherwise
          error('distributed computing system not supported')
      end
      
      % state current time as que submission time
      disp(['submitting to queue on: ' datestr(now)])
      
      % distribute function calls to peers
      switch model
        case 'parafac'
          cellcompmodeskey  = repmat({'compmodes'},[nrand 1]);
          cellcompmodesval  = repmat({compmodes},[nrand 1]);
          [cellcomp,cellssqres,cellexpvar,cellscaling,celltuckcongr] = feval(distcompfun,['nwaydecomp_' model ],celldat,cellncomp,cellcompmodeskey,cellcompmodesval,opt{:},distcompopt{:},'memreq', memreq);
        case 'parafac2'
          cellcompmodeskey  = repmat({'compmodes'},[nrand 1]);
          cellcompmodesval  = repmat({compmodes},[nrand 1]);
          cellspecmodes     = repmat({specmodes},[nrand 1]);
          [cellcomp,dum,cellssqres,cellexpvar,cellscaling,celltuckcongr] = feval(distcompfun,['nwaydecomp_' model ],celldat,cellncomp,cellspecmodes,cellcompmodeskey,cellcompmodesval,opt{:},distcompopt{:},'memreq', memreq);
        case 'parafac2cp'
          cellcompmodeskey   = repmat({'compmodes'},[nrand 1]);
          cellcompmodesval   = repmat({compmodes},[nrand 1]);
          cellspecmodes      = repmat({specmodes},[nrand 1]);
          cellssqdatnoncpkey = repmat({'ssqdatnoncp'},[nrand 1]);
          cellssqdatnoncpval = repmat({ssqdatnoncp},[nrand 1]);
          [cellcomp,dum,cellssqres,cellexpvar,cellscaling,celltuckcongr] = feval(distcompfun,['nwaydecomp_' model ],celldat,cellncomp,cellspecmodes,cellcompmodeskey,cellcompmodesval,cellssqdatnoncpkey,cellssqdatnoncpval,opt{:},distcompopt{:},'memreq', memreq);
        case 'spacetime'
          cellfreq = repmat({freq},[nrand 1]);
          cellDmodekey  = repmat({'Dmode'},[nrand 1]);
          cellDmodeval  = repmat({Dmode},[nrand 1]);
          [cellcomp,dum,cellssqres,cellexpvar,cellscaling,celltuckcongr] = feval(distcompfun,['nwaydecomp_' model ],celldat,cellncomp,cellfreq,cellDmodekey,cellDmodeval,opt{:},distcompopt{:},'memreq', memreq);
        case 'spacefsp'
          cellDmodekey  = repmat({'Dmode'},[nrand 1]);
          cellDmodeval  = repmat({Dmode},[nrand 1]);
          [cellcomp,dum,cellssqres,cellexpvar,cellscaling,celltuckcongr] = feval(distcompfun,['nwaydecomp_' model ],celldat,cellncomp,cellDmodekey,cellDmodeval,opt{:},distcompopt{:},'memreq', memreq);
        otherwise
          error('model not yet supported in automatic random starting')
      end
      % delete temporary copy of data if it was created
      if isnumeric(dat) && ~isempty(distcomp.inputsaveprefix) && ischar(distcomp.inputsaveprefix)
        delete(filename)
      end
      
      % format output to ssqres with output below
      randssqres    = [cellssqres{:}];
      randcomp      = cellcomp(:).';
      randexpvar    = [cellexpvar{:}];
      randscaling   = cellscaling(:).';
      randtuckcongr = celltuckcongr(:).';
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
    case 'matlabpct'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% Matlab PARFOR distribution of random initializations
      % pool is opened and closed here to hopefully circumvent memory leak in matlab's PCT
      % first, check wether parfeval is available
      useparfeval = ~verLessThan('matlab','8.2'); % (<2013b) this is a safe gamble, it might be available earlier though
      
      %%% open pool
      % open pool, always create a new pool to not interfere with concurrent work
      haspoolsize = ~isempty(distcomp.mpctpoolsize);
      hascluster  = ~isempty(distcomp.mpctcluster);
      if haspoolsize && ~hascluster
        if ~verLessThan('matlab','8.2') % parpool is a 2013b function
          poolobj = parpool(distcomp.mpctpoolsize,'IdleTimeout',inf);
        else
          matlabpool(distcomp.mpctpoolsize);
        end
      elseif ~haspoolsize && hascluster
        if ~verLessThan('matlab','8.2') % parpool is a 2013b function
          poolobj = parpool(distcomp.mpctcluster,'IdleTimeout',inf);
        else
          matlabpool(distcomp.mpctcluster);
        end
      elseif haspoolsize && hascluster
        if ~verLessThan('matlab','8.2') % parpool is a 2013b function
          poolobj = parpool(distcomp.mpctcluster,distcomp.mpctpoolsize,'IdleTimeout',inf);
        else
          matlabpool(distcomp.mpctcluster,distcomp.mpctpoolsize);
        end
      else
        if ~verLessThan('matlab','8.2')  % parpool is a 2013b function
          poolobj = parpool('IdleTimeout',inf);
        else
          matlabpool;
        end
      end
      %%% open pool
      
      %%% The below code is identical to local computation specified below except for parfor, keep it as such
      % allocate
      randssqres    = zeros(1,nrand);
      randcomp      = cell(1,nrand);
      randexpvar    = zeros(1,nrand);
      randscaling   = cell(1,nrand);
      randtuckcongr = cell(1,nrand);
      disp([dispprefix 'random start: determining optimal starting values from random initialization']);
      disp([dispprefix 'random start: decomposition of ' num2str(nrand) ' random initializations will be distributed using MATLABs parallel computing toolbox']);
      % first, gather inputs for each mode (necessary for parfor)
      % set up general options
      opt = {'niter', niter, 'convcrit', convcrit};
      switch model
        case 'parafac'
          modelinput = {dat, ncomp, 'compmodes', compmodes, opt{:}};
        case 'parafac2'
          modelinput = {dat, ncomp, specmodes, 'compmodes', compmodes, opt{:}};
        case 'parafac2cp'
          modelinput = {dat, ncomp, specmodes, 'compmodes', compmodes, opt{:}, 'ssqdatnoncp', ssqdatnoncp};
        case 'spacetime'
          modelinput = {dat, ncomp, freq, 'Dmode', Dmode, opt{:}};
        case 'spacefsp'
          modelinput = {dat, ncomp, 'Dmode', Dmode, opt{:}};
        otherwise
          error('model not yet supported in automatic random starting')
      end
      % start distribution
      if ~useparfeval
        
        % PARFOR way of handling parallelization
        parfor irand = 1:nrand
          disp([dispprefix 'parallel random initialization ' num2str(irand) ' of ' num2str(nrand) ' started']);
          % set up general options
          opt = {'dispprefix', [dispprefix 'parallel random start ' num2str(irand) ': ']};
          switch model
            case 'parafac'
              [randcomp{irand},    randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = feval(['nwaydecomp_' model], modelinput{:}, opt{:});
            case 'parafac2'
              [randcomp{irand},dum,randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = feval(['nwaydecomp_' model], modelinput{:}, opt{:});
            case 'parafac2cp'
              [randcomp{irand},dum,randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = feval(['nwaydecomp_' model], modelinput{:}, opt{:});
            case 'spacetime'
              [randcomp{irand},dum,randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = feval(['nwaydecomp_' model], modelinput{:}, opt{:});
            case 'spacefsp'
              [randcomp{irand},dum,randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = feval(['nwaydecomp_' model], modelinput{:}, opt{:});
            otherwise
              error('model not yet supported in automatic random starting')
          end
          disp([dispprefix 'parallel random initialization ' num2str(irand) ' of ' num2str(nrand) ' finished:  exp. var. = ' num2str(randexpvar(irand),'%-2.2f') '%']);
        end
        
      else
        
        % PARFEVAL way of handling parallelization
        % explicitly clear jobid
        clear jobid
        % start timing
        stopwatch  = tic;
        timeused   = NaN(1,nrand);
        for irand = 1:nrand
          disp([dispprefix 'parallel random initialization ' num2str(irand) ' of ' num2str(nrand) ' submitted']);
          % set up general options
          opt = {'niter', niter, 'convcrit', convcrit, 'dispprefix',[dispprefix 'parallel random start ' num2str(irand) ': ']};
          switch model
            case 'parafac'
              jobid(irand) = parfeval(['nwaydecomp_' model], 5, dat, ncomp, 'compmodes', compmodes, opt{:});
            case 'parafac2'
              jobid(irand) = parfeval(['nwaydecomp_' model], 6, dat, ncomp, specmodes, 'compmodes', compmodes, opt{:});
            case 'parafac2cp'
              jobid(irand) = parfeval(['nwaydecomp_' model], 6, dat, ncomp, specmodes, 'compmodes', compmodes, opt{:}, 'ssqdatnoncp', ssqdatnoncp);
            case 'spacetime'
              jobid(irand) = parfeval(['nwaydecomp_' model], 6, dat, ncomp, freq, 'Dmode', Dmode, opt{:});
            case 'spacefsp'
              jobid(irand) = parfeval(['nwaydecomp_' model], 6, dat, ncomp, 'Dmode', Dmode, opt{:});
            otherwise
              error('model not yet supported in automatic random starting')
          end
        end
        % fetch outputs
        for irand = 1:nrand
          switch model
            case 'parafac'
              [currid, randcomp{irand},    randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = fetchNext(jobid);
            case 'parafac2'
              [currid, randcomp{irand},dum,randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = fetchNext(jobid);
            case 'parafac2cp'
              [currid, randcomp{irand},dum,randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = fetchNext(jobid);
            case 'spacetime'
              [currid, randcomp{irand},dum,randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = fetchNext(jobid);
            case 'spacefsp'
              [currid, randcomp{irand},dum,randssqres(irand),randexpvar(irand),randscaling{irand},randtuckcongr{irand}] = fetchNext(jobid);
            otherwise
              error('model not yet supported in automatic random starting')
          end
          % extract execution time from diary
          diaryout = jobid(currid).Diary;
          exectime = diaryout(strfind(diaryout,'execution took'):end);
          nind     = regexp(exectime,'[0-9]');
          if ~isempty(nind)
            exectime = exectime(nind(1):nind(end));
            exectime = str2double(exectime);
            timeused(currid) = exectime;
          end
          disp([dispprefix 'parallel random initialization ' num2str(currid) ' (' num2str(irand) '/' num2str(nrand) ') returned and took ' num2str(timeused(currid),'%.1f') ' sec | exp. var. = ' num2str(randexpvar(irand),'%-2.2f') '%']);
        end
        % display time used (from qsubcellfun)
        fprintf('computational time = %.1f sec, elapsed = %.1f sec, speedup %.1f x\n', nansum(timeused), toc(stopwatch), nansum(timeused)/toc(stopwatch));
      end
      
      %%% close pool
      % if pool is present, close it, no longer needed
      if ~verLessThan('matlab','8.2') % parpool is a 2013b function
        delete(poolobj);
      else
        matlabpool close
      end
      %%% close pool
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
      error('specified distributed computing system not supported')
  end
  
else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% LOCAL compution of random initializations
  % allocate
  randssqres    = zeros(1,nrand);
  randcomp      = cell(1,nrand);
  randexpvar    = zeros(1,nrand);
  randscaling   = cell(1,nrand);
  randtuckcongr = cell(1,nrand);
  disp([dispprefix 'random start: determining optimal starting values from random initialization']);
  for irand = 1:nrand
    disp([dispprefix 'random start ' num2str(irand) ': decomposition of random initialization ' num2str(irand) ' of ' num2str(nrand) ' started']);
    % set up general options
    opt = {'niter', niter, 'convcrit', convcrit, 'dispprefix',[dispprefix 'random start ' num2str(irand) ': ']};
    switch model
      case 'parafac'
        [rndcomp,rndssqres,rndexpvar,rndscaling,rndtuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, 'compmodes', compmodes, opt{:});
      case 'parafac2'
        [rndcomp,dum,rndssqres,rndexpvar,rndscaling,rndtuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, specmodes, 'compmodes', compmodes, opt{:});
      case 'parafac2cp'
        [rndcomp,dum,rndssqres,rndexpvar,rndscaling,rndtuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, specmodes, 'compmodes', compmodes, opt{:}, 'ssqdatnoncp', ssqdatnoncp);
      case 'spacetime'
        [rndcomp,dum,rndssqres,rndexpvar,rndscaling,rndtuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, freq, 'Dmode', Dmode, opt{:});
      case 'spacefsp'
        [rndcomp,dum,rndssqres,rndexpvar,rndscaling,rndtuckcongr] = feval(['nwaydecomp_' model], dat, ncomp, 'Dmode', Dmode, opt{:});
      otherwise
        error('model not yet supported in automatic random starting')
    end
    randssqres(irand)    = rndssqres;
    randcomp{irand}      = rndcomp;
    randexpvar(irand)    = rndexpvar;
    randscaling{irand}   = rndscaling;
    randtuckcongr{irand} = rndtuckcongr;
    disp([dispprefix 'random start ' num2str(irand) ': decomposition of random initialization ' num2str(irand) ' of ' num2str(nrand) ' finished: exp. var. = ' num2str(rndexpvar,'%-2.2f') '%']);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


% NaN failsafe
if any(isnan(randexpvar))
  error('at least one random start result in NaN explained variance, likely cause is input numerical array')
end

% Create randomstat structure and set best possible start-values
% sort output by ssqres
[randssqres,sortorder] = sort(randssqres);
randexpvar          = randexpvar(sortorder);
randscaling         = randscaling(sortorder);
randcomp            = randcomp(sortorder);
randtuckcongr       = randtuckcongr(sortorder);


% find possible degenerate solutions and create degeneracy index
if ncomp==1 % if ncomp is 1, there can be no degeneracy
  maxtg     = zeros(1,nrand);
  degeninit = [];
else
  maxtg     = cellfun(@max,(randtuckcongr));
  degeninit = find(maxtg >= degencrit);
end
% create index of non-degenerate initializations
if ~(length(degeninit)==nrand)
  initindex = 1:nrand;
  initindex(degeninit) = [];
else
  initindex = [];
end

% prepare for computing congruence between random initializations
% compute congruence coeffcients for all possible pairs of components from between all possible random starts
% using all possible pairs is necessary, because the order of components might not be exactly similar over random starts
% (this is the case when some components have nearly equal 'strength')
% computing between all possible random starts is not necessary, but done out of convenience
nparam     = numel(randcomp{1});
congrallrp = zeros(nrand,nrand,ncomp,nparam);
for irand1 = 1:nrand
  for irand2 = 1:nrand
    if irand1 == irand2
      % set congr to NaN when computing congr between the same initializations (to avoid contaminating the average later)
      congrallrp(irand1,irand2,:,:) = NaN(ncomp,nparam);
    else
      % set current estcomp
      currcomp = cell(1,2);
      currcomp{1} = randcomp{irand1};
      currcomp{2} = randcomp{irand2};
      % compute component congruence for all possible pairs between selected rands
      compcongr = zeros(ncomp,ncomp,length(currcomp{1}));
      for iparam = 1:nparam
        % perform model specific stuff
        switch model
          case {'parafac','parafac2','parafac2cp'}
            paramc1 = currcomp{1}{iparam};
            paramc2 = currcomp{2}{iparam};
            % normalize
            paramc1 = bsxfun(@rdivide,paramc1,sqrt(sum(abs(paramc1).^2,1)));
            paramc2 = bsxfun(@rdivide,paramc2,sqrt(sum(abs(paramc2).^2,1)));
            % put in compsrc
            compcongr(:,:,iparam) = abs(paramc1' * paramc2);
          case {'spacetime','spacefsp'}
            switch iparam
              case {1,2}
                paramc1 = currcomp{1}{iparam};
                paramc2 = currcomp{2}{iparam};
                % normalize
                paramc1 = bsxfun(@rdivide,paramc1,sqrt(sum(abs(paramc1).^2,1)));
                paramc2 = bsxfun(@rdivide,paramc2,sqrt(sum(abs(paramc2).^2,1)));
                % put in compsrc
                compcongr(:,:,iparam) = abs(paramc1' * paramc2);
              case 3
                paramc1 = currcomp{1}{iparam};
                paramc2 = currcomp{2}{iparam};
                % normalize
                paramc1 = bsxfun(@rdivide,paramc1,sqrt(sum(abs(paramc1).^2,1)));
                paramc2 = bsxfun(@rdivide,paramc2,sqrt(sum(abs(paramc2).^2,1)));
                % put in compsrc
                if size(paramc1,1) == size(paramc2,1)
                  compcongr(:,:,iparam) = abs(paramc1' * paramc2);
                else
                  compcongr(:,:,iparam) = 0; % congruence can't be computed, set to maximally incongruent (0)
                end
              case 4
                switch model
                  case 'spacetime'
                    % create frequency specific phases weighted by spatial maps and frequency profiles
                    A1 = currcomp{1}{1};
                    A2 = currcomp{2}{1};
                    B1 = currcomp{1}{2};
                    B2 = currcomp{2}{2};
                    S1 = currcomp{1}{4};
                    S2 = currcomp{2}{4};
                    % normalize
                    A1 = bsxfun(@rdivide,A1,sqrt(sum(abs(A1).^2,1)));
                    A2 = bsxfun(@rdivide,A2,sqrt(sum(abs(A2).^2,1)));
                    B1 = bsxfun(@rdivide,B1,sqrt(sum(abs(B1).^2,1)));
                    B2 = bsxfun(@rdivide,B2,sqrt(sum(abs(B2).^2,1)));
                    % construct spatial phase maps
                    Scomp1 = ipermute(exp(1i*2*pi*bsxfun(@times,permute(S1,[3 1 2]),freq')),[3 1 2]);
                    Scomp2 = ipermute(exp(1i*2*pi*bsxfun(@times,permute(S2,[3 1 2]),freq')),[3 1 2]);
                    % scale with A
                    Scomp1 = bsxfun(@times,Scomp1,A1);
                    Scomp2 = bsxfun(@times,Scomp2,A2);
                    % compute splitrelcoef over freqs,
                    srcoverfreq = zeros(ncomp,ncomp,size(B1,1));
                    for ifreq = 1:size(B1,1)
                      srcoverfreq(:,:,ifreq) = abs(Scomp1(:,:,ifreq)'*Scomp2(:,:,ifreq));
                    end
                    % weight with average B and combine over freq
                    Bweight    = bsxfun(@times,permute(B1',[1 3 2]),permute(B2',[3 1 2]));
                    Bweight    = bsxfun(@rdivide,Bweight,sum(Bweight,3));
                    srcsumfreq = sum(srcoverfreq .* Bweight,3);
                    % put in compsrc
                    compcongr(:,:,iparam) = srcsumfreq;
                  case 'spacefsp'
                    % create frequency specific phases weighted by spatial maps and frequency profiles
                    A1 = currcomp{1}{1};
                    A2 = currcomp{2}{1};
                    B1 = currcomp{1}{2};
                    B2 = currcomp{2}{2};
                    L1 = currcomp{1}{4};
                    L2 = currcomp{2}{4};
                    % normalize
                    A1 = bsxfun(@rdivide,A1,sqrt(sum(abs(A1).^2,1)));
                    A2 = bsxfun(@rdivide,A2,sqrt(sum(abs(A2).^2,1)));
                    B1 = bsxfun(@rdivide,B1,sqrt(sum(abs(B1).^2,1)));
                    B2 = bsxfun(@rdivide,B2,sqrt(sum(abs(B2).^2,1)));
                    % construct spatial phase maps
                    Lcomp1 = exp(1i*2*pi*permute(L1,[1 3 2]));
                    Lcomp2 = exp(1i*2*pi*permute(L2,[1 3 2]));
                    % scale with A
                    Lcomp1 = bsxfun(@times,Lcomp1,A1);
                    Lcomp2 = bsxfun(@times,Lcomp2,A2);
                    % compute splitrelcoef over freqs,
                    srcoverfreq = zeros(ncomp,ncomp,size(B1,1));
                    for ifreq = 1:size(B1,1)
                      srcoverfreq(:,:,ifreq) = abs(Lcomp1(:,:,ifreq)'*Lcomp2(:,:,ifreq));
                    end
                    % weight with average B and combine over freq
                    Bweight    = bsxfun(@times,permute(B1',[1 3 2]),permute(B2',[3 1 2]));
                    Bweight    = bsxfun(@rdivide,Bweight,sum(Bweight,3));
                    srcsumfreq = sum(srcoverfreq .* Bweight,3);
                    % put in compsrc
                    compcongr(:,:,iparam) = srcsumfreq;
                end
              case 5
                switch Dmode
                  case 'identity'
                    % D is fixed with arbitrary order, make its congruence coefficient irrelevant
                    compcongr(:,:,iparam) = 1;
                  case 'kdepcomplex'
                    B1 = currcomp{1}{2};
                    B2 = currcomp{2}{2};
                    D1 = currcomp{1}{5};
                    D2 = currcomp{2}{5};
                    % weight with B
                    D1 = bsxfun(@times,D1,permute(B1,[1 3 2]));
                    D2 = bsxfun(@times,D2,permute(B2,[1 3 2]));
                    % vectorize
                    paramc1 = reshape(permute(D1,[3 1 2]),[ncomp numel(B1)]).';
                    paramc2 = reshape(permute(D2,[3 1 2]),[ncomp numel(B2)]).';
                    % normalize
                    paramc1 = bsxfun(@rdivide,paramc1,sqrt(sum(abs(paramc1).^2,1)));
                    paramc2 = bsxfun(@rdivide,paramc2,sqrt(sum(abs(paramc2).^2,1)));
                    % put in compsrc
                    compcongr(:,:,iparam) = abs(paramc1' * paramc2);
                end
            end
        end
      end
      % get cong coefficients by selecting most-similair unique pairings
      compcongrsel = zeros(ncomp,nparam);
      congrsum     = sum(compcongr,3);
      % match from perspective of first rand (i.e. find components of rand 2 that match those of rand 1)
      % do so by starting from the component-pair with the highest similarity, then the next most similar, etc.
      r1ind   = zeros(1,ncomp);
      r2ind   = zeros(1,ncomp);
      for icomp = 1:ncomp
        [dum, r1ind(icomp)] = max(max(congrsum,[],2));
        [dum, r2ind(icomp)] = max(congrsum(r1ind(icomp),:));
        congrsum(r1ind(icomp),:) = 0;
        congrsum(:,r2ind(icomp)) = 0;
      end
      % sanity check
      if any(diff(sort(r1ind))==0) || any(diff(sort(r2ind))==0)
        error('some components were selected multiple times')
      end
      % sort for convenience
      [r1ind, sortind] = sort(r1ind);
      r2ind = r2ind(sortind);
      for iparam = 1:nparam
        compcongrsel(:,iparam) = diag(compcongr(r1ind,r2ind,iparam));
      end
      % save compsrc
      congrallrp(irand1,irand2,:,:) = compcongrsel;
    end
  end % irand2
end % irand1

% get congruence coeficient over ALL initializations (including possible degenerate ones)
congrall = squeeze(nanmean(nanmean(congrallrp,1),2));


% get congruence coeficient over SELECTION of initializations, disregarding those that might be degenerate and those not thought to be at the global minimum
% select which initilizations to calculate congruence over
if ~isempty(initindex)
  nanexpvar = randexpvar; % create a randexpvar vector with NaNs for degenerate initializations
  nanexpvar(degeninit) = NaN;
  expvarchange = abs(nanexpvar - nanexpvar(initindex(1)));
  globmininit  = find((expvarchange < .1) & (maxtg < degencrit)); % criterion means initializations are included if they differ less than 0.1% in expvar to the first (sorted) non-degenerate expvar
  % calculate per component per param the abs(inner product) (if only degenerate solutions are found than matrix will contain only NaNs)
  if numel(globmininit) ~= 1
    congrglobmin = squeeze(nanmean(nanmean(congrallrp(globmininit,globmininit,:,:),1),2));
  else
    congrglobmin = ones(ncomp,nrand);
  end
else
  globmininit = [];
  congrglobmin = [];
end

% compute cumulative component congruence, between the first random start, the first two, the first three, etc
% this is done using non-denegerate initializations only
congrcumul = zeros(numel(initindex),ncomp,nparam);
for iset = 1:numel(initindex)
  ind = initindex(1:iset);
  congrcumul(iset,:,:) = squeeze(nanmean(nanmean(congrallrp(ind,ind,:,:),1),2));
end

% put together randomstat
randomstat.expvar         = randexpvar;
randomstat.error          = randssqres;
randomstat.globmininit    = globmininit;
randomstat.congrall       = congrall;
randomstat.congrglobmin   = congrglobmin;
randomstat.congrcumul     = congrcumul;
randomstat.tuckcongr      = randtuckcongr;
randomstat.degeninit      = degeninit;
randomstat.scaling        = randscaling;
% add settings for randstart
randomstat.nrand          = nrand;
randomstat.convcrit       = convcrit;
randomstat.niter          = niter;



% put the scaling coefficients back in to make the startvalues proper startvalues
for irand = 1:nrand
  currrndcomp = randcomp{irand};
  currrndscal = randscaling{irand};
  
  switch model
    
    case {'parafac','parafac2','parafac2cp'}
      % magnitude scaling
      for icomp = 1:ncomp
        % set mode1
        mode1     = currrndcomp{1}(:,icomp);
        mode1norm = currrndscal{1}(icomp);
        mode1     = mode1 .* mode1norm;
        currrndcomp{1}(:,icomp) = mode1;
      end
      % phase scaling
      compmodesindex = find(compmodes==1);
      if sum(compmodes)~=0 % only do this if there is a complex mode
        for icomp = 1:ncomp
          % set mode1
          mode1      = currrndcomp{compmodesindex(1)}(:,icomp);
          phaseshift =  currrndscal{2}(icomp);
          mode1      = mode1 ./ phaseshift;
          currrndcomp{compmodesindex(1)}(:,icomp) = mode1;
        end
      end
      
    case {'spacetime','spacefsp'}
      % magnitude scaling
      for icomp = 1:ncomp
        % set mode1
        param3     = currrndcomp{3}(:,icomp);
        param3norm = currrndscal(icomp);
        param3     = param3 .* param3norm;
        currrndcomp{3}(:,icomp) = param3;
      end
      
    otherwise
      error('unsupported model')
  end
  
  % put back in randcomp
  randcomp{irand} = currrndcomp;
end

% get best possible solution and display number of expected degenerate solutions
disp([dispprefix 'random start: ' num2str(length(degeninit)) ' out of ' num2str(nrand) ' initializations expected to be degenerate'])
if ~isempty(globmininit)
  startval = randcomp{globmininit(1)};
  disp([dispprefix 'random start: lowest error-term from procedure = ' num2str(randssqres(globmininit(1)))])
else % if only degenerate solutions are found, than better pick the "best" one
  startval = randcomp{1};
  disp([dispprefix 'random start: lowest error-term from procedure = ' num2str(randssqres(1))])
  disp([dispprefix 'random start: warning: best random initialization might be degenerate'])
end

% save startvalues
randomstat.startvalglobmin = startval;
randomstat.startvalall     = randcomp;
























