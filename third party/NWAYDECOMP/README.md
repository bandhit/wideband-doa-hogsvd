**NWAYDECOMP** is a MATLAB toolbox intended for N-way/multi-way decompositions of neuroscience data. It contains several models that decompose numerical arrays, which are complex-valued or real-valued, into components. What kind of neuroscience data each decomposition requires, how a component is modeled, what a components represents, as well as the required structure of the numerical array, depends on the model. The toolbox depends on functions from other toolboxes/releases (most notably FieldTrip), which are described at the end of this document. 

The toolbox contains the following functions:
-	nd_nwaydecomposition
-	nwaydecomp_spacefsp
-	nwaydecomp_spacetime
-	nwaydecomp_parafac
-	nwaydecomp_parafac2
-	nwaydecomp_parafac2cp

The first function, *nd_nwaydecomposition* is the main interface function, and is an intelligent wrapper that repeatedly calls the individual model functions. It takes care of randomly initializing the models to avoid local minima of their loss functions, and to avoid degenerate solutions. It additionally contains algorithms for estimating the number of components to extract from the numerical array, which needs to be determined empirically (similar to ICA/PARAFAC).  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Nd_nwaydecomposition requires ‘FieldTrip style’ input. That is, a configuration structure (MATLAB structure array) which contains input options in its fields, and a FieldTrip style data structure (MATLAB structure array). The numerical array to be decomposed should be present in a field of this data structure, as well as at least a ‘label’ and ‘dimord’ field.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The output of this function is a MATLAB structure, which contains a ‘comp’ field. This field is a MATLAB cell-array which contains an extracted component in each cell, which then contains the component-specific parameters in separate lower-level cells. Additional output fields are described in the documentation of each individual function (at the top).  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For those users novice to FieldTrip, the ‘label’ field is a MATLAB cell-array of strings, where each cell contains the name for each site/channel in the data. The ‘dimord’ field is a string, which contains an abbreviation of the names of each of the dimensions of the numerical array. See below for an example.

The other functions are ‘low-level’ model functions that extract components according to their models. These take as input the numerical array to be decomposed, a set of required options, and additional options as key-value pairs. The functions are repeatedly called by the high-level main function, but can also be called individually.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The output of these functions is different from the high-level function, and is separated into multiple output arguments. The most important one, by default ‘comp’, contains the extracted components. It is different from the ‘comp’ field described above. It is a MATLAB cell-array, but now it contains in each cell the individual parameter matrices of the model that was used. With a few exceptions, these matrices contain in their columns the individual component loading-vectors. Output-specifics are described in each of the functions’ documentation (at the top).  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The models and their usage are briefly discussed below, with the appropriate references. Please cite these when you use them.
<br><br><br>
**SPACE-FSP/SPACE-time**  
These two models decompose Fourier coefficients obtained from electrophysiological recordings at multiple sites, multiple frequencies, multiple epochs/trials, and using multiple tapers. Importantly, the models extract components according to a plausible model of a neurobiological rhythm: a spatially distributed oscillation with energy in a range of frequencies and between-site phase relations for each frequency. The first model, SPACE-FSP extracts between-site phase relations for each frequency. The second model, SPACE-time, fits between-site time delays to these phase relations over frequencies. The models, their components, and the rationale behind them, are described in detail in the publications below (please cite this when used):  
* van der Meij R, Jacobs J, Maris E (2015). *Uncovering phase-coupled oscillatory networks in electrophysiological data.* Human Brain Mapping  
* van der Meij R, Van Ede F, Maris E (2016). *Rhythmic Components in Extracranial Brain Signals Reveal Multifaceted Task Modulation of Overlapping Neuronal Activity.* PLOS One 
  
  
The two models require specific input. The numerical array needs to be complex-valued and contain Fourier coefficients over sites, frequencies, trials/epochs, and tapers, in dimensions *in that order*, in a *4-way numerical array*. As such, the ‘dimord’ field in the data-structure needs to be ‘chan_freq_epoch_tap’. Importantly, the term tapers here is loosely defined. It can refer to e.g. Welch windows in a single trial, multiple Slepian tapers, or a combination of both. For example, if you use ‘mtmconvol’ as a spectral analysis method in FieldTrip to obtain temporally resolved Fourier coefficients, and use of these time-point to compute a controlled estimate of the between-site cross-products, then each of the time-points is a Welch-window (which describes the spectral content of a time-window surrounding it). The dimension containing tapers can also contain NaNs that depend on trial/epoch and frequency; i.e. the models can deal with a variable number of tapers. Input-specifics are described in the documentation of the functions.
<br><br><br>
**PARAFAC/PARAFAC2/PARAFAC2CP**  
These functions extract components according to the PARAFAC/2 models, and can be used on any N-way array (complex- and real-valued; N>2). These models are modified in one way: they allow for parameter matrices to be real-valued when the numerical array is complex-valued. The PARAFAC model was used to extracted components from a complex-valued 4-way array which described phase-amplitude coupling between oscillations at amplitude-providing and phase-providing sites (first 2 dimensions), and at amplitude-providing and phase-providing frequencies (last 2 dimensions). Only the first two parameter matrices, corresponding to the first two dimensions, were complex-valued, and described all phase information in the array. This was published in the publication below (please cite this when you use this model):  
* van der Meij R, Kahana M, Maris E (2012). *Phase-Amplitude Coupling in Human Electrocorticography Is Spatially Distributed and Phase Diverse.* Journal of Neuroscience 32: 111-123.  

The PARAFAC2 models have been modified in the same way, and allow for some matrices to be real-valued (not all of them, see function documentation). When using them, please cite the publication above. The PARAFAC2CP model is nearly identical to the PARAFAC2 model, except that it allows for a more computational efficient usage (with the consequence that the final ‘P’ needs to be computed by the user). The models will be merged in the future. The models are expanded in one of two possible ways when they are used on N-way arrays where N>3 (the original model was specified for N=3). This is specified in the function documentation.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;An extremely valuable resource on PARAFAC/2 and other N-way models is written by Rasmus Bro, which can be found at http://www.models.kvl.dk/users/rasmus/ (the monograph, highly recommended). This site also contains a very valuable code base, which, amongst others, contains the N-way Toolbox (a very important inspiration for my implementation of the PARAFAC models described above). Other important work which has been of great importance to the development of the above PARAFAC models, but also for both SPACE models, is that of Kiers & Bro on PARAFAC2 and Sidiropoulos, Giannakis & Bro on complex-valued PARAFAC. The appropriate references (and others) are provided below.  
<br><br><br>
**Dependencies**  
This toolbox depends on *FieldTrip*, an open-source toolbox for advanced electrophysiological data-analysis. FieldTrip can be downloaded from:  
*The FieldTrip wiki*: http://fieldtriptoolbox.org  
*The FieldTrip github repository*: https://github.com/fieldtrip/fieldtrip  

*Nwaydecomp* also depends on the following contribution from the MATLAB file exchange:
Khatri-Rao product (fast): http://www.mathworks.com/matlabcentral/fileexchange/28872  
(*(used in all models)*
This files is included in the repository as external.
IMPORTANT, if the N-way toolbox from Rasmus Bro (fantastic work) is also installed, then there might be a function naming conflict, as it also contains a function called kr.m. The function described above is a faster implementation, can perform successive khatri-rao products, and is required for the toolbox described here, nwaydecomp. )
<br><br><br>
**Reference literature on PARAFAC/2**
* Bro R (1998) Multi-way Analysis in the Food Industry. Models, Algorithms, and Applications. In: 
	Universiteit van Amsterdam.
* Carrol JD, Chang J (1970) Analysis of individual differences in multidimensional scaling via an N-way
 generalization of "Eckart-Young" decomposition. Psychometrika 35.
* Harshman RA (1970) Foundations of the PARAFAC procedure: model and conditions for an 'explanatory' 
	multi-mode factor analysis. UCLA Working Papers in Phonetics 16:1-84.
* Kiers HAL, Ten Berge JMF, Bro R (1999) PARAFAC2 - Part I. A direct fitting algorithm for the PARAFAC2 
	model. Journal of Chemometrics 13:275-294.
* Sidiropoulos ND, Giannakis GB, Bro R (2000) Blind PARAFAC receivers for DS-CDMA systems. Ieee T 
	Signal Proces 48:810-823. 

<br><br>
**Copyright statement**  
Copyright (C) 2012-2016 (or present), Roemer van der Meij, roemervandermeij AT gmail DOT com

This file is part of Nwaydecomp, see https://github.com/roemervandermeij/nwaydecomp
Nwaydecomp is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

Nwaydecomp is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. You should have received a copy of the GNU General Public License along with Nwaydecomp. If not, see http://www.gnu.org/licenses/.
