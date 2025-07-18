%% Advanced Transfer Functions synthesizer
% Interface script to Advanced Transfer Functions synthesizer

%% Path definition and variable declaration
warning off all;
clearvars;
%close all;
%clc;

% Definition of the working directory path
% Path    =   'C:\Users\1424366\OneDrive - UAB\MATLAB\Global-Tools-GitHub\Non_Standard_Functions_Synthesizer';
% addpath(genpath(Path));

% Variable pre-declaration
Topology    =   struct;
Polyn       =   struct;

%% Filter Specification
% Type of Transfer Function
    % Standard Transfer Function
        %   0 -> Generalized Chebyshev (GC)
    % Advanced Transfer Functions
        %   1 -> Non-Equiripple (NE)
        %   2 -> Reduced Chebyshev (RC): Complex Reflection Zeros (RZs)
        %   3 -> Bounded Chebyshev (BC)
        %   4 -> Non-Equiripple & Reduced Chebyshev (NE-RC)
Topology.TransferFunction       =   2;

% Polynomial Synthesis Methodology
    %   0 -> Recursive Method (based on Cameron's book)
    %   1 -> Iterative Method (based on Remez Algorithm)
Topology.SynthesisMethodology   =   1;

% Nature of Band
    %   0 -> LowPass (LP), HighPass (HP), BandPass (BP)
    %   1 -> BandStop (BS)
Topology.Nature     =   0;

% Multi-Band Response
    %   0 -> Single-Band
    %   1 -> Multi-Band
Topology.MultiBand  =   1;

%% Standard Parameters Definition
% Filter order
Topology.N      =  6;

% Transmission Zeros (TZs) in w-domain
Topology.TZw    =   [1.12];

Topology.TZw    =   [-Topology.TZw , 0, Topology.TZw];

% Return Loss (RL) level [dB]
Topology.RL     =   18;

% Input Phase [deg]
Topology.InputPhase     =   0;
Topology.OutputPhase     =   0;
%% MultiBand Parameters Definition
% Number of bands
Topology.M      =   2;

% Order of each band
Topology.Nk     =   [Topology.N/2, Topology.N/2];

% Corners of each band

%Topology.wlim   =   [-1.0, -0.55;
%                     +0.55, +1.0];
FBW = 0.4;
% second option to choose FBW and to have the limits established
fL=roots([1 -(2+FBW^2) 1]);
fL=fL(fL<1);
Topology.wlim   =   [-1, -fL; fL, 1];




%% Non-Standard Parameters Definition
% Non-Equiripple parameter:
% Set the scaling factors 'alpha' of the characteristic function to modify
% the ripple of the passband.
%   - Note 1: alpha's should be equal to or less than 1.
%   - Note 2: A row-vector is defined for a single band response while a
%             matrix is defined with M-rows for a multi-band response.
% Examples:
% - Single-band: [1 0.4 0.2 1 0.2 0.4 1];
% - Dual-band: [1 0.5 1 0.2 1; 1 0.2 1 0.5 1];
% Topology.alpha  =   [1 0.4 0.2 0.4 1];
Topology.alpha  =   [1 1 1 1 1];
%Topologty.alpha     =   [-Topology.alpha , Topology.alpha]

% Reduced Chebyshev parameters:
% Set the location of RZs in w-domain (attention with 'j', s = jw).
% Examples:
%   - On the real axis: [-3j, 3j];
%   - Out-of-Band RZs: [5];
%   - Complex RZs: [5 - 3j];
%Topology.RZw    =   [-0.8480-0.01j 0.8480-0.01j 0.7287-0.01j -0.7287-0.01j];%, 0.5+2j, 0.2+3j];
%1 tz at origin
Topology.RZw    =   [3 -3];%3 odd couplings

Topology.RZw    =   [];

% Initialize an optimization that recomputes the scaling factor 'alpha'
% with the objective of recovering the equiripple in the passband.
Topology.OptimizeRipple     =   1;

% Set whether or not the scaling factor 'alpha' will maintain or not 
% symmetry in the equiripple recovery optimization.
% - Note: This flag is enabled by default for Multi-Band.
Topology.SymmetryAround0_RZ =   0;

% Set the number of bins into which the passband is divided to evaluate the
% ripple performance.
% - Note: 2001 points by default.
Topology.bins   =   2001;

% Bounded Chebyshev parameter:
% Set the maximum Return Loss (RLmax) level to bound the depth of RZs.
%   - RLmax must be higher than RL, at least in 0.5 dB.
Topology.RLmax  =   25;

%% Synthesizer Summon
%DataConsistency(Topology);
[ Polyn, ABCD ]     =   Compute_Polynomials(Topology);

% Plot of the S-Parameters (in dB)
%Plot_Polyn(Polyn);
%Plot_ABCD(ABCD);

[CM.Trans]      =   Get_CM_orig(Polyn,ABCD, Topology);
CM.Trans.M
%going to the inner
NetMap                      =   NetMapDescription2(CM.Trans);
[CM.Inner, NetMap]          =   MatrixReduction(CM.Trans, NetMap, 1, -1);
[CM.InnerTrans, NetMap]     =   DiagCM(CM.Inner, NetMap);
[CM.Order, eigenmodes]      =   ModesSortingCM(CM.InnerTrans);
CM.Order.M
%draw_CouplingMatrix(CM.Order,length(CM.Order.M))
values = (sign(CM.Order.M(1,:)).*sign(CM.Order.M(end,:)))
real_range=[0 2.5];
imag_range=[-0.6 0.6];
Tech=CM;
%Topology.RZw_out=Topology.RZw;

%Sweep_AWTech_Uniformity(real_range, imag_range, Topology, [], []);

[CM]   =   rotations_third_shunt(CM.Order);
[Elements_norm]   =   Elements_third_shunt(CM)
