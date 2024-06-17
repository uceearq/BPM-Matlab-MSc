clear all
close all

P = BPMmatlab.model;

% This example shows that modes can be found and stored in the parameters
% struct using the findModes function, then used as initial E fields. Here,
% we use the LP31o-like mode of a non-trivial refractive index profile.

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.

%% Resolution-related parameters (check for convergence)
P.Lx_main = 200e-6;        % [m] x side length of main area
P.Ly_main = 200e-6;        % [m] y side length of main area
P.Nx_main = 200;          % x resolution of main area
P.Ny_main = 200;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
% P.dz_target = 5e-7; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1550e-9; % [m] Wavelength
P.n_background = 1.435; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.46; % [] reference refractive index
% P.Lz = 2e-3; % [m] z propagation distances for this segment

%% no taper
[X,Y] = ndgrid(single(P.x),single(P.y));
P.n.n = single(calcRI_7bundle(X,Y,P.n_background,1));

P.n.Lx = P.dx*size(P.n.n,1);
P.n.Ly = P.dy*size(P.n.n,2);
P.n.xSymmetry = P.xSymmetry;
P.n.ySymmetry = P.ySymmetry;


figure()
mesh(P.n.n)

% We search for the 10 modes with effective refractive index closest to
% n_0:
P = findModes(P,10);

figure(123)
for k1 = 1:10
    subplot(5,5,k1)
    mesh(abs(P.modes(k1).field)); %view([0 0 1])
end

neff(:,1) = [P.modes.neff];

%% explaining how to taper
nBefore = single(calcRI_7bundle(X,Y,P.n_background,1));
nAfter  = single(calcRI_7bundle(X,Y,P.n_background,2));

figure()
subplot(1,2,1)
mesh(X,Y,nBefore); view([0 0 1])
subplot(1,2,2)
mesh(X,Y,nAfter); view([0 0 1])

%% tapered by a factor of 2
[X,Y] = ndgrid(single(P.x),single(P.y));
P.n.n = single(calcRI_7bundle(X,Y,P.n_background,2));

P.n.Lx = P.dx*size(P.n.n,1);
P.n.Ly = P.dy*size(P.n.n,2);
P.n.xSymmetry = P.xSymmetry;
P.n.ySymmetry = P.ySymmetry;


figure()
mesh(P.n.n)

% We search for the 10 modes with effective refractive index closest to
% n_0:
P = findModes(P,10);

figure(256)
for k1 = 1:9
    subplot(3,3,k1)
    mesh(abs(P.modes(k1).field)); view([0 0 1])
end

neff(:,2) = [P.modes.neff];

%% modes effective index
figure()
plot(1:20,neff)

figure()
plot([1 2],neff.')
