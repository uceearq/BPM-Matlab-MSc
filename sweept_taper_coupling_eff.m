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
% P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.

%% Resolution-related parameters (check for convergence)
P.Lx_main = 200e-6;        % [m] x side length of main area
P.Ly_main = 200e-6;        % [m] y side length of main area
P.Nx_main = 100;          % x resolution of main area
P.Ny_main = 100;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
% P.dz_target = 5e-7; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area


%% first calculate properties of the multimode fibre
P.lambda = 1550e-9; % [m] Wavelength
P.n_background = 1.435; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.4586; % [] reference refractive index
% P.Lz = 2e-3; % [m] z propagation distances for this segment

[X,Y] = ndgrid(single(P.x),single(P.y));
P.n.n = single(calcRI_mmf(X,Y,P.n_background));

P.n.Lx = P.dx*size(P.n.n,1);
P.n.Ly = P.dy*size(P.n.n,2);
P.n.xSymmetry = P.xSymmetry;
P.n.ySymmetry = P.ySymmetry;

figure(123)
mesh(P.n.n)

P = findModes(P,10);

neff_mmf = [P.modes.neff];
for k1 = 1:10
    field_mmf(:,:,k1) = real([P.modes(k1).field]);
end

figure(256)
for k1 = 1:9
    subplot(3,3,k1)
    mesh(field_mmf(:,:,k1)); view([0 0 1])
end
title('multimode fibre modes')

dx_mmf = P.dx; 
dy_mmf = P.dy;
X_mmf = X;
Y_mmf = Y;

for i1 = 1:10
    for i2 = 1:10
        eta(i1,i2) = sum(sum((field_mmf(:,:,i1).*conj(field_mmf(:,:,i2)))));
    end
end
disp('eta is approximately a identity matrix because the modes are the fundamental solutions of the waveguide, thereby orthogonal')

%% Problem definition for photonic lantern
P.lambda = 1550e-9; % [m] Wavelength
P.n_background = 1.435; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.4586; % [] reference refractive index
% P.Lz = 2e-3; % [m] z propagation distances for this segment

taper = [1:30];

for k11 = 1:length(taper)
    %% tapered by a factor of 2
    [X,Y] = ndgrid(single(P.x),single(P.y));
    P.n.n = single(calcRI_7bundle(X,Y,P.n_background,taper(k11)));

    P.n.Lx = P.dx*size(P.n.n,1);
    P.n.Ly = P.dy*size(P.n.n,2);
    P.n.xSymmetry = P.xSymmetry;
    P.n.ySymmetry = P.ySymmetry;


    figure(123)
    mesh(P.n.n)

    % We search for the 10 modes with effective refractive index closest to
    % n_0:
    P = findModes(P,10);

    figure(123)
    for k1 = 1:9
        subplot(3,3,k1)
        mesh(abs(P.modes(k1).field)); view([0 0 1])
    end

    neff_pl(:,k11) = [P.modes.neff];
    for k1 = 1:10
        field_pl(:,:,k1,k11) = real([P.modes(k1).field]);
    end


    [X,Y] = ndgrid(single(P.x),single(P.y));
    for i1 = 1:10
        field_mmf_resized(:,:,i1) = interpn(X_mmf,Y_mmf,field_mmf(:,:,i1),X,Y);
    end


    figure(512)
    for i1 = 1:9
        subplot(3,3,i1)
        mesh((field_mmf_resized(:,:,i1))); view([0 0 1])
    end


    eta(1:10,1:10,k11) = zeros(10);
    for i1 = 1:10
        for i2 = 1:10
            if neff_mmf(i1) < 1.444
                continue;
            end
            if neff_pl(i2) < P.n_background+1e-4
                continue;
            end
            eta(i1,i2,k11) = sum(sum((field_mmf_resized(:,:,i1).*conj(field_pl(:,:,i2,k11)))));
        end
    end


end

for k11 = 1:length(taper)
    aux = abs(eta(:,:,k11));
    coupeff(k11) = sum(sum(aux(1:6,1:6)));
end

%% photonic lantern coupling efficiency for different taper ratios
figure()
plot(taper,coupeff)
xlabel('taper ratio')
ylabel('Coupling Efficiency')


%% modes effective index
figure()
plot(taper,neff_pl)
ylabel('neff')
xlabel('taper ratio')

error('increase number of points to get smoother curves')