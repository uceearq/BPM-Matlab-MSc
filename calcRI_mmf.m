%% USER DEFINED RI FUNCTIONS
function n_profile = calcRI_mmf(X,Y,n_background)

% Define parameters
core_radius = 17.2 / 2 * 1e-6; % Core radius in micrometers
cladding_radius = 125 / 2 * 1e-6; % Cladding radius in micrometers
delta_n = 0.00466; % Relative refractive index difference
n_cladding = 1.444; % Cladding refractive index
% n_background = 1.43; % Background refractive index

n_core = n_cladding * (1 + delta_n); % Core refractive index
fiber_radius = 3 * cladding_radius; % Estimate for the whole fiber radius
% fiber_diameter = 2 * fiber_radius;

% % Create a grid
% grid_size = 2000; % Number of grid points for better resolution
% x = linspace(-fiber_radius, fiber_radius, grid_size);
% y = linspace(-fiber_radius, fiber_radius, grid_size);
% [X, Y] = meshgrid(x, y);

% Initialize refractive index profile
n_profile = n_background * ones(size(X));

% Create refractive index profile for each core and cladding

    % Set core region
    core_mask = ((X ).^2 + (Y ).^2) <= core_radius^2;
    n_profile(core_mask) = n_core;

    % Set cladding region
    cladding_mask = ((X ).^2 + (Y ).^2) <= cladding_radius^2 & ~core_mask;
    n_profile(cladding_mask) = n_cladding;







end
