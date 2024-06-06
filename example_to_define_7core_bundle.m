% Define parameters
core_radius = 8.2 / 2; % Core radius in micrometers
cladding_radius = 80 / 2; % Cladding radius in micrometers
delta_n = 0.01; % Relative refractive index difference
n_cladding = 1.444; % Cladding refractive index
n_background = 1.43; % Background refractive index

n_core = n_cladding * (1 + delta_n); % Core refractive index
fiber_radius = 3 * cladding_radius; % Estimate for the whole fiber radius
fiber_diameter = 2 * fiber_radius;

% Create a grid
grid_size = 2000; % Number of grid points for better resolution
x = linspace(-fiber_radius, fiber_radius, grid_size);
y = linspace(-fiber_radius, fiber_radius, grid_size);
[X, Y] = meshgrid(x, y);

% Initialize refractive index profile
n_profile = n_background * ones(grid_size);

% Define core positions for hexagonal arrangement
d = 2 * cladding_radius; % Center-to-center distance between adjacent cores
core_positions = [
    0, 0;
    d, 0;
    -d, 0;
    d/2, d * sqrt(3)/2;
    d/2, -d * sqrt(3)/2;
    -d/2, d * sqrt(3)/2;
    -d/2, -d * sqrt(3)/2
];

% Create refractive index profile for each core and cladding
for i = 1:size(core_positions, 1)
    cx = core_positions(i, 1);
    cy = core_positions(i, 2);
    
    % Set core region
    core_mask = ((X - cx).^2 + (Y - cy).^2) <= core_radius^2;
    n_profile(core_mask) = n_core;
    
    % Set cladding region
    cladding_mask = ((X - cx).^2 + (Y - cy).^2) <= cladding_radius^2 & ~core_mask;
    n_profile(cladding_mask) = n_cladding;
end

% Plot the refractive index profile
figure;
imagesc(x, y, n_profile);
axis equal;
colorbar;
title('Refractive Index Profile of 7-Core Hexagonal Fiber');
xlabel('x (µm)');
ylabel('y (µm)');
colormap(jet);
