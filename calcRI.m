%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,taper)
    X = X.*taper;
    Y = Y.*taper;
    % n may be complex
    n = n_background*ones(size(X)); % Start by setting all pixels to n_background
    % Core 1 is step index:
    n((X + 2.5e-6).^2 + Y.^2 < 5e-6^2) = 1.46;
    
    % Core 2 is graded index:
    corepos = [5e-6 0];
    r = 2.5e-6;
    R = sqrt((X - corepos(1)).^2 + (Y - corepos(2)).^2);
    n(R < r) = n_background + (1.47 - n_background)*(1 - (R(R < r)/r).^2); % Equation for parabolic graded index core
end
