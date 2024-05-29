function P = findModes(P,nModes,varargin)
singleCoreModes = false;
sortByLoss = false;
plotModes = true;
for i=1:2:numel(varargin)
  switch varargin{i}
    case 'singleCoreModes'
      singleCoreModes = varargin{i+1};
    case 'sortByLoss'
      sortByLoss = varargin{i+1};
    case 'plotModes'
      plotModes = varargin{i+1};
    otherwise
      error('Argument #%d for findModes not recognized. The syntax for findModes has recently changed. See the examples and the readme.',i+2);
  end
end

if P.xSymmetry ~= 0 && ~isinf(P.bendingRoC) && sind(P.bendDirection) || ...
   P.ySymmetry ~= 0 && ~isinf(P.bendingRoC) && cosd(P.bendDirection)
  error('The specified bending direction is inconsistent with the symmetry assumption');
end

dx = P.dx;
dy = P.dy;
Nx = P.Nx;
Ny = P.Ny;
Lx = P.Lx;
Ly = P.Ly;
x = P.x;
y = P.y;

N = Nx*Ny;   % N*N - size of sparse matrices
if nModes >= N - 1
  error('Error: The number of modes requested must be less than the pixels in the full simulation window minus one (roughly Nx_main*padfactor*Ny_main*padfactor - 1)');
end

[X,Y] = ndgrid(x,y);

k_0 = 2*pi/P.lambda;  % [m^-1] Wavenumber

fprintf('Finding modes...\n');
tic

[Nx_source,Ny_source,Nz_source] = size(P.n.n);
dx_source = P.n.Lx/Nx_source;
dy_source = P.n.Ly/Ny_source;
x_source = getGridArray(Nx_source,dx_source,P.n.ySymmetry);
y_source = getGridArray(Ny_source,dy_source,P.n.xSymmetry);
[x_source,y_source,n_source] = calcFullRI(x_source,y_source,P.n.n);
if Nz_source == 1 % If P.n.n is 2D, interpolate it to the simulation grid
  n = interpn(x_source,y_source,n_source,x,y.','linear',P.n_background);
else % Otherwise it's 3D so we take the first slice and interpolate to the simulation grid
  n = interpn(x_source,y_source,n_source(:,:,1),x,y.','linear',P.n_background);
end

anycomplex = ~isreal(n);

if singleCoreModes
  coreIdxs = findCores(n,P.n_background);
else
  coreIdxs = ones(size(n));
end
nCores = max(coreIdxs(:));

% P.modes = []; % Clear any existing modes
iMode = 0;
for iCore = 1:nCores
  n_core = n;
  n_core(coreIdxs ~= iCore) = P.n_background;
  
  if ~isfinite(P.bendingRoC)
    [radiallySymmetric,xC,yC] = testRadialSymmetry(X,Y,n_core,P.n_background,P.xSymmetry,P.ySymmetry); % Numerically estimates whether this core is radially symmetric, and finds the centroid coordinates (xC, yC)
  else
    radiallySymmetric = false;
  end
  
  n_core = double(n_core); % Needed for sparse operations
  
  n_core_bent = real(n_core).*(1-(real(n_core).^2.*(X*cosd(P.bendDirection) + Y*sind(P.bendDirection))*P.rho_e/(2*P.bendingRoC))).*exp((X*cosd(P.bendDirection) + Y*sind(P.bendDirection))/P.bendingRoC);
  
  delta_n_2 = real(n_core_bent).^2 - P.n_0^2;                              %delta_n^2 in the expression for FD BPM

  dz = 1e-10;
  
  xEdge = P.Lx_main*(1 + (P.ySymmetry ~= 0))/2;
  yEdge = P.Ly_main*(1 + (P.xSymmetry ~= 0))/2;
  absorber = exp(-dz*(max(0,max(0,max(abs(Y) - yEdge,abs(X) - xEdge))).^2*P.alpha + 2*pi*imag(n_core)/P.lambda)); % First part is edge absorber, second is absorption from the imaginary part of the refractive indices
  ax = 1.00001*dz/(dx^2*2i*k_0*P.n_0);
%   ax = dz/(dx^2*2i*k_0*P.n_0);
  ay = dz/(dy^2*2i*k_0*P.n_0);

  M_rhs = sparse(1:N,1:N,absorber(1:N) + delta_n_2(1:N)*dz*k_0/(2i*P.n_0),N,N) + ...
    sparse(1:N-1,2:N,[repmat([repmat(ax,1,Nx-1) 0],1,Ny-1) repmat(ax,1,Nx-1)],N,N) + ...
    sparse(2:N,1:N-1,[repmat([repmat(ax,1,Nx-1) 0],1,Ny-1) repmat(ax,1,Nx-1)],N,N) + ...
    sparse(1:N-Nx,1+Nx:N,ay,N,N) + ...
    sparse(1+Nx:N,1:N-Nx,ay,N,N);
  M_rhs(1:N+1:N*N) = M_rhs(1:N+1:N*N) - repmat([ax repmat(2*ax,1,Nx-2) ax],1,Ny);
  M_rhs(1:N+1:N*N) = M_rhs(1:N+1:N*N) - [repmat(ay,1,Nx) repmat(2*ay,1,Nx*(Ny-2)) repmat(ay,1,Nx)];
  absorberdiag = sparse(1:N,1:N,absorber(1:N),N,N);
  M_rhs = M_rhs*absorberdiag;
  if P.xSymmetry == 2
    M_rhs(1:N+1:N*Nx) = 0;
    M_rhs(N*Nx:N+1:2*N*Nx) = 0;
  end
  if P.ySymmetry == 2
    M_rhs(1:((N+1)*Nx):end) = 0;
    M_rhs(N:((N+1)*Nx):end) = 0;
  end

  [V,D] = eigs(M_rhs,ceil(nModes/nCores),1,'Display',false,'SubspaceDimension',min(N,ceil(nModes/nCores)*10));
  D = diag(D);
  
  kappa = (1-real(D))/dz*P.lambda/(4*pi);
  realn = sqrt(P.n_0^2 - 2*P.n_0*imag(D/(dz*k_0)));

  neff = realn;%realn(realn > P.n_background) + 1i*kappa((realn > P.n_background));
  % V = V(:,realn > P.n_background);
  % D = D(realn > P.n_background);
  
  for iCoreMode = 1:numel(D)
    iMode = iMode + 1;
    P.modes(iMode).Lx = Lx;
    P.modes(iMode).Ly = Ly;
    P.modes(iMode).xSymmetry = P.xSymmetry;
    P.modes(iMode).ySymmetry = P.ySymmetry;
    E = reshape(V(:,iCoreMode),[Nx Ny]);
    E = E.*exp(-1i*angle(max(E(:)))); % Shift phase arbitrarily so that the resulting modes are (nearly) real
    P.modes(iMode).field = E;
    if anycomplex
      P.modes(iMode).neff = neff(iCoreMode);
    else
      P.modes(iMode).neff = real(neff(iCoreMode)); % If the user did not specify any complex refractive indices, then the only contribution to kappa would be from the edge absorber, and showing a complex neff would just confuse people and not be very physically meaningful
    end
    if radiallySymmetric
      [x_full,y_full,E_full] = calcFullField(x,y,E);
      [X_full,Y_full] = ndgrid(x_full,y_full);
      [~,iMax] = max(E_full(:));
      xMax = X_full(iMax);
      yMax = Y_full(iMax);
      theta = atan2(yMax - yC,xMax - xC);
      radialE = interpn(X_full,Y_full,E_full,xC + linspace(0,max(Lx,Ly),1000)*cos(theta),yC + linspace(0,max(Lx,Ly),1000)*sin(theta));
      radialEpruned = radialE(abs(radialE) > 0.1*max(abs(radialE)));
      m = sum(abs(diff(angle(radialEpruned*exp(1i*pi/2)) > 0))) + 1;

      R = sqrt((xMax - xC)^2 + (yMax - yC)^2);
      azimuthalE = interpn(X_full,Y_full,E_full,xC + R*cos(theta + linspace(0,2*pi,1000)),yC + R*sin(theta + linspace(0,2*pi,1000)));
      azimuthalEpruned = azimuthalE(abs(azimuthalE) > 0.1*max(abs(azimuthalE)));
      l = sum(abs(diff(angle(azimuthalEpruned*exp(1i*pi/2)) > 0)))/2;

      if l > 0
        Emaxmirrored = interpn(X_full,Y_full,E_full,xMax,2*yC - yMax);
        if real(E_full(iMax)/Emaxmirrored) < 0
          parity = 'o';
        else
          parity = 'e';
        end
      else
        parity = '';
      end
      P.modes(iMode).label = [', LP' num2str(l) num2str(m) parity];
    else
      P.modes(iMode).label = '';
    end
  end
end

if isempty(P.modes)
  fprintf('\b Done, %.1f seconds elapsed.\nNo guided modes found.\n',toc);
  return
end

if sortByLoss
  [~,sortedidxs] = sort(imag([P.modes.neff]),'ascend');
else
  [~,sortedidxs] = sort(real([P.modes.neff]),'descend');
end
P.modes = P.modes(sortedidxs(1:min(numel(P.modes),nModes)));

fprintf('\b Done, %.1f seconds elapsed.\n%d guided modes found.\n',toc,numel(P.modes));

% % % for iMode = 1:numel(P.modes)
% % %   P.modes(iMode).label = ['Mode ' num2str(iMode) P.modes(iMode).label];
% % %   if plotModes
% % %     E = P.modes(iMode).field;
% % %     h_f = figure(100+iMode);
% % %     h_f.WindowStyle = 'docked';
% % %     subplot(1,2,1);
% % %     imagesc(x,y,abs(E.').^2);
% % %     axis equal; axis tight; axis xy;
% % %     setColormap(gca,P.intensityColormap);
% % %     subplot(1,2,2);
% % %     imagesc(x,y,angle(E.'),'AlphaData',max(0,(1+log10(abs(E.'/max(E(:))).^2)/3)));
% % %     set(gca,'Color',0.7*[1 1 1]);  % To set the color corresponding to phase outside the cores where there is no field at all
% % %     caxis([-pi pi]);
% % %     axis equal; axis tight; axis xy;
% % %     setColormap(gca,P.phaseColormap);
% % %     neff = P.modes(iMode).neff;
% % %     if ~verLessThan('matlab','9.5')
% % %       if anycomplex
% % %         sgtitle({[P.modes(iMode).label ', n_{eff} = ' num2str(real(neff),'%.6g') ' + ' num2str(imag(neff),'%.3g') 'i'],['rough loss estimate: ' num2str(imag(neff)*4*pi/P.lambda,'%.3g') ' m^{-1} (' num2str(-10*log10(exp(-1))*imag(neff)*4*pi/P.lambda,'%.3g') ' dB/m)']});
% % %       else
% % %         sgtitle({[P.modes(iMode).label ', n_{eff} = ' num2str(real(neff),'%.6g')],['rough loss estimate: ' num2str(imag(neff)*4*pi/P.lambda,'%.3g') ' m^{-1} (' num2str(-10*log10(exp(-1))*imag(neff)*4*pi/P.lambda,'%.3g') ' dB/m)']});
% % %       end
% % %     end
% % %   end
% % % end
% % % drawnow;
% % % end
% % % 
% % % function setColormap(gca,colormapType)
% % % switch colormapType
% % %   case 'GPBGYR'
% % %     colormap(gca,GPBGYRcolormap);
% % %   case 'HSV'
% % %     colormap(gca,hsv/1.5);
% % %   case 'Parula'
% % %     colormap(gca,parula);
% % %   case 'Gray'
% % %     colormap(gca,gray);
% % %   case 'Cividis'
% % %     colormap(gca,cividisColormap);
% % % end
% % % end