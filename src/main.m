<<<<<<< Updated upstream:src/main.m
=======
clear; close all; clc;

%% SET MODEL PARAMETERS

runID = 'test_DoubleDiff_TK1';
nop   = 20;       % print output every 'nop' steps
lvplt = 0;        % plot figures live (1) or in background (0)     
svfig = 1;        % save figures to file (1)

% set domain parameters
N     = 200;      % num. grid size
D     = 1e3;      % phys. domain depth [m]
h     = D./N;     % grid spacing [m] 

% set physical parameters

mu    = 1e-3;     % pore fluid viscosity (water) [Pa s]
a     = 5e-3;     % grain size of matrix (sandstone) [m]
b     = 100;      % geom. factor for permeability [1]
n     = 3;        % permeability powerlaw [1]
rhol0 = 1000;     % fluid density [kg/m3]
grav  = 9.81;     % gravity [m/s2]
kC    = 1e-8;     % chemical diffusivity [m2/s]  
kT    = 1e-6;     % thermal diffusivity [m2/s]
aT    = 2e-4;     % thermal expansivity [1/K]
gC    = 2;        % chemical density contrast [1/wt]

% set initial condition parameters
finit = 'linear'; % initial condition: 'linear' or 'layer'
f0    = 0.1;      % background porosity [vol]
f1    = 0.1;      % base porosity [vol]  
df    = 0.001;    % perturbation amplitude [vol]
Tinit = 'layer';  % initial condition: 'linear' or 'layer'
T0    = 50;       % top temperature [C]
T1    = 100;      % base temperature [C]
dT    = 0.1;      % perturbation amplitude [C]
Cinit = 'layer';  % initial condition: 'linear' or 'layer'
C0    = 0.005;    % top concentration  [C]
C1    = 0.01;     % base concentration [C]
dC    = 0.0;      % perturbation amplitude [C]
zlay  = 0.5;      % relative depth of layer boundary
smth  = (N/30)^2; % smoothness of random noise

% set model timing parameters
tend  = 1e11;     % model stopping time [s]

% set numerical solver parameters
CFL   = 0.25;     % Courant number to limit time step size
tol   = 1e-8;     % residual tolerance for iterative solver
alpha = 0.99;     % step size for iterative solver
beta  = 0.95;     % damping parameter for iterative solver

% create output directory
if ~isfolder(['./out/',runID])
    mkdir(['./out/',runID]);
end

>>>>>>> Stashed changes:darcyDoubleDiffusion_TK1.m
%% INITIAL MODEL SETUP

% initialise grid coordinates
h = D./N;     % grid spacing
x = linspace(-h/2,D+h/2,N+2);
z = linspace(-h/2,D+h/2,N+2);
[X,Z] = meshgrid(x,z);

% initialise smooth random noise
rng(5);
rn = rand(N+2,N+2) - 0.5;
for i=1:smth
   rn(2:end-1,2:end-1) = rn(2:end-1,2:end-1) ...
                       + diff(rn(:,2:end-1),2,1)./8 ...
                       + diff(rn(2:end-1,:),2,2)./8;
    rn([1 end],:) = rn([2 end-1],:);
    rn(:,[1 end]) = rn(:,[2 end-1]);
end
rn = rn./max(abs(rn(:)));

% set initial condition
switch finit  % initial porosity
    case 'linear'
        f = f0 + (f1-f0) .* Z/D + df.*rn;
    case 'layer'
        f = f0 + (f1-f0) .* (1+erf(25*(Z/D-zlay)))/2 + df.*rn;
end
switch Tinit  % initial temperature
    case 'linear'
        T = T0 + (T1-T0) .* Z/D + dT.*rn;
    case 'layer'
        T = T0 + (T1-T0) .* (1+erf(25*(Z/D-zlay)))/2 + dT.*rn;
end
switch Cinit  % initial concentration
    case 'linear'
        C = C0 + (C1-C0) .* Z/D + dC.*rn;
    case 'layer'
        C = C0 + (C1-C0) .* (1+erf(25*(Z/D-zlay)))/2 + dC.*rn;
end

DTDt = 0.*T(2:end-1,2:end-1);
DCDt = 0.*C(2:end-1,2:end-1);
dt   = 1e-6;

fh1 = figure(1);
subplot(3,1,1)
imagesc(x,z,f); axis equal tight; colorbar;
title('Initial Porosity [vol]')
subplot(3,1,2)
imagesc(x,z,T); axis equal tight; colorbar;
title('Initial Temperature [C]')
subplot(3,1,3)
imagesc(x,z,C); axis equal tight; colorbar;
title('Initial Concentration')
drawnow

% prepare solution & residual arrays for VP solver
w = zeros(N+1,N+2);  % vertical Darcy speed
u = zeros(N+2,N+1);  % horizontal Darcy speed
p = zeros(N+2,N+2);  % pore fluid pressure
F = zeros(N+2,N+2);  % residual for pressure equation


%% MAIN TIME STEPPING LOOP
m    = 0;
time = 0;
while time <= tend

    fprintf(1,'\n\n*****  step = %d,  time = %1.3e [s],  step = %1.3e \n\n',m,time,dt);

    % store previous rates of change
    To    = T;
    Co    = C;
    DTDto = DTDt;
    DCDto = DCDt;
    
    % calculate permeability [m2]
    k = a^2/b * f.^n;  % Kozeny-Carman relationship
    
    % calculate Darcy coefficient [m2/Pas]
    K = k/mu;
    
    % calculate iterative step size
    dtau = (h/2)^2./K;
        
    % update density difference
    Drho = - rhol0.*(- aT.*(T-mean(T,2)) + gC.*(C-mean(C,2)));
    
    
    % UPDATE VELOCITY-PRESSURE SOLUTION (PSEUDO-TRANSIENT SOLVER)
    Fnorm = 1e6;
    pi    = p;
    it    = 0;
    while Fnorm >= tol || it < 100
        
        % store previous iterative solution guesses
        pii = pi; pi = p;
        
        % calculate pressure gradient [Pa/m]
        gradPz = diff(p,1,1)./h;  % vertical gradient
        gradPx = diff(p,1,2)./h;  % horizontal gradient
        
        % calculate Darcy segregation speed [m/s]
        w = -(K(1:end-1,:)+K(2:end,:))./2 .* (gradPz + (Drho(1:end-1,:)+Drho(2:end,:))./2.*grav);
        u = -(K(:,1:end-1)+K(:,2:end))./2 .* (gradPx + 0                                       );
        
        % calculate residual of pressure equation
        F(2:end-1,2:end-1) = diff(w(:,2:end-1),1,1)./h + diff(u(2:end-1,:),1,2)./h;

        % update pressure solution
        p = pi - alpha.*F.*dtau + beta.*(pi-pii);
        
        % apply pressure boundary conditions
        p(:,1  ) = p(:,2    );
        p(:,end) = p(:,end-1);
        p(1  ,:) = p(2    ,:) + (Drho(1  ,:)+Drho(2    ,:))./2.*grav.*h;
        p(end,:) = p(end-1,:) - (Drho(end,:)+Drho(end-1,:))./2.*grav.*h;
        
        % get preconditioned residual norm to monitor convergence
        Fnorm = norm(F(:).*dtau(:),2)./norm(p(:),2);
        
        % print convergence
        if ~mod(it,100)
            
            % UPDATE TEMPERATURE SOLUTION (EXPLICIT SOLVER)
            dt = CFL .* min([(h/2)/max(abs(w(:))) , (h/2)/max(abs(u(:))) , (h/2)^2./kT]);  % diffusive timestep
            
            % calculate temperature advection
            wp = max(0, (w(1:end-1,2:end-1)+w(2:end,2:end-1))./2);
            wm = min(0, (w(1:end-1,2:end-1)+w(2:end,2:end-1))./2);
            up = max(0, (u(2:end-1,1:end-1)+u(2:end-1,2:end))./2);
            um = min(0, (u(2:end-1,1:end-1)+u(2:end-1,2:end))./2);
            
            grdTzp = diff(T(2:end  ,2:end-1),1,1)./h;
            grdTzm = diff(T(1:end-1,2:end-1),1,1)./h;
            grdTxp = diff(T(2:end-1,2:end  ),1,2)./h;
            grdTxm = diff(T(2:end-1,1:end-1),1,2)./h;
            
            DTDt = - (wp.*grdTzm + wm.*grdTzp + up.*grdTxm + um.*grdTxp) ...
                 + kT.* (diff(T(:,2:end-1),2,1)./h^2 + diff(T(2:end-1,:),2,2)./h^2);
            
            T(2:end-1,2:end-1) = To(2:end-1,2:end-1) + (DTDt + DTDto)/2 .* dt;
            
            % apply temperature boundary conditions
            T(:,1  ) = T(:,2    );  % left boundary: insulating
            T(:,end) = T(:,end-1);  % right boundary: unsulating
            T(1  ,:) = T0;          % top boundary: isothermal
            T(end,:) = T1;          % bottom boundary: isothermal
            
            
            % UPDATE CONCENTRATION SOLUTION (EXPLICIT SOLVER)
            
            % calculate concentration advection
            grdCzp = diff(C(2:end  ,2:end-1),1,1)./h;
            grdCzm = diff(C(1:end-1,2:end-1),1,1)./h;
            grdCxp = diff(C(2:end-1,2:end  ),1,2)./h;
            grdCxm = diff(C(2:end-1,1:end-1),1,2)./h;
            
            DCDt = - (wp.*grdCzm + wm.*grdCzp + up.*grdCxm + um.*grdCxp) ...
                + kC.* (diff(C(:,2:end-1),2,1)./h^2 + diff(C(2:end-1,:),2,2)./h^2);
            
            C(2:end-1,2:end-1) = Co(2:end-1,2:end-1) + (DCDt + DCDto)/2 .* dt;
            
            % apply concentration boundary conditions
            C(:,1  ) = C(:,2    );  % left boundary: closed
            C(:,end) = C(:,end-1);  % right boundary: closed
            C(1  ,:) = C0;          % top boundary: isochemical
            C(end,:) = C1;          % bottom boundary: isochemical
            
            % report convergence
            fprintf(1,'---  %d,  %e\n',it+1,Fnorm);
        end
        
        % increment iteration count
        it = it+1;
    end

    % plot solution
    if ~mod(m,nop)
        if lvplt
            fh2 = figure(2);
        else
            fh2=figure('Visible','off');
        end
        
        sgtitle(sprintf('Time elapsed %.1f years', time/31557600))
        subplot(2,3,1);
        imagesc(x,z,-w.*3600*24*365.25); axis equal tight; colorbar;
        title(['Segregation','z-speed [m/yr]'])
        subplot(2,3,2);
        imagesc(x,z,u.*3600*24*365.25); axis equal tight; colorbar;
        title(['Segregation','x-speed [m/yr]'])
        subplot(2,3,3);
        imagesc(x,z,p); axis equal tight; colorbar;
        title(['Dynamic fluid','pressure [Pa]'])
        subplot(2,3,4);
        imagesc(x,z,T); axis equal tight; colorbar;
        title('Temperature [C]')
        subplot(2,3,5);
        imagesc(x,z,C); axis equal tight; colorbar;
        title('Concentration [wt]')
        subplot(2,3,6);
        imagesc(x,z,f); axis equal tight; colorbar;
        title('Porosity [vol]')
        drawnow % draw figure now (avoids Matlab deciding to delay)
        
        % print figure to file
        if svfig
            print(fh2,['../out/',runID,'/',runID,'_',int2str(m/nop)],'-dpng','-r200')
        end
    end
    
    % update time and step count
    m    = m + 1;
    time = time + dt;  

end 
