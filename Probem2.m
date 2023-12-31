
% Number of vertices
N = 21;
% Time step size
dt = 0.01; % second
% Rod length
RodLength = 0.1; % meter
% Discrete length
deltaL = RodLength / (N-1);
% Radius of spheres
R = zeros(N,1);
R(:) = deltaL / 10;
midNode = (N+1)/2;
R(midNode) = 0.025;
% Density
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;
% Rod radius
r0 = 0.001;
% Young's modulus
Y = 1e9; % Using Y instead of E to avoid ambiguity
% Gravity
g = 9.8; % m/s^2
% Viscosity
visc = 1000; % Pa-s
% Total time
totalTime = 10; % seconds
% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;
% Geometry
nodes = zeros(N, 2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
    %     nodes(c,2) = 0;
end
% Mass matrix
M = zeros(2*N,2*N);
for k=1:N
    M(2*k-1,2*k-1) = 4/3*pi*R(k)^3*rho_metal;
    M(2*k,2*k) = M(2*k-1,2*k-1);
end
% Viscous damping matrix
C = zeros(2*N,2*N);
for k=1:N
    C(2*k-1,2*k-1) = 6*pi*visc*R(k);
    C(2*k,2*k) = C(2*k-1,2*k-1);
end
% Gravity
W = zeros(2*N,1);
for k=1:N
W(2*k-1) = 0;
    W(2*k) = -4/3*pi*R(k)^3*rho*g;
end
% Initial DOF vector
q0 = zeros(2*N,1);
for c=1:N
q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
    q0 ( 2*c ) = nodes(c,2); % y coordinate
end
% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector
% Number of time steps
Nsteps = round( totalTime / dt );
all_mid_y = zeros( Nsteps, 1); % y-position of R2 
all_mid_v = zeros( Nsteps, 1); % y-velocity of R2
all_mid_y(1) = q(2*midNode);
all_mid_v(1) = u(2*midNode);
% Tolerance
tol = EI / RodLength^2 * 1e-3;
% Time marching scheme
for c=2:Nsteps
    fprintf('Time = %f\n', (c-1) * dt );
    q = q0; % Guess
     % Newton Raphson
     err = 10 * tol;
     while err > tol
% Inertia
         f = M / dt * ( (q-q0) / dt - u );
         J = M / dt^2;
         %
         % Elastic forces
         %
         % Linear spring 1 between nodes 1 and 2
         for k=1:N-1
             xk = q(2*k-1);
             yk = q(2*k);
             xkp1 = q(2*k+1);
             ykp1 = q(2*k+2);
             dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
             dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
             ind = [2*k-1; 2*k; 2*k+1; 2*k+2];
             f(ind) = f(ind) + dF;
             J(ind,ind) = J(ind,ind) + dJ;
end

% Bending spring between nodes 1, 2, and 3
for k=2:N-1
xkm1 = q(2*k-3);
ykm1 = q(2*k-2);
xk = q(2*k-1);
yk = q(2*k);
xkp1 = q(2*k+1);
ykp1 = q(2*k+2);
curvature0 = 0;
dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
curvature0, deltaL, EI);
dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                 curvature0, deltaL, EI);
             ind = [2*k-3; 2*k-2; 2*k-1; 2*k; 2*k+1; 2*k+2];
             f(ind) = f(ind) + dF;
             J(ind,ind) = J(ind,ind) + dJ;
end
    % Viscous force
    f = f + C * ( q - q0 ) / dt;
    J = J + C / dt;
% Weight
f = f - W;
% Update
q = q - J \ f;
err = sum( abs(f) );
end
    % Update
    u = (q - q0) / dt; % Velocity
    q0 = q; % Old position
figure(1);
plot( q(1:2:end), q(2:2:end), 'ro-'); axis equal
drawnow
% Store
    all_mid_y(c) = q(22);
    all_mid_v(c) = u(22);
end
figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, 'k-'); 
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]');
figure(3);
plot(timeArray, all_mid_y, 'k-');
xlabel('Time, t [sec]');
ylabel('Position of mid-node, y [meter]');





%% This section is to plot terminal velocity vs. number of nodes
% Time step size
dt = 0.01; % second
% Rod length
RodLength = 0.1; % meter
% Density
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;
% Rod radius
r0 = 0.001;
% Young's modulus
Y = 1e9; % Using Y instead of E to avoid ambiguity
% Gravity
g = 9.8; % m/s^2
% Viscosity
visc = 1000; % Pa-s
% Total time
totalTime = 10; % seconds

for N=3:2:31
% Discrete length
deltaL = RodLength / (N-1);
% Radius of spheres
R = zeros(N,1);
R(:) = deltaL / 10;
midNode = (N+1)/2;
R(midNode) = 0.025;

% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;
% Geometry
nodes = zeros(N, 2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
    %     nodes(c,2) = 0;
end
% Mass matrix
M = zeros(2*N,2*N);
for k=1:N
    M(2*k-1,2*k-1) = 4/3*pi*R(k)^3*rho_metal;
    M(2*k,2*k) = M(2*k-1,2*k-1);
end
% Viscous damping matrix
C = zeros(2*N,2*N);
for k=1:N
    C(2*k-1,2*k-1) = 6*pi*visc*R(k);
    C(2*k,2*k) = C(2*k-1,2*k-1);
end
% Gravity
W = zeros(2*N,1);
for k=1:N
W(2*k-1) = 0;
    W(2*k) = -4/3*pi*R(k)^3*rho*g;
end
% Initial DOF vector
q0 = zeros(2*N,1);
for c=1:N
q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
    q0 ( 2*c ) = nodes(c,2); % y coordinate
end
% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector
% Number of time steps
Nsteps = round( totalTime / dt );
all_mid_y = zeros( Nsteps, 1); % y-position of R2 
all_mid_v = zeros( Nsteps, 1); % y-velocity of R2
all_mid_y(1) = q(2*midNode);
all_mid_v(1) = u(2*midNode);
% Tolerance
tol = EI / RodLength^2 * 1e-3;
% Time marching scheme
for c=2:Nsteps
    %fprintf('Time = %f\n', (c-1) * dt );
    q = q0; % Guess
     % Newton Raphson
     err = 10 * tol;
     while err > tol
% Inertia
         f = M / dt * ( (q-q0) / dt - u );
         J = M / dt^2;
         %
         % Elastic forces
         %
         % Linear spring 1 between nodes 1 and 2
         for k=1:N-1
             xk = q(2*k-1);
             yk = q(2*k);
             xkp1 = q(2*k+1);
             ykp1 = q(2*k+2);
             dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
             dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
             ind = [2*k-1; 2*k; 2*k+1; 2*k+2];
             f(ind) = f(ind) + dF;
             J(ind,ind) = J(ind,ind) + dJ;
end

% Bending spring between nodes 1, 2, and 3
for k=2:N-1
xkm1 = q(2*k-3);
ykm1 = q(2*k-2);
xk = q(2*k-1);
yk = q(2*k);
xkp1 = q(2*k+1);
ykp1 = q(2*k+2);
curvature0 = 0;
dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
curvature0, deltaL, EI);
dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                 curvature0, deltaL, EI);
             ind = [2*k-3; 2*k-2; 2*k-1; 2*k; 2*k+1; 2*k+2];
             f(ind) = f(ind) + dF;
             J(ind,ind) = J(ind,ind) + dJ;
end
    % Viscous force
    f = f + C * ( q - q0 ) / dt;
    J = J + C / dt;
% Weight
f = f - W;
% Update
q = q - J \ f;
err = sum( abs(f) );
end
    % Update
    u = (q - q0) / dt; % Velocity
    q0 = q; % Old position
% Store
    all_mid_y(c) = q(N+1);
    all_mid_v(c) = u(N+1);
end
for i=1:numel(3:2:31)
all_mid_v_terminal(i) = all_mid_v(end);
end
end
node_num=[3:2:31];
figure(1);
plot(node_num, all_mid_v_terminal, 'ko-');
xlabel('Number of nodes');
ylabel('Terminal velocity of mid-node, v [meter/s]');
ylim([-0.008 0]);




%% This section is to plot terminal velocity vs. time step size
% Assume number of nodes is fixed
N = 21;

% Rod length
RodLength = 0.1; % meter
% Discrete length
deltaL = RodLength / (N-1);
% Radius of spheres
R = zeros(N,1);
R(:) = deltaL / 10;
midNode = (N+1)/2;
R(midNode) = 0.025;
% Density
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;
% Rod radius
r0 = 0.001;
% Young's modulus
Y = 1e9; % Using Y instead of E to avoid ambiguity
% Gravity
g = 9.8; % m/s^2
% Viscosity
visc = 1000; % Pa-s
% Total time
totalTime = 10; % seconds
% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;
% Geometry
nodes = zeros(N, 2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
    %     nodes(c,2) = 0;
end
% Mass matrix
M = zeros(2*N,2*N);
for k=1:N
    M(2*k-1,2*k-1) = 4/3*pi*R(k)^3*rho_metal;
    M(2*k,2*k) = M(2*k-1,2*k-1);
end
% Viscous damping matrix
C = zeros(2*N,2*N);
for k=1:N
    C(2*k-1,2*k-1) = 6*pi*visc*R(k);
    C(2*k,2*k) = C(2*k-1,2*k-1);
end
% Gravity
W = zeros(2*N,1);
for k=1:N
W(2*k-1) = 0;
    W(2*k) = -4/3*pi*R(k)^3*rho*g;
end
% Initial DOF vector
q0 = zeros(2*N,1);
for c=1:N
q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
    q0 ( 2*c ) = nodes(c,2); % y coordinate
end
% New position and velocity
q = q0; % DOF vector

for dt=0.0001:0.001:0.01
u = (q - q0) / dt; % Velocity vector
% Number of time steps
Nsteps = round( totalTime / dt );
all_mid_y = zeros( Nsteps, 1); % y-position of R2 
all_mid_v = zeros( Nsteps, 1); % y-velocity of R2
all_mid_y(1) = q(2*midNode);
all_mid_v(1) = u(2*midNode);
% Tolerance
tol = EI / RodLength^2 * 1e-3;
% Time marching scheme
for c=2:Nsteps
    q = q0; % Guess
     % Newton Raphson
     err = 10 * tol;
     while err > tol
% Inertia
         f = M / dt * ( (q-q0) / dt - u );
         J = M / dt^2;
         %
         % Elastic forces
         %
         % Linear spring 1 between nodes 1 and 2
         for k=1:N-1
             xk = q(2*k-1);
             yk = q(2*k);
             xkp1 = q(2*k+1);
             ykp1 = q(2*k+2);
             dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
             dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
             ind = [2*k-1; 2*k; 2*k+1; 2*k+2];
             f(ind) = f(ind) + dF;
             J(ind,ind) = J(ind,ind) + dJ;
end

% Bending spring between nodes 1, 2, and 3
for k=2:N-1
xkm1 = q(2*k-3);
ykm1 = q(2*k-2);
xk = q(2*k-1);
yk = q(2*k);
xkp1 = q(2*k+1);
ykp1 = q(2*k+2);
curvature0 = 0;
dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
curvature0, deltaL, EI);
dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                 curvature0, deltaL, EI);
             ind = [2*k-3; 2*k-2; 2*k-1; 2*k; 2*k+1; 2*k+2];
             f(ind) = f(ind) + dF;
             J(ind,ind) = J(ind,ind) + dJ;
end
    % Viscous force
    f = f + C * ( q - q0 ) / dt;
    J = J + C / dt;
% Weight
f = f - W;
% Update
q = q - J \ f;
err = sum( abs(f) );
end
    % Update
    u = (q - q0) / dt; % Velocity
    q0 = q; % Old position

% Store
    all_mid_y(c) = q(22);
    all_mid_v(c) = u(22);
end
for i=1:numel(0.0001:0.001:0.01)
    all_mid_v_terminal(i) = all_mid_v(end);
end
end
figure(1);
time_step = [0.0001:0.001:0.01];
plot(time_step, all_mid_v_terminal, 'ko-'); 
xlabel('Time step, t [sec]');
ylabel('Terminal velocity of mid-node, v [meter/sec]');
ylim([-0.008 0]);
