% Number of nodes
N = 50;
ndof = N * 2; % number of degrees of freedom
dt = 0.01; % second - Time step size
RodLength = 1; %
deltaL = RodLength / (N-1);
midNode = (N+1)/2;

% Radii of spheres
R = 0.013; 
r = 0.011;

% External force matrix, P
P_location=round(0.75/deltaL);
P = zeros(ndof,1);
P(P_location*2,1)=2000;

% Density
rho_metal = 2700; % kg/m^3
r0 = 0.001; % meter - rod radius
Y = 7e10; % Young's modulus (Y instead of E for clarity)
g = 9.8; % m/s^2 - gravity
m_beam=pi*(R^2-r^2)*RodLength*rho_metal;

totalTime = 1; % second - total simulation time

% Utility parameter
ne = N - 1; % number of edges
EI = Y * pi * r0^4 / 4; % Nm^2 - bending stiffness
EA = Y * pi * r0^2; % Newton

% Geometry - initial configuration
nodes = zeros(N,2);
for c=1:N % Loop over all the nodes
nodes(c,1) = (c-1) * deltaL; % x coordinates
    nodes(c,2) = 0;
end

% Mass, M
M = zeros(ndof, ndof);
for k=1:N
M(2*k-1, 2*k-1) = m_beam/(N-1); % Mass for x_k
M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
end

% Weight vector, W
W = zeros(ndof, 1);
for k=1:N
    W(2*k-1) = 0; % weight along x is zero
    W(2*k) = m_beam*g/N ;
end


% Initial DOF
q0 = zeros(ndof, 1);
for c=1:N % loop over nodes
    q0( 2*c-1 ) = nodes(c,1); % x1, x2, xn
    q0( 2*c )   = nodes(c,2); % y1, y2, yn
end
u0 = zeros(ndof, 1); % old velocity (initial velocity)

% tolerance
tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected
% Time marching scheme
Nsteps = round(totalTime/dt);
% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);

% Fixed and free DOFs
fixedDOF = [1;2; ndof];
freeDOF = 3:ndof-1;
boundaryConditionVector = [0;0;0];

for c = 2:Nsteps
    fprintf('Time = %f\n', (c-1) * dt);

% Guess
q = q0; % New DOFs are initialized to be equal to old DOFs 
q(fixedDOF) = boundaryConditionVector; 
    % Newton Raphson
    err = 10 * tol;
    while err > tol
        
        q_free = q(freeDOF);
        f = M / dt * ( (q-q0)/dt - u0 );
        J = M / dt^2;
         %
         % Elastic forces
         %
        % Linear spring 1 between nodes k and k+1
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
       % Bending spring between nodes k, k+1, and k+2
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
         
         f = f - W - P;
         
         % At this point, we have f and J
         f_free = f(freeDOF);
         J_free = J(freeDOF, freeDOF);
        % Update
         dq_free = J_free \ f_free;
         q_free = q_free - dq_free;
         err = sum ( abs(f_free) );
         q(freeDOF) = q_free;
        
end
    % New velocity
    u = (q - q0) / dt;
    % Store some information
    all_mid_v(c) = u(2*midNode);

% Update (new becomes old)
q0 = q;
u0 = u;

% Plot
figure(1);
plot( q(1:2:end), q(2:2:end), 'ro-'); axis equal
xlabel('x [meter]');
ylabel('y [meter]');
drawnow


end
% Plot middle node downward velocity
figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity (vertical) of middle node, v [m/s]');

