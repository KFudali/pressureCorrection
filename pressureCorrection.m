% description: Cavity flow solution using finite difference dscretization
% author: kfudali
% date: 14.03.2023

clear all
clc

%create domain mesh
Lx = 1;
Ly = 1;
nx = 50;
ny = 50;

x = [0:(Lx/(nx-1)):Lx];
hx = Lx / nx;
y = [0:(Ly/(ny-1)):Lx];
hy = Ly / ny;

[x,y] = meshgrid(x,y);

ux = zeros((nx-2)*(ny-2),1);
uy = zeros((nx-2)*(ny-2),1);

u = zeros((nx-2)*(ny-2),1);

a0 = 3/2;
a1 = 2;
a2 = -1/2;

p = zeros((nx)*(ny),1);
phi = zeros(size(p));
p_star = zeros(size(p));

[Dx, Dy]  = assembleGrad(nx-2,ny-2,hx,hy);
[Dx_p, Dy_p]  = assembleGrad(nx,ny,hx,hy);

L = assembleLaplace(nx,ny,hx,hy);
L_p = assembleLaplace(nx+2,ny+2,hx,hy);

L_p_inv = L_p^-1;


%Boundary conditions
bottom_ids = 1:nx;
top_ids = nx*(ny-1) + 1:nx*ny;
left_ids = 1:nx:nx*ny;
right_ids = nx:nx:nx*ny;

boundary_ids = unique([bottom_ids, top_ids, left_ids, right_ids]);
interior_ids = setdiff(1:nx*ny,boundary_ids);

%Dirichlet conditions on velocity
%     top_velocity = cos(2*pi*[0:(Lx/(nx-1)):Lx]) - 1;
top_velocity = ones(1,nx);
    L_rhs = zeros(size(ux));
    
%Neumann conditions on pressure

    neumann_ids_top = [1:nx; nx+1:2*nx];
    neumann_ids_bottom  = [nx*(nx - 1):nx^2; nx*(nx - 2):nx*(nx - 1)];
    neumann_ids_right = [1:nx:nx^2; [1:nx:nx^2] + 1];
    neumann_ids_left = [nx:nx:nx^2; [nx:nx:nx^2] - 1];

    neumann_ids = [neumann_ids_top, neumann_ids_bottom, neumann_ids_right, neumann_ids_left]';
    interior_ids = ~ismember(1:(nx*ny),neumann_ids(:,1));


    ids = sub2ind(size(L_p),neumann_ids(:,1), neumann_ids(:,2));
    L_p(ids) = 2 * L_p(ids);


%add time discretization
    nu = 0.1;
    dt = 0.01;
    t_end = 3;
    t = 0;
    n_time_steps = t_end / dt;

    ux_old = zeros(size(ux,1), n_time_steps);
    uy_old = zeros(size(uy,1), n_time_steps);
    p_old = zeros(size(p,1), n_time_steps);
    phi_old = zeros(size(p,1), n_time_steps);
    t_step = 1;

    x_vec = reshape(y,nx*ny,1);
    y_vec = flip(reshape(x,nx*ny,1));
    xlabel('y');
    ylabel('x');

    iter = 1;
%main loop
while t < t_end

    if t_step < 3
        ux_old(:,t_step) = ux;
        uy_old(:,t_step) = uy;
        p_old(:,t_step) = p;
        phi_old(:,t_step) = phi;    
        t_step = t_step + 1;
        t = t +dt;
    else

    p = p_old(:,t_step-1) + 4/3 * phi_old(:,t_step-1) - 1/3 * phi_old(:,t_step-2);

    %extrapolate convection terms

    L_rhs(1:nx-2) = top_velocity(2:end-1)/hy^2; 

%     add BCs to velocity field in order to calculate pressure correction
    ux_full = [zeros(1,nx-2);  reshape(ux,(nx-2),(ny-2)); zeros(1,nx-2)];
    ux_full = [top_velocity', ux_full, zeros(ny,1)];
    ux_full = reshape(ux_full, nx * ny, 1);

    uy_full = [zeros(ny-2,1), reshape(uy,(nx-2),(ny-2)), zeros(ny-2,1)];
    uy_full = [zeros(1,nx); uy_full; zeros(1,nx)];
    uy_full = reshape(uy_full, nx * ny, 1);

    ux_star = ux_full;
    uy_star = uy_full;

    term_1 = 0.5 * diag(Dx_p * ux_star + Dy_p * uy_star);
    term_2 = ux_star.* Dx_p + uy_star.* Dy_p;

    term_1(unique([1:nx, (nx-1)*ny+1:nx*ny, 1:nx:nx*ny, nx:nx:nx*ny]),:) = [];
    term_1(:,unique([1:nx, (nx-1)*ny+1:nx*ny, 1:nx:nx*ny, nx:nx:nx*ny])) = [];
    term_2(unique([1:nx, (nx-1)*ny+1:nx*ny, 1:nx:nx*ny, nx:nx:nx*ny]),:) = [];
    term_2(:,unique([1:nx, (nx-1)*ny+1:nx*ny, 1:nx:nx*ny, nx:nx:nx*ny])) = [];

    
    %solve velocity
    lhs = a0/dt * eye(size(L)) - nu * L+ (term_1 + term_2);
    %

    Dxp = Dx_p * p;
    Dyp = Dy_p * p;

    rhs_x = - Dxp(interior_ids) + 1/dt * (a1 * ux_old(:,t_step- 1) + a2 * ux_old(:, t_step - 2))  + nu*L_rhs;
    rhs_y = - Dyp(interior_ids) + 1/dt * (a1 * uy_old(:,t_step- 1) + a2 * uy_old(:, t_step - 2));


    [ux, ~, ~, iter] = pcg(lhs,rhs_x);
    [uy, ~, ~, iter] = pcg(lhs,rhs_y);
%     ux = lhs\rhs_x;
%     uy = lhs\rhs_y;

    %add BCs to velocity field in order to calculate pressure correction\
    ux_full = [zeros(1,nx-2);  reshape(ux,(nx-2),(ny-2)); zeros(1,nx-2)];
    ux_full = [top_velocity', ux_full, zeros(ny,1)];
    ux_full = reshape(ux_full, nx * ny, 1);

    uy_full = [zeros(ny-2,1), reshape(uy,(nx-2),(ny-2)), zeros(ny-2,1)];
    uy_full = [zeros(1,nx); uy_full; zeros(1,nx)];
    uy_full = reshape(uy_full, nx * ny, 1);

    %solve phi to update pressure
    
    lhs = L_p;
    rhs = a0/dt * (Dx_p * ux_full + Dy_p * uy_full);

    phi = L_p_inv*rhs;
%     phi = pcg(L_p, rhs);

    p = p_old(:,t_step - 1) + phi - nu * (Dx_p * ux_full + Dy_p * uy_full);

    ux_old(:,t_step) = ux;
    uy_old(:,t_step) = uy;
    p_old(:,t_step) = p;
    phi_old(:,t_step) = phi;

    vel = (ux_full.^2 + uy_full.^2).^(0.5);
            contourf(x,y,reshape(p,nx,ny))
%  contourf(x,y,reshape(ux_full,nx,ny))
%  contourf(x,y,reshape(p,nx,ny))           
    title('vel')
    colorbar
    drawnow

    t = t + dt;
    t_step = t_step + 1;
    end
end

results.u = ux_full;
results.v = uy_full;
results.x = x_vec;
results.y = y_vec;
results.p = p;

save("results_stokes/results_nu_0001.mat","results")

