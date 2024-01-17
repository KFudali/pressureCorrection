% description: 2D Poisson equation solution using finite difference dscretization
% author: kfudali
% date: 14.03.2023

clear all
clc


%create domain mesh
Lx = 1;
Ly = 1;
nx = 4;
ny = 4;


x = [0:(Lx/(nx-1)):Lx];
hx = Lx / nx;

y = [0:(Ly/(ny-1)):Lx];
hy = Ly / ny;

[x,y] = meshgrid(x,y);
points = [x(:), y(:)];


u = zeros((nx)*(ny),1);



%Neumann
% nx = nx + 1;
% u = zeros((nx-2)*(ny-2),1);
% 
 f = zeros(size(u));
% boundary_value = 0;
% f(1:(nx-2))  = 2 * boundary_value / hy;
% neumann_ids = [1:(nx-2); ny:1:ny+(nx-3)]';
% 
% % f = zeros(size(u));
% f(1:(nx-2):end) = -1/hy^2;
L = assembleLaplace(nx+2,ny+2,hx,hy);
f = ones(size(u));

    neumann_ids_top = [1:nx; nx+1:2*nx];
    neumann_ids_bottom  = [nx*(nx - 1):nx^2; nx*(nx - 2):nx*(nx - 1)];
    neumann_ids_right = [1:nx:nx^2; [1:nx:nx^2] + 1];
    neumann_ids_left = [nx:nx:nx^2; [nx:nx:nx^2] - 1];

    neumann_ids = [neumann_ids_top, neumann_ids_bottom, neumann_ids_right, neumann_ids_left]';
    interior_ids = ~ismember(1:(nx*ny),neumann_ids(:,1));

    
    ids = sub2ind(size(L),neumann_ids(:,1), neumann_ids(:,2));
    L(ids) = 2 * L(ids);

% f(1:(nx-2)) = -1/hx^2;
% u = L\f;


%add time discretization
dt =0.01;
t_end = 1;
t = 0;
n_time_steps = t_end / dt;
figure;
clf

u_old = u;


while t < t_end
    u = u_old(:,end) + dt * L\f;
    U = reshape(u,nx-2,ny-2);
    contourf(U)
    drawnow
    t = t + dt;
end


% 
% U = reshape(u,nx-2,ny-2);
% contourf(U)
% colorbar


