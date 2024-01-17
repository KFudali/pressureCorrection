function [Dx, Dy] = assembleGrad(Nx, Ny, dx, dy)
%     nx = Nx+2;
%     ny = Ny+2;
%     Dx = (diag(ones((nx-2)*(ny-2),1)) - diag(ones((nx-2)*(ny-2)-1,1),-1)) / (dx);
%     Dy = (diag(ones((nx-2)*(ny-2),1)) - diag(ones((nx-2)*(ny-2) -(nx-2) ,1), -(nx-2))) / (dy);

    u_xd = diag(ones(1,Nx*Nx-1),1);
    u_xd(Nx:Nx:end) = 0;
    
    b_xd = -diag(ones(1,Nx*Ny-1),-1);
    b_xd((Nx+1):Nx:end) = 0;
    
    Hx = 1/(2*dx)*(u_xd + b_xd);
    
    
    u_yd = diag(ones(1,Nx*Ny-Ny),Ny);
    b_yd = -diag(ones(1,Nx*Ny-Ny),-Ny);
    
    Hy = 1/(2*dy)*(u_yd + b_yd);

    Dx = Hx;
    Dy = Hy;
end