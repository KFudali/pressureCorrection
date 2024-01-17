function L = assembleLaplace(nx, ny, hx, hy)
    diag_x_upper = diag(ones((nx-2)*(ny-2) - 1,1),1);
    diag_x_upper(nx-2:nx-2:end) = 0;
    diag_x_bottom = diag(ones((nx-2)*(ny-2) - 1,1),-1);
    diag_x_bottom(1:nx-2:end) = 0;
    lx = -2 * eye((nx-2)*(ny-2)) + diag_x_upper + diag_x_bottom;

    diag_y_upper = diag(ones((nx-2)*(ny-2) - (nx-2) ,1),nx-2);
    diag_y_bottom =  diag(ones((nx-2)*(ny-2) -(nx-2) ,1), -(nx-2));
    ly = -2 * eye((nx-2)*(ny-2)) + diag_y_upper + diag_y_bottom;
 
    L = (1/hx^2 * lx +  1/hy^2 * ly);
end