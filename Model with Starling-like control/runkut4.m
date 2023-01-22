function value = runkut4(passo,A,x,B)
    u=1;
    xdot = (A*x) + (B*u);
    kx1 = passo*xdot;

    x1 = x + 0.5*kx1;
    xdot = (A*x1) + (B*u);
    kx2 = passo*xdot;

    x1 = x + 0.5*kx2;
    xdot = (A*x1) + (B*u);
    kx3 = passo*xdot;

    x1 = x + kx3;
    xdot = (A*x1) + (B*u);
    kx4 = passo*xdot;

    value = x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6;
end
    