function dy = ei_modelA_sys_rhs()
      
    xv = [12., 12.;
          1., 0.9;
          15., 13.;
          0.5, 1.1];
  
    C     = 1000;
    gamma = 0.25;
    tau   = 0.02;
    trp   = 2e-3;
    vl    = 0.;
    vr    = 10.;
    vt    = 20.;
    u0    = 0.3;
    j     = 0.1;
    g     = 5;
    r     = 2;
    D     = 1.5e-3;
  
    pars = [C,gamma,tau,trp,vl,vr,vt,u0,j,g,r,D];
    dy = ei_modelA(xv,pars);
endfunction  