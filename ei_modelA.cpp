/*
Compile by:

mkoctfile -g -v -O CFLAGS="\$CFLAGS -std=c++11" -lm -lgsl -lgslcblas -I/home/maurizio/Dropbox/Ongoing.Projects/DePitta.PNAS/Stability.Analysis -I/home/maurizio/Dropbox/Ongoing.Projects/c_libraries -I/home/maurizio/Dropbox/Ongoing.Projects/c_libraries/models /home/maurizio/Dropbox/Ongoing.Projects/c_libraries/dcomplex_libc.c /home/maurizio/Dropbox/Ongoing.Projects/c_libraries/special_functions_libc.c /home/maurizio/Dropbox/Ongoing.Projects/c_libraries/prob_dist_libc.c /home/maurizio/Dropbox/Ongoing.Projects/c_libraries/models/freq_cv_libc.c asn_meanfield_libc.c ei_modelA.cpp -output ei_modelA
mkoctfile -g -v -O -I/home/maurizio/Dropbox/Ongoing.Projects/DePitta.PNAS/Software/Stability.Analysis -I/home/maurizio/Dropbox/Ongoing.Projects/c_libraries -I/home/maurizio/Dropbox/Ongoing.Projects/c_libraries/models /home/maurizio/Dropbox/Ongoing.Projects/c_libraries/dcomplex_libc.c /home/maurizio/Dropbox/Ongoing.Projects/c_libraries/special_functions_libc.c /home/maurizio/Dropbox/Ongoing.Projects/c_libraries/prob_dist_libc.c /home/maurizio/Dropbox/Ongoing.Projects/c_libraries/models/freq_cv_libc.c /home/maurizio/Dropbox/Ongoing.Projects/DePitta.PNAS/Software/Stability.Analysis/asn_meanfield_libc.c ei_modelA.cpp -output ei_modelA -lm -lgsl -lgslcblas

*/

#include <octave/oct.h>

#include "asn_meanfield_libc.h"

DEFUN_DLD(ei_modelA,args, , "EI MF network model A with delays (Brunel, JCN 2000)"){

    // Input arguments
    NDArray xvar = args(0).array_value();
    NDArray pars = args(1).array_value();

    // Output arguments
    double dy[4];
    dim_vector dv (4,1);
    NDArray dxvars(dv);

    // Assignment to individual variables for readability
    double m_e    = xvar(0);
    double s2_e   = xvar(1);
    double m_i    = xvar(2);
    double s2_i   = xvar(3);
    double m_e_D  = xvar(4);
    double s2_e_D = xvar(5);
    double m_i_D  = xvar(6);
    double s2_i_D = xvar(7);

    // Assignment of parameters;
    int    C     = (int)pars(0);
    double gamma = pars(1);
    double tau   = pars(2);
    double trp   = pars(3);
    double vl    = pars(4);
    double vr    = pars(5);
    double vt    = pars(6);
    double u0    = pars(7);
    double j     = pars(8);
    double g     = pars(9);
    double r     = pars(10);
    double D     = pars(11);

    // Computation of derivatives
    rhs_ei_modelA(dy,
                  m_e, s2_e, m_i, s2_i,
                  m_e_D, s2_e_D, m_i_D, s2_i_D,
                  C, gamma,
                  tau, trp, vl, vr, vt,
                  u0, j, g, r, D);

    for(int i=0; i<4; i++) dxvars(i) = dy[i];

    // Must cast output array as octave_value
    return octave_value (dxvars);
}