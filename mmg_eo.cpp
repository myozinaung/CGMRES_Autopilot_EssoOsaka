#include "mmg_eo.hpp"

void StateODE::MMG_EO_ODE(const double t, const double* x, const double* u, double* f)
{
    // get variables

    double psi   = x[2];
    double u_vel = x[3];
    double v_vel = x[4];
    double r     = x[5];

    // Call Physical Component functions
    double tau_hull[3], tau_prop[3], tau_rudd[3]; // 3 DOF
    hull_EO(t, x, u, tau_hull);
    prop_EO(t, x, u, tau_prop);
    rudd_EO(t, x, u, tau_rudd);

    double X = tau_hull[0] + tau_prop[0] + tau_rudd[0];
    double Y = tau_hull[1] + tau_prop[1] + tau_rudd[1];
    double N = tau_hull[2] + tau_prop[2] + tau_rudd[2]; 

     // state odes
     f[0] = u_vel*cos(psi) - v_vel*sin(psi);
     f[1] = u_vel*sin(psi) + v_vel*cos(psi);
     f[2] = r;
     f[3] = (X + m*v_vel*r)/(m+mx);
     f[4] = (Y - m*u_vel*r)/(m+my);
     f[5] = (N - xG*Y)/(I_Jzz);

    //f[0] = u_vel*cos(psi) - v_vel*sin(psi);
    //f[1] = u_vel*sin(psi) + v_vel*cos(psi);
    //f[2] = r;
    //f[3] = (X + (m+my)*v_vel*r +xG*m*r*r)/(m+mx);
    //f[4] = (Y - m*u_vel*r)/(m+my);
    //f[5] = (N - xG*Y)/(I_Jzz);
}

// MMG Hull, Propeller, Rudder Components Functions required by StateFunc
void StateODE::hull_EO(const double t, const double* x, const double* u, double* tau_hull)
{
    // get variables
    double u_vel = x[3];
    double v_vel = x[4];
    double r     = x[5];

    if (abs(u_vel) < 1e-8) // to avoid singularity
    {
        u_vel = 1e-8;
    }
    double U    = sqrt(u_vel*u_vel + v_vel*v_vel);
    if (abs(U) < 1e-8) {
        U = 1e-8;
    }
    double beta;
    if (abs(u_vel) < 1e-8) {// to avoid singularity
        beta = 0;
    } else {
        beta = -atan(v_vel/u_vel);
    }

    // Yoshimura(2009) Low Speed Model
    int i_max = 100;
    double Y_SUM = 0.0;
    double N_SUM = 0.0;
    for (int i = 1; i <= i_max; ++i) {
        double x0 = -0.5 + (i-1)/i_max;
        double x1 = -0.5 + i/i_max;

        double comp0 = v_vel + CrY*r*Lpp*x0;
        double comp1 = v_vel + CrY*r*Lpp*x1;
        double comp2 = v_vel + CrN*r*Lpp*x0;
        double comp3 = v_vel + CrN*r*Lpp*x1;

        Y_SUM = Y_SUM + 0.5*(abs(comp0)*comp0 +    abs(comp1)*comp1)/i_max;
        N_SUM = N_SUM + 0.5*(abs(comp2)*comp2*x0 + abs(comp3)*comp3*x1)/i_max;
    }

    double YHN_ND = -CD*Y_SUM;
    double NHN_ND = -CD*N_SUM;

    double XH = 0.5*rho_w*Lpp     *d* ((X0F_ND + (X0A_ND-X0F_ND)*(abs(beta)/pi))*u_vel*U  + (Xvr_ND+my_ND)*v_vel*r*Lpp);
    double YH = 0.5*rho_w*Lpp     *d* (Yv_ND*v_vel*abs(u_vel) + (Yr_ND-mx_ND)*r*Lpp*u_vel + YHN_ND);
    double NH = 0.5*rho_w*Lpp*Lpp *d* (Nv_ND*v_vel*u_vel      + Nr_ND*r*Lpp*abs(u_vel)    + NHN_ND);

    tau_hull[0] = XH;
    tau_hull[1] = YH;
    tau_hull[2] = NH;

}

void StateODE::prop_EO(const double t, const double* x, const double* u, double* tau_prop)
{
    // get state and control variables
    double u_vel = x[3];
    double v_vel = x[4];
    double r     = x[5];

    //double n     = u[1];
	double n = 10;

    if (abs(u_vel) < 1e-8)
    {
        u_vel = 1e-8;
    }

    double U    = sqrt(u_vel*u_vel + v_vel*v_vel);
    double v_ND, r_ND;
    if (abs(U) < 1e-5) {
        v_ND = 1e-8;
        r_ND = 1e-8;
    } else {
        v_ND = v_vel/U;
        r_ND = r*Lpp/U;
    }
    
    double wP = wP0 - tau*abs(v_ND + xP_ND*r_ND) + CP_ND*pow((v_ND + xP_ND*r_ND),2);
    double Js;
    if (n == 0) {
        Js = 1.0e10;
    } else {
        Js = u_vel/(Dp*n);
    }

    double XP, YP, NP;
    if (n >= 0) {
        double J  = Js*(1-wP);
        double KT = a0 + a1*J + a2*J*J;

        XP = rho_w*pow(Dp,4)*(n*n)*(1-tP)*KT;
        YP = 0;
        NP = 0;
    } else {
        if (Js >= C10)
        {
            XP = rho_w*pow(Dp,4)*(n*n)*(C6+C7*Js);
        } else {
            XP = rho_w*pow(Dp,4)*(n*n)*C3;
        }
        
        if (Jsyn <= Js && Js <= Jsyn0)
        {
            YP = 0.5*rho_w*Lpp    *d*pow(n*Dp,2)*(A1+A2*Js);
            NP = 0.5*rho_w*Lpp*Lpp*d*pow(n*Dp,2)*(B1+B2*Js);

        } else if (Js < Jsyn) {
            YP = 0.5*rho_w*Lpp    *d*pow(n*Dp,2)*(A3+A4*Js);
            NP = 0.5*rho_w*Lpp*Lpp*d*pow(n*Dp,2)*(B3+B4*Js);

        } else {
            YP = 0.5*rho_w*Lpp    *d*pow(n*Dp,2)*A5;
            NP = 0.5*rho_w*Lpp*Lpp*d*pow(n*Dp,2)*B5;
        }
    }

    tau_prop[0] = XP;
    tau_prop[1] = YP;
    tau_prop[2] = NP;
}

void StateODE::rudd_EO(const double t, const double* x, const double* u, double* tau_rudd)
{
    // get state and control variables
    double u_vel = x[3];
    double v_vel = x[4];
    double r     = x[5];

    double delta = u[0];

	delta = remainder(delta + ((delta > 0) - (delta < 0)) * pi, 2 * pi) - ((delta > 0) - (delta < 0))* pi;

	//double delta = 0;
    //double n     = u[1];
	double n = 10;

    if (abs(u_vel) < 1e-8)
    {
        u_vel = 1e-8;
    }

    // Propeller related calculations
    double U    = sqrt(u_vel*u_vel + v_vel*v_vel);
    double v_ND, r_ND;
    if (abs(U) < 1e-5) {
        v_ND = 1e-8;
        r_ND = 1e-8;
    } else {
        v_ND = v_vel/U;
        r_ND = r*Lpp/U;
    }
    
    double wP = wP0 - tau*abs(v_ND + xP_ND*r_ND) + CP_ND*pow((v_ND + xP_ND*r_ND),2);
    double wPdash = 1-wP;    
    if (wPdash > 1) {// 2018/1/23 Bug fix (wP exceeded hundred or over)
        wPdash = 1;
    } else if (wPdash < 0) {
        wPdash = 0;
    }

    double uP = u_vel*wPdash;
    double Js;
    if (n == 0) {
        Js = 1.0e10;
    } else {
        Js = u_vel/(Dp*n);
    }

    // Rudder calculation
    double UR, aR;
    if (n > 0) {
        double J = Js*(1-wP);
        double KT = a0+a1*J+a2*J*J;
        double Ep = epsilon+kappa*(sqrt(1+8*KT/(pi*J*J))-1);
        double uR = uP*Ep;
        
        double vR;
        if (v_vel+xR*r >= 0) {// r added by Myo
            vR = gammaP*(v_vel+lR*r);
        } else {
            vR = gammaS*(v_vel+lR*r);
        }

        UR = sqrt(uR*uR+vR*vR);
        aR = delta + atan(vR/uR);
    } else {  // n < 0 (Reverse) 
        UR = sqrt(u_vel*u_vel+pow(-v_vel+xR*r,2));
        aR = delta + atan((-v_vel+xR*r)/abs(u_vel));
    }

    double FN;
    if (u_vel >= 0) {
        FN = 0.5*rho_w*AR*fa*UR*UR*sin(aR); 
    } else {
        FN = -0.5*rho_w*AR*fa*UR*UR*sin(aR);
    }

    double XR = -(1-tR)*    FN*sin(delta);
    double YR = -(1+aH)*    FN*cos(delta);
    double NR = -(xR+aH*xH)*FN*cos(delta);

    tau_rudd[0] = XR;
    tau_rudd[1] = YR;
    tau_rudd[2] = NR;
}