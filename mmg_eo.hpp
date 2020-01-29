#ifndef MMG_EO_H
#define MMG_EO_H
#define _USE_MATH_DEFINES

#include <cmath>

class StateODE{
private:

// Constant Parameters
    static constexpr double pi       = 3.14159265359;
    static constexpr double rho_w    = 1000/9.80665; // why divided by 9.81

    static constexpr double Lpp      = 3.0;
    static constexpr double B        = 0.48925;
    static constexpr double d        = 0.20114;

    static constexpr double m_ND     = 0.2453;
    static constexpr double mx_ND    = 0.0467;
    static constexpr double my_ND    = 0.2579;
    static constexpr double I_Jzz_ND = 0.0286;
    static constexpr double xG_ND    = 0.03119;

    static constexpr double m     = m_ND      * (0.5*rho_w*(Lpp*Lpp)*d);
    static constexpr double mx    = mx_ND     * (0.5*rho_w*(Lpp*Lpp)*d);
    static constexpr double my    = my_ND     * (0.5*rho_w*(Lpp*Lpp)*d);
    static constexpr double I_Jzz = I_Jzz_ND  * (0.5*rho_w*(Lpp*Lpp*Lpp*Lpp)*d);
    static constexpr double xG    = xG_ND     * Lpp;

    // Hull Reaction Forces Coefficients
    static constexpr double Xuu_ND = -0.02139;
    static constexpr double Xvr_ND =  0.4585;
    static constexpr double Yv_ND  = -0.37283;
    static constexpr double Yr_ND  =  0.1160;
    static constexpr double Nv_ND  = -0.14575;
    static constexpr double Nr_ND  = -0.04849;

    static constexpr double CD     = -0.0591*Lpp/d + 1.848;
    static constexpr double CrY    =  0.5200*Lpp/B - 1.062;
    static constexpr double CrN    =  0.0742*Lpp/d - 0.297;

    static constexpr double X0F_ND = Xuu_ND;
    static constexpr double X0A_ND = -0.03189;

    // Propeller related Parameters
    static constexpr double Dp    =  0.084;
    static constexpr double tP    =  0.220;
    static constexpr double wP0   =  0.614;
    static constexpr double tau   =  0.871;
    static constexpr double CP_ND = -0.359;
    static constexpr double xP_ND = -0.517;

    static constexpr double a0    =  0.3278;
    static constexpr double a1    = -0.3223;
    static constexpr double a2    = -0.1560;

    static constexpr double A1    = -7.90e-5;
    static constexpr double A2    =  7.99e-3;
    static constexpr double A3    = -4.93e-3;
    static constexpr double A4    = -5.87e-3;
    static constexpr double A5    = -5.58e-4;
  
    static constexpr double B1    =  3.50e-5;
    static constexpr double B2    = -3.17e-3;
    static constexpr double B3    =  1.96e-3;
    static constexpr double B4    =  2.33e-3;
    static constexpr double B5    =  2.25e-4;

    static constexpr double C3    = -0.251;
    static constexpr double C6    = -0.175;
    static constexpr double C7    =  0.33;
    static constexpr double C10   = -0.233;

    static constexpr double Jsyn  = -0.35;
    static constexpr double Jsyn0 = -0.06;

    // Rudder related parameters
    static constexpr double AR      = 0.01063;
    static constexpr double lambda  = 1.539;
    static constexpr double tR      = 0.19;
    static constexpr double aH      = 0.393;
    static constexpr double xR      = -1.5;
    static constexpr double xH_ND   = -0.45;
    static constexpr double kappa   = 0.288;
    static constexpr double epsilon = 1.42;
    static constexpr double lR_ND   = -1.08;
    static constexpr double gammaS  = 0.4406;
    static constexpr double gammaP  = 0.3506;

    static constexpr double xH = xH_ND*Lpp;
    static constexpr double lR = lR_ND*Lpp;
    static constexpr double fa = 6.13/(2.25+lambda);

public:
    void MMG_EO_ODE(const double t, const double* x, const double* u, double* f);
    void hull_EO(const double t, const double* x, const double* u, double* tau_hull);
    void prop_EO(const double t, const double* x, const double* u, double* tau_prop);
    void rudd_EO(const double t, const double* x, const double* u, double* tau_rudd);
};

#endif // !MMG_EO_H