/*
Functions in this file take care of all the black-scholes pricing/implied volatility calculations.

@author: Justin Clark

This software is free and open to the public, but if used in a commercial setting, please cite this project.
*/

#include <math.h>
#include <string.h>

// directives that calculate a normal pdf or normal cdf of a given value
const double INV_SQRT_2PI = 1.0 / sqrt(2.0 * M_PI); // constant so we don't have to calculate each time
#define normCDF(d) (0.5 * erfc(-d * M_SQRT1_2))
#define normPDF(d) (INV_SQRT_2PI * exp(-0.5 * d * d))

// struct that contains all the calculated greeks
struct Greeks
{
    double delta;
    double gamma;
    double theta;
    double vega;
    double rho;
    Greeks() : delta(0), gamma(0), theta(0), vega(0), rho(0) {}
};

/////////////////// BLACK SCHOLES PRICING MODEL ////////////////////
double blsPrice(double S, double K, double T, double v, const char* right, double r=0.0, double d=0.0)
{
    double d1 = (log(S/K) + ((r-d)+(v*v/2.0)) * T) / (v*sqrt(T));
    double d2 = d1 - v*sqrt(T);
    
    if (strcmp(right,"C") == 0)
        return S * exp(-d*T) * normCDF(d1) - K * exp(-r*T) * normCDF(d2);
    else if (strcmp(right,"P") == 0)
        return K * exp(-r*T) * normCDF(-d2) - S * exp(-d*T) * normCDF(-d1);
    else {
        fprintf(stderr, "Right must be either call or put -- exiting\n");
        exit(-1);
    }
}

//////////////////// BLACK SCHOLES GREEKS MODEL /////////////////////
Greeks blsGreeks(double S, double K, double T, double v, const char* right, double r=0.0, double d=0.0)
{
    Greeks g;
    if (v == 0) return g;
    double d1  = (log(S/K) + ((r-d)+(v*v/2.0)) * T) / (v*sqrt(T));
    double d2  = d1 - v*sqrt(T);
    
    // set greeks depending on the right (gamma and vega are right-neutral)
    g.gamma    = normPDF(d1) / (S*v*sqrt(T));
    g.vega     = S*sqrt(T)*normPDF(d1);

    double t_0 = -(S*v*normPDF(d1))/(2*sqrt(T)); // same term of theta for both right of C and P
    if (strcmp(right,"C") == 0){
        g.delta = normCDF(d1);        
        g.theta = t_0 - r*K*exp(-r*T)*normCDF(d2);
        g.rho   = K*T*exp(-r*T)*normCDF(d2);
        
    } else if (strcmp(right,"P") == 0){
        g.delta = normCDF(-d1);
        g.theta = t_0 + r*K*exp(-r*T)*normCDF(-d2);
        g.rho   = -K*T*exp(-r*T)*normCDF(-d2);
        
    } else {
        fprintf(stderr, "Right must be either call or put -- exiting\n");
        exit(-1);
    }

    return g;
}

////////////////////// BLACK SCHOLES IMPLIED VOLATILITY //////////////////
double blsImpv(double S, double K, double T, double option_price, const char* right, double r=0.0, double d=0.0)
{
    double bm_iterations = 10; // number of iterations used by the bisection method
    double nm_iterations = 10; // maximum number of iterations used by newtons method
    double precision     = 1.0e-6;
    double v             = 0.0;

    // check no-arbitrage conditions
    if (strcmp(right,"C") == 0 and option_price < S-K) return v;
    else if (strcmp(right,"P") == 0 and option_price < K-S) return v;

    // get starting point with bisection method
    double a = 0.0001;
    double b = 100.0;
    double diff;
    double p;
    double vega;

    for (unsigned int i = 0; i < bm_iterations; i++){
        v    = sqrt(a*b);
        p    = blsPrice(S, K, T, v, right, r, d);
        diff = p - option_price;
        if      (diff < 0) a = v;
        else if (diff > 0) b = v;
    }
    
    // finish with newton's method until convergence
    for (unsigned int i = 0; i < nm_iterations; i++){
        p    = blsPrice(S, K, T, v, right, r, d);
        vega = blsGreeks(S, K, T, v, right, r, d).vega;
        diff = option_price - p;
        if (abs(diff) < precision) break;
        v    = v + (diff/vega);
        if (v >= 10 or v <= 0) v = 1.0; // checking for ridiculous values
    }

    return v;
}





