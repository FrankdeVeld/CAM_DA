#include <dace/dace.h>
#include <cmath>
#include <iomanip>
#include <typeinfo>
using namespace std;
using namespace DACE;

//#warning This example only works if the DACE was built with support for the (non-default) AlgebraicMatrix type!

double min( double a, double b ) {
    
    double res = a;
    if(a>b){res = b;}
    return res;
}


double max( double a, double b ) {
    
    double res = a;
    if(b>a){res = b;}
    return res;
}


double normtmp( int N, vector<double> X ) {
    
    int I;
    double res = 0.0;
    
    for (I=0; I<N; I++) {
        if (X[I]<0) { X[I] = -X[I]; }
        res = max(res,X[I]);
    }

    return res;
}

template <typename T> T checkatan2(T angle, T a, T b) {
    if ( DACE::cons(a) > 0 && DACE::cons(b) < 0 ) return angle + M_PI;
    else if ( DACE::cons(a) < 0 && DACE::cons(b) < 0 ) return angle - M_PI;
    else return angle;
}

template <typename T> T atan2_mod(T a, T b) {
    T angle = atan(a/b);
    checkatan2(angle, a, b);
    return angle;
}

double Initialtca(int Scenario,double MuEarth) {
    double tca; 
    switch(Scenario) {
        case 1:
            {
            tca = 7.68061e3;
            break;
            }
        case 2:
            {
            tca = 7.68061e3;
            break;
            }
        case 3:
            {
            tca = 7.68061e3;
            break;
            }
    }
    return tca;
}

AlgebraicVector<double> InitialXs(int Scenario,double MuEarth) {
    AlgebraicVector<double> xs_t0(6);
    switch(Scenario) {
        case 1:
            {
            // Set initial conditions
            const double eccs = 0.5;

            const double sec  = -20035.5;               
            
            // xs0[0] = sec    + DA(1); // sec km altitude
            // xs0[1] = 0.0    + DA(2);
            // xs0[2] = 0.0    + DA(3);
            // xs0[3] = 0.0    + DA(4);
            // xs0[4] = 0.0    + DA(5);  // at apogee
            // xs0[5] = -sqrt(mu/abs(sec))*sqrt(1-eccs) + DA(6);

            // Secondary initial position, guaranteeing a collision at tca
            xs_t0[0] = -20035.49963;     // sec km altitude
            xs_t0[1] = 0.0;          
            xs_t0[2] = -2.72058944;    
            xs_t0[3] = 0.0008565374959; 
            xs_t0[4] = 0.0;             
            xs_t0[5] = -3.153940884;      
            // Nominal miss distance: 500 m
            break;
            }
        case 2:
            {
            // Set initial conditions
            const double eccs = 0.5;

            const double sec  = 20035.5;               
            
            // xs0[0] = sec    + DA(1); // sec km altitude
            // xs0[1] = 0.0    + DA(2);
            // xs0[2] = 0.0    + DA(3);
            // xs0[3] = 0.0    + DA(4);
            // xs0[4] = 0.0    + DA(5);  // at apogee
            // xs0[5] = -sqrt(mu/abs(sec))*sqrt(1-eccs) + DA(6);

            // Secondary initial position, guaranteeing a collision at tca
            xs_t0[0] = 20035.49963;     // sec km altitude
            xs_t0[1] = 0.0;          
            xs_t0[2] = 2.72058944;    
            xs_t0[3] = -0.0008565374959; 
            xs_t0[4] = 0.0;             
            xs_t0[5] = 3.153940884;      
            // Nominal miss distance: 500 m
            break;
            }
        case 3:
            {
            // Set initial conditions
            const double eccs = 0.5;

            const double sec  = 20034;              

            // Secondary initial position, guaranteeing a collision at tca
            xs_t0[0] = sec;                              // prim km altitude
            xs_t0[1] = 0.0;                            
            xs_t0[2] = 0.0;                             
            xs_t0[3] = 0.0;                            
            xs_t0[4] = sqrt(MuEarth/abs(sec))*sqrt(1-eccs);  // at apogee
            xs_t0[5] = 0.0;
            break;      
            }
    }
    return xs_t0;
}

AlgebraicVector<double> InitialXp(int Scenario,double MuEarth) {
    AlgebraicVector<double> xp_t0(6);
    switch(Scenario) {
        case 1:
            {
            // Set initial conditions
            const double eccp = 0.5;

            const double prim = -20034;
            
            // Some random primary
            xp_t0[0] = prim;                              // prim km altitude
            xp_t0[1] = 0.0;                            
            xp_t0[2] = 0.0;                             
            xp_t0[3] = 0.0;                            
            xp_t0[4] = -sqrt(MuEarth/abs(prim))*sqrt(1-eccp);  // at apogee
            xp_t0[5] = 0.0;
            break;   
            }
        case 2:
            {
            // Set initial conditions
            const double eccp = 0.5;

            const double prim = 20034;
            
            // Some random primary
            xp_t0[0] = prim;                              // prim km altitude
            xp_t0[1] = 0.0;                            
            xp_t0[2] = 0.0;                             
            xp_t0[3] = 0.0;                            
            xp_t0[4] = sqrt(MuEarth/abs(prim))*sqrt(1-eccp);  // at apogee
            xp_t0[5] = 0.0;
            break;   
            }   
        case 3:
            {
            // Set initial conditions
            xp_t0[0] = 20035.49963;     // prim km altitude
            xp_t0[1] = 0.0;          
            xp_t0[2] = 2.72058944;    
            xp_t0[3] = -0.0008565374959; 
            xp_t0[4] = 0.0;             
            xp_t0[5] = 3.153940884;  
            break;   
            }                     

    }
    return xp_t0;
}


template<typename T, typename U>
DACE::AlgebraicVector<T> RK78(const int N, DACE::AlgebraicVector<T> Y0, DACE::AlgebraicVector<T> U0,
	const U X0, const U X1, AlgebraicVector<T> (*dyn)(AlgebraicVector<T>,AlgebraicVector<T>,double,double,double), const double mu, const double Lsc,  
    const bool returnIntermediatePoints = false, const double tolerance = 1.e-11){

    double ERREST;
    double H0_i;
    double HS_i;
    double H1_i;
    const double EPS = tolerance;
    const double BS = 20. * EPS;
	
    H0_i=1e-14;
    HS_i=1e+5;
    H1_i=1e+1;	
	
	const double H0=H0_i;
	const double HS=HS_i;
	const double H1=H1_i;
	
    T Z[N][16];

    DACE::AlgebraicVector<T> Yout = Y0;
    Yout.push_back(X0);
    if (returnIntermediatePoints)
    {
        Yout.reserve(128);
    }

    DACE::AlgebraicVector<T> Y1(N);
    std::vector<double> Y1cons(N);

    double VIHMAX = 0.0;
    U X, H;
    double RFNORM, HH0, HH1;

    const double HSQR = 1.0 / 9.0;
    double A[13], C[13], D[13];
    double B[13][12];



    A[0] = 0.0; A[1] = 1.0 / 18.0; A[2] = 1.0 / 12.0; A[3] = 1.0 / 8.0; A[4] = 5.0 / 16.0; A[5] = 3.0 / 8.0;
    A[6] = 59.0 / 400.0; A[7] = 93.0 / 200.0; A[8] = 5490023248.0 / 9719169821.0; A[9] = 13.0 / 20.0; A[10] = 1201146811.0 / 1299019798.0; A[11] = 1.0;
    A[12] = 1.0;

    B[0][0] = 0.0; B[0][1] = 0.0; B[0][2] = 0.0; B[0][3] = 0.0; B[0][4] = 0.0;
    B[0][5] = 0.0; B[0][6] = 0.0; B[0][7] = 0.0; B[0][8] = 0.0; B[0][9] = 0.0;
    B[0][10] = 0.0; B[0][11] = 0.0;

    B[1][0] = 1.0 / 18.0; B[1][1] = 0.0; B[1][2] = 0.0; B[1][3] = 0.0; B[1][4] = 0.0;
    B[1][5] = 0.0; B[1][6] = 0.0; B[1][7] = 0.0; B[1][8] = 0.0; B[1][9] = 0.0;
    B[1][10] = 0.0; B[1][11] = 0.0;

    B[2][0] = 1.0 / 48.0; B[2][1] = 1.0 / 16.0; B[2][2] = 0.0; B[2][3] = 0.0; B[2][4] = 0.0;
    B[2][5] = 0.0; B[2][6] = 0.0; B[2][7] = 0.0; B[2][8] = 0.0; B[2][9] = 0.0;
    B[2][10] = 0.0; B[2][11] = 0.0;

    B[3][0] = 1.0 / 32.0; B[3][1] = 0.0; B[3][2] = 3.0 / 32.0; B[3][3] = 0.0; B[3][4] = 0.0;
    B[3][5] = 0.0; B[3][6] = 0.0; B[3][7] = 0.0; B[3][8] = 0.0; B[3][9] = 0.0;
    B[3][10] = 0.0; B[3][11] = 0.0;

    B[4][0] = 5.0 / 16.0; B[4][1] = 0.0; B[4][2] = -75.0 / 64.0; B[4][3] = 75.0 / 64.0; B[4][4] = 0.0;
    B[4][5] = 0.0; B[4][6] = 0.0; B[4][7] = 0.0; B[4][8] = 0.0; B[4][9] = 0.0;
    B[4][10] = 0.0; B[4][11] = 0.0;

    B[5][0] = 3.0 / 80.0; B[5][1] = 0.0; B[5][2] = 0.0; B[5][3] = 3.0 / 16.0; B[5][4] = 3.0 / 20.0;
    B[5][5] = 0.0; B[5][6] = 0.0; B[5][7] = 0.0; B[5][8] = 0.0; B[5][9] = 0.0;
    B[5][10] = 0.0; B[5][11] = 0.0;

    B[6][0] = 29443841.0 / 614563906.0; B[6][1] = 0.0; B[6][2] = 0.0; B[6][3] = 77736538.0 / 692538347.0; B[6][4] = -28693883.0 / 1125000000.0;
    B[6][5] = 23124283.0 / 1800000000.0; B[6][6] = 0.0; B[6][7] = 0.0; B[6][8] = 0.0; B[6][9] = 0.0;
    B[6][10] = 0.0; B[6][11] = 0.0;

    B[7][0] = 16016141.0 / 946692911.0; B[7][1] = 0.0; B[7][2] = 0.0; B[7][3] = 61564180.0 / 158732637.0; B[7][4] = 22789713.0 / 633445777.0;
    B[7][5] = 545815736.0 / 2771057229.0; B[7][6] = -180193667.0 / 1043307555.0; B[7][7] = 0.0; B[7][8] = 0.0; B[7][9] = 0.0;
    B[7][10] = 0.0; B[7][11] = 0.0;

    B[8][0] = 39632708.0 / 573591083.0; B[8][1] = 0.0; B[8][2] = 0.0; B[8][3] = -433636366.0 / 683701615.0; B[8][4] = -421739975.0 / 2616292301.0;
    B[8][5] = 100302831.0 / 723423059.0; B[8][6] = 790204164.0 / 839813087.0; B[8][7] = 800635310.0 / 3783071287.0; B[8][8] = 0.0; B[8][9] = 0.0;
    B[8][10] = 0.0; B[8][11] = 0.0;

    B[9][0] = 246121993.0 / 1340847787.0; B[9][1] = 0.0; B[9][2] = 0.0; B[9][3] = -37695042795.0 / 15268766246.0; B[9][4] = -309121744.0 / 1061227803.0;
    B[9][5] = -12992083.0 / 490766935.0; B[9][6] = 6005943493.0 / 2108947869.0; B[9][7] = 393006217.0 / 1396673457.0; B[9][8] = 123872331.0 / 1001029789.0; B[9][9] = 0.0;
    B[9][10] = 0.0; B[9][11] = 0.0;

    B[10][0] = -1028468189.0 / 846180014.0; B[10][1] = 0.0; B[10][2] = 0.0; B[10][3] = 8478235783.0 / 508512852.0; B[10][4] = 1311729495.0 / 1432422823.0;
    B[10][5] = -10304129995.0 / 1701304382.0; B[10][6] = -48777925059.0 / 3047939560.0; B[10][7] = 15336726248.0 / 1032824649.0; B[10][8] = -45442868181.0 / 3398467696.0; B[10][9] = 3065993473.0 / 597172653.0;
    B[10][10] = 0.0; B[10][11] = 0.0;

    B[11][0] = 185892177.0 / 718116043.0; B[11][1] = 0.0; B[11][2] = 0.0; B[11][3] = -3185094517.0 / 667107341.0; B[11][4] = -477755414.0 / 1098053517.0;
    B[11][5] = -703635378.0 / 230739211.0; B[11][6] = 5731566787.0 / 1027545527.0; B[11][7] = 5232866602.0 / 850066563.0; B[11][8] = -4093664535.0 / 808688257.0; B[11][9] = 3962137247.0 / 1805957418.0;
    B[11][10] = 65686358.0 / 487910083.0; B[11][11] = 0.0;

    B[12][0] = 403863854.0 / 491063109.0; B[12][1] = 0.0; B[12][2] = 0.0; B[12][3] = -5068492393.0 / 434740067.0; B[12][4] = -411421997.0 / 543043805.0;
    B[12][5] = 652783627.0 / 914296604.0; B[12][6] = 11173962825.0 / 925320556.0; B[12][7] = -13158990841.0 / 6184727034.0; B[12][8] = 3936647629.0 / 1978049680.0; B[12][9] = -160528059.0 / 685178525.0;
    B[12][10] = 248638103.0 / 1413531060.0; B[12][11] = 0.0;

    C[0] = 14005451.0 / 335480064.0; C[1] = 0.0; C[2] = 0.0; C[3] = 0.0; C[4] = 0.0; C[5] = -59238493.0 / 1068277825.0;
    C[6] = 181606767.0 / 758867731.0; C[7] = 561292985.0 / 797845732.0; C[8] = -1041891430.0 / 1371343529.0; C[9] = 760417239.0 / 1151165299.0; C[10] = 118820643.0 / 751138087.0; C[11] = -528747749.0 / 2220607170.0;
    C[12] = 1.0 / 4.0;

    D[0] = 13451932.0 / 455176623.0; D[1] = 0.0; D[2] = 0.0; D[3] = 0.0; D[4] = 0.0; D[5] = -808719846.0 / 976000145.0;
    D[6] = 1757004468.0 / 5645159321.0; D[7] = 656045339.0 / 265891186.0; D[8] = -3867574721.0 / 1518517206.0; D[9] = 465885868.0 / 322736535.0; D[10] = 53011238.0 / 667516719.0; D[11] = 2.0 / 45.0;
    D[12] = 0.0;

    for (int i = 0; i < N; i++)
    {
        Z[i][0] = Y0[i];
        Z[i][1] = 0.0;
    }

    double sig=1.;
    if(X1<X0){
	sig=-1.;
    }
    H=sig*abs(HS);
    HH0 = abs(H0); HH1 = abs(H1);
    X = X0; RFNORM = 0.0; ERREST = 0.0;
    
    while (std::abs(DACE::cons(X)-DACE::cons(X1)) > 0.0)
	{
        // compute new stepsize
		if (RFNORM > 0)
        {
            H = H*std::min(4.0, exp(HSQR*log(EPS / RFNORM)));
        }
        if (abs(H)>abs(HH1))
        {
            H = sig*HH1;
        }
        else if (abs(H)<abs(HH0)*0.99)
        {
            H = sig*HH0;
            std::cout << "--- WARNING, MINIMUM STEPSIZE REACHED IN RK" << std::endl;
        }
        if ((DACE::cons(X) + DACE::cons(H) - DACE::cons(X1))*DACE::cons(H)>0)
        {
            H = X1 - X;
        }
		for (int j = 0; j<13; j++) {
            for (int i = 0; i<N; i++)
            {
                Y0[i] = 0.0; // EVALUATE RHS AT 13 POINTS

                for (int k = 0; k<j; k++)
                {
                    Y0[i] = Y0[i] + Z[i][k + 3] * B[j][k];
                }

                Y0[i] = H*Y0[i] + Z[i][0];
            }
            Y1 = dyn(Y0, U0, X + H*A[j], mu, Lsc);

            for (int i = 0; i<N; i++)
            {
                Z[i][j + 3] = Y1[i];
            }
        }
	
        for (int i = 0; i<N; i++) {

            Z[i][1] = 0.0; Z[i][2] = 0.0; // EXECUTE 7TH,8TH ORDER STEPS

            for (int j = 0; j<13; j++)
            {
                Z[i][1] = Z[i][1] + Z[i][j + 3] * D[j];
                Z[i][2] = Z[i][2] + Z[i][j + 3] * C[j];
            }

            Y1[i] = (Z[i][2] - Z[i][1])*H;
            Z[i][2] = Z[i][2] * H + Z[i][0];
        }



        for (int i = 0; i<N; i++)
        {
            Y1cons[i] = DACE::cons(Y1[i]);
        }
	
        RFNORM = normtmp(N, Y1cons); // ESTIMATE ERROR AND DECIDE ABOUT BACKSTEP

        if ((RFNORM>BS) && (abs(H / H0)>1.2))
        {
            H = H / 3.0;
            RFNORM = 0;
        }
        else
        {
            for (int i = 0; i<N; i++)
            {
                Z[i][0] = Z[i][2];
            }
            X = X + H;
            VIHMAX = std::max(VIHMAX, DACE::cons(H));
            ERREST = ERREST + RFNORM;

            if (returnIntermediatePoints)
            {
                for (int i = 0; i<N; i++)
                {
                    Yout.push_back(Z[i][0]);
                }
                Yout.push_back(X);
            }
        }
    }

    if (returnIntermediatePoints)
    {
        return Yout;
    }

    // Return final state and time
    for (int i = 0; i<N; i++)
    {
        Y1[i] = Z[i][0];
    }
   // Y1.push_back(X);

    return Y1;

}

template<typename T, typename U> DACE::AlgebraicVector<T> 
	KeplerProp(const DACE::AlgebraicVector<T>& rv, const U& t, const double mu){

    int ord =  DACE::DA::getMaxOrder();

    DACE::AlgebraicVector<T> rv_fin(6);
    DACE::AlgebraicVector<T> rr0(3), vv0(3);
    for (int i=0; i<3; i++) {
        rr0[i] = rv[i];
        vv0[i] = rv[i+3];
    }
    DACE::AlgebraicVector<T> hh = DACE::cross(rr0,vv0);
    T h = DACE::vnorm(hh);
    T r0 = DACE::vnorm(rr0);
    T v0 = DACE::vnorm(vv0);

    T a = mu/(2*mu/r0 -v0*v0);
    T p = h*h/mu;
    T sigma0 = DACE::dot(rr0,vv0)/std::sqrt(mu);

    double tol = 1.0;
    int iter = 0;
    double scl = 1e-4;

    T F, Ft, G, Gt;

    if (DACE::cons(a)>0) {

        T MmM0 = t * sqrt(mu/a/a/a);
        T EmE0 = DACE::cons(MmM0);

        while ( tol>1e-13 || iter < ord + 1) {
			iter++;
			T fx0 = -(MmM0) + (EmE0) + (sigma0)/sqrt((a))*
				(1 - cos((EmE0))) - (1-(r0)/(a)) * sin((EmE0));
			T fxp = 1 + (sigma0)/sqrt((a)) * sin((EmE0)) - 
				(1-(r0)/(a)) * cos((EmE0));
			tol = std::fabs(DACE::cons(fx0/fxp));
			EmE0 = EmE0 - fx0/fxp;
        }


        T theta = 2*atan2_mod(sqrt(a*p)*tan(EmE0/2), r0 + sigma0*sqrt(a)*tan(EmE0/2));
        T r = p*r0 / (r0 + (p-r0)*cos(theta) - sqrt(p)*sigma0*sin(theta));

        //{compute the Lagrangian coefficients}
        F = 1 - a/r0 * (1 - cos(EmE0)) ;
        G = a*sigma0/sqrt(mu)*(1 - cos(EmE0)) + r0 * sqrt(a/mu) * sin(EmE0);
        Ft = - sqrt(mu*a)/(r*r0) * sin(EmE0);
        Gt = 1 - a/r * (1-cos(EmE0));
    }
    else {
        std::cout << "Negative semimajor axis" << std::endl;
        T NmN0 = t*sqrt(-mu/a/a/a);
        T HmH0 = 0.0;


        while(tol>1e-14 || iter < ord + 1){
            iter ++;
            T fx0 = - (NmN0) - (HmH0) + (sigma0)/sqrt((-a)) * 
				(-1 + cosh((HmH0))) + (1-(r0)/(a)) * sinh((HmH0)) ;
            T fxp = -1 + (sigma0)/sqrt((-a))*sinh((HmH0)) + 
				(1-(r0)/(a))*cosh((HmH0));
            tol = std::abs(DACE::cons(fx0/fxp));
            HmH0 = HmH0 - fx0/fxp;
        }


        for( iter=0; iter<ord; iter++){
            T fx0 = - (NmN0) - HmH0 + (sigma0)/sqrt((-a)) * 
				(-1 + cosh(HmH0)) + (1-(r0)/(a)) * sinh(HmH0) ;
            T fxp = -1 + (sigma0)/sqrt((-a))*sinh(HmH0) + 
				(1-(r0)/(a))*cosh(HmH0);
            HmH0 = HmH0 - fx0/fxp;
        }

        //{DACE::DA expansion of HmH0 parameter}
        double Htemp, DE = 1.0;

        F = 1 - a/r0 * (1 - cosh(HmH0));
        G = a*sigma0/sqrt(mu)*(1 - cosh(HmH0)) + r0 * sqrt(-a/mu) * sinh(HmH0);

        DACE::AlgebraicVector<T> rv_temp(3);
        for (int i=0; i<3; i++) {
            rv_temp[i] = F * rr0[i] + G * vv0[i];
        }
        T r = DACE::vnorm(rv_temp);
        Ft = - sqrt(mu*(-a))/(r*r0) * sinh(HmH0);
        Gt = 1 - a/r*(1-cosh(HmH0));
    }

    for (int i=0 ; i<3; i++) {
        rv_fin[i] = F * rr0[i] + G * vv0[i];
        rv_fin[i+3] = Ft * rr0[i] + Gt * vv0[i];
    }
    return rv_fin;
}

template <typename T> AlgebraicMatrix<T> rtn2eci(const AlgebraicVector<T> & rv)
{
  AlgebraicVector<T> r(3), v(3), Rvec(3), Tvec(3), Nvec(3);
  AlgebraicMatrix<T> dcm(3,3);
  r[0] = rv[0];
  r[1] = rv[1];
  r[2] = rv[2];
  v[0] = rv[3];
  v[1] = rv[4];
  v[2] = rv[5];

  Rvec = r; 
  Rvec = Rvec/vnorm(Rvec);

  Nvec = DACE::cross(r,v);
  Nvec = Nvec/vnorm(Nvec);

  Tvec =  DACE::cross(Nvec, Rvec);         

  dcm.at(0,0) = Rvec[0];   dcm.at(0,1) = Tvec[0];  dcm.at(0,2) = Nvec[0];
  dcm.at(1,0) = Rvec[1];   dcm.at(1,1) = Tvec[1];  dcm.at(1,2) = Nvec[1];
  dcm.at(2,0) = Rvec[2];   dcm.at(2,1) = Tvec[2];  dcm.at(2,2) = Nvec[2];


  return dcm;
}

template<typename T> AlgebraicVector<T> TBAcc( AlgebraicVector<T> x, AlgebraicVector<T> uRTN, double t, double mu, double Lsc )
{
    
    AlgebraicVector<T> pos(3), res(6);
    
    pos[0] = x[0]; pos[1] = x[1]; pos[2] = x[2];
    
    T r = pos.vnorm();
    
    //const double mu = 398600; // km^3/s^2
    
    AlgebraicMatrix<double> R2E = rtn2eci(cons(x));
    AlgebraicVector<T>     uECI = R2E*uRTN;

    res[0] = x[3];
    res[1] = x[4];
    res[2] = x[5];
    
    res[3] = -mu*pos[0]/(r*r*r) + uECI[0];
    res[4] = -mu*pos[1]/(r*r*r) + uECI[1];
    res[5] = -mu*pos[2]/(r*r*r) + uECI[2];
    
    return res;
    
}

DA findTCA(const AlgebraicVector<DA> xrel, const int nvar){

  AlgebraicVector<DA> rr(3), vv(3), MAPD(nvar), MAPI(nvar), dx(nvar);

  //relative velocities
  rr[0] = xrel[0]; rr[1] = xrel[1]; rr[2] = xrel[2];
  vv[0] = xrel[3]; vv[1] = xrel[4]; vv[2] = xrel[5];

  // we want rr to be orthogonal to vv, so dot(rr,vv) = 0
  DA rvdot = dot(rr,vv);
  double rvdotcons = cons(rvdot);

  // MAPD = (dot(rr,vv), dxx) = f(dxx, dt)
  MAPD[0] = rvdot-rvdotcons;
  for (int i = 1; i < nvar ; i++) {
    MAPD[i] = DA(i);
  }
  // MAPI = (dxx, dt) = f(dot(rr,vv), dxx)
  MAPI = MAPD.invert();

  // we need to evaluate the map in -cons(dot(rr,vv)), dxx
  dx[0] = - rvdotcons + 0*DA(1);
  for (int i = 1; i < nvar ; i++){
      dx[i] = DA(i);}
  MAPD = MAPI.eval(dx);

  // dt is the last row of the MAPD
  DA tca = MAPD[nvar - 1];
  return tca;
}



DA tcaInversion(int ControlCases, AlgebraicVector<DA> xp_tn_DA,AlgebraicVector<double> xs_tn_Vec, DA tCA_tn, double tCA_Nom, double MuEarth, double Lsc){ 
    switch(ControlCases){
        case 1:
            AlgebraicVector<DA>     xp_tnp1_DA(6); 
            AlgebraicVector<double> xs_tnp1_Vec(6);
            xp_tnp1_DA  = RK78(6, xp_tn_DA,  {DA(7),DA(8),DA(9)}, 0.0, tCA_Nom*0.1,TBAcc,MuEarth,Lsc); // Propagation from t(N=9) to t(N=10) of DA object
            xs_tnp1_Vec = RK78(6, xs_tn_Vec, {0.0, 0.0, 0.0},     0.0, tCA_Nom*0.1,TBAcc,MuEarth,Lsc); // Propagation from t(N=9) to t(N=10) of vector object
            // Variable tCA at t(N=9) due to control between t(N=9) and t(N=10)
            tCA_tn = 0.0 + DA(10);
            AlgebraicVector<DA> xp_tCA_DA(6), xs_tCA_DA(6), xrel_tCA_DA(6);

            xp_tCA_DA        = KeplerProp(xp_tnp1_DA,  tCA_tn, MuEarth);
            xs_tCA_DA        = KeplerProp(xs_tnp1_Vec, tCA_tn, MuEarth);
            xrel_tCA_DA      = xp_tCA_DA  - xs_tCA_DA;
            DA tCA_tn        = findTCA(xrel_tCA_DA, 10);
            break;
        //case 2: TBD, what is written is outdated
            // AlgebraicVector<DA> xptN10(6), xstN10(6);
            // xptN10 = KeplerProp(6, xptN9DA, utN9, 0.0, pi*sqrt(ap*ap*ap/mu)*0.1 + DeltatCAtN9,TBAcc); // Propagation from t(N=9) to t(N=10) + delta tCA
            // xstN10 = KeplerProp(6, xsf, utN9, 0.0, DeltatCAtN9,TBAcc); // Propagation from t(N=10) to t(N=10) + delta tCA

            // Eliminate tCA as a DA variable, formulate Taylor polynomial in terms of dr, dv and du as DA variables
            // tCAtN10 = findTCA(XReltN10, 11)
    };
    return tCA_tn;
}