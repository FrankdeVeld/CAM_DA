#include <dace/dace.h>
#include "DABackwardFunctionsIterative.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <typeinfo>
using namespace std;
using namespace DACE;

int main( void )
{   
    // Specifically for case N=100;
    string SaveName = "N50_Euc";
    int N=50;
    int i;
    int j;

    double ThrustMagnitude = 1e-7;                                      // Thrust magnitude                         [km/s^3]
    AlgebraicMatrix<double> xp_Nom(N,6); 
    AlgebraicMatrix<double> u_Opt(N,3); 
    AlgebraicVector<double> tCA_Vec(N,3); 

    int Scenario = 1;
    double MuEarth = 398600;
    double Lsc  = 1;
    double tCA_Nom;                                                      // Initialise the time of closest approach for these objects
    tCA_Nom  = Initialtca(Scenario, MuEarth);                            // Obtain initial tCA

    // Read the xp_Nom matrix
    std::ifstream xp("./write_read/xp_" + SaveName + ".dat");
    if (!xp.is_open()) {
        std::cerr << "Failed to open xp file" << std::endl;
        return 1;
    }
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (!(xp >> xp_Nom.at(i, j))) {
                std::cerr << "Error reading xp_Nom matrix at (" << i << ", " << j << ")" << std::endl;
                return 1;
            }
        }
    }
    xp.close();

    // Read the u_Opt matrix
    std::ifstream u("./write_read/u_" + SaveName + ".dat");
    if (!u.is_open()) {
        std::cerr << "Failed to open u file" << std::endl;
        return 1;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (!(u >> u_Opt.at(i, j))) {
                std::cerr << "Error reading u_Opt matrix at (" << i << ", " << j << ")" << std::endl;
                return 1;
            }
        }
    }
    u.close();

    // Read the tCA vector
    std::ifstream tCA("./write_read/tCA_" + SaveName + ".dat");
    if (!tCA.is_open()) {
        std::cerr << "Failed to open tCA file" << std::endl;
        return 1;
    }

    for (int i = 0; i < N; ++i) {
        if (!(tCA >> tCA_Vec[i])) {
            std::cerr << "Error reading tCA vector at (" << i << ")" << std::endl;
            return 1;
        }
    }
    tCA.close();

    double dtcaf = tCA_Vec[0];
    double dtca0 = tCA_Vec[N-1];

    AlgebraicVector<double> xs_t0(6);
    AlgebraicVector<double> xs_tf_new(6);
    xs_t0    = InitialXs(Scenario, MuEarth, Lsc);
    xs_tf_new = RK78(6, xs_t0, {0.0, 0.0, 0.0}, 0.0, tCA_Nom + dtcaf,TBAcc,MuEarth,Lsc);      // Propagation till nominal tca (secondary)

    AlgebraicMatrix<double> xnp1_save(N+1,6);
    AlgebraicMatrix<double> xnfull_save(N+1,6);

    for(j=0;j<6;j++)
    {
        xnp1_save.at(0,j)    = xp_Nom.at(0,j);
    }
    for(i=0; i<N; i++){ // xnadj error
        AlgebraicVector<double> xn(6); 
        AlgebraicVector<double> xnp1(6); 
        for(j=0; j<6; j++){
            xn[j] = xp_Nom.at(i,j);
        }
        AlgebraicVector<double> un(3); 
        for(j=0; j<3; j++){
            un[j] = u_Opt.at(i,j);
        }
        double StepSizeN = static_cast<double>(1)/N;
        if (i==N-1)
        {
            xnp1 = RK78(6, xn, un*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN+dtca0,TBAcc,MuEarth,Lsc); 
        } else {
            xnp1 = RK78(6, xn, un*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); 
        }
         
        for(j=0;j<6;j++)
        {
            xnp1_save.at(i+1,j)    = xnp1[j];
        } 
    }
    
    AlgebraicVector<double> xn(6); 
    AlgebraicVector<double> un(3);

    AlgebraicVector<double> x1(6);
    AlgebraicVector<double> xnext(6);

    for(j=0;j<6;j++)
    {
        xnfull_save.at(0,j)    = xp_Nom.at(0,j);
    }
    
    for(i=0;  i<N; i++){ // Full propagation
        double StepSizeN = static_cast<double>(1)/N;
        
        if (i==0){
            for(j=0; j<6; j++){
                xn[j] = xp_Nom.at(0,j);
            }
            for(j=0; j<3; j++){
                un[j] = u_Opt.at(0,j);
            }
            x1 = RK78(6, xn, un*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc);  
            xnext = x1;
        } else if(i==N-1)
        {
            for(j=0; j<3; j++){
                un[j] = u_Opt.at(i-1,j);
            }
            xnext = RK78(6, xnext, un*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN+dtcaf,TBAcc,MuEarth,Lsc); 
        }
        else {
            for(j=0; j<3; j++){
                un[j] = u_Opt.at(i-1,j);
            }
            xnext = RK78(6, xnext, un*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); 
        }
        for(j=0;j<6;j++)
        {
            xnfull_save.at(i+1,j)    = xnext[j];
        } 
        
    }

    AlgebraicVector<double> u_Val_R(3);                                 // Validation thrust: pure radial
    AlgebraicVector<double> u_Val_T(3);                                 // Validation thrust: pure tangential

    AlgebraicVector<double> un_T(3);
    AlgebraicVector<double> un_R(3);

    AlgebraicVector<double> x1_T(6);
    AlgebraicVector<double> xnext_T(6);
    AlgebraicVector<double> x1_R(6);
    AlgebraicVector<double> xnext_R(6);

    AlgebraicMatrix<double> xR_save(N+1,6);
    AlgebraicMatrix<double> xT_save(N+1,6);

    for (i=0; i<3; i++){                                                // Initialise pure radial, tangential thrust (zero values)
        u_Val_R[i]   = 0;
        u_Val_T[i]   = 0;
    }
    u_Val_R[0] =  1;                            // Initialise pure radial, tangential thrust (non-zero values)
    u_Val_T[1] =  1;
    for(j=0;j<6;j++)
    {
        xR_save.at(0,j)    = xp_Nom.at(0,j);
        xT_save.at(0,j)    = xp_Nom.at(0,j);
    }

    for(i=0;  i<N; i++){ // Pure R, T thrust validation
        double StepSizeN = static_cast<double>(1)/N;
        if (i==0){
            for(j=0; j<6; j++){
                xn[j] = xp_Nom.at(0,j);
            }
            for(j=0; j<3; j++){
                un_T[j] = u_Val_T[j];
                un_R[j] = u_Val_R[j];
            }
            x1_R = RK78(6, xn, un_R*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc);  
            x1_T = RK78(6, xn, un_T*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); 
            xnext_R = x1_R;
            xnext_T = x1_T;
        } else if(i==N-1) {
            for(j=0; j<3; j++){
                un_T[j] = u_Val_T[j];
                un_R[j] = u_Val_R[j];
            }
            xnext_R = RK78(6, xnext_R, un_R*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN+dtcaf,TBAcc,MuEarth,Lsc); 
            xnext_T = RK78(6, xnext_T, un_T*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN+dtcaf,TBAcc,MuEarth,Lsc);  
        } else {
            for(j=0; j<3; j++){
                un_T[j] = u_Val_T[j];
                un_R[j] = u_Val_R[j];
            }
            xnext_R = RK78(6, xnext_R, un_R*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); 
            xnext_T = RK78(6, xnext_T, un_T*ThrustMagnitude, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc);  
        }
        for(j=0;j<6;j++)
        {
            xR_save.at(i+1,j)    = xnext_R[j];
            xT_save.at(i+1,j)    = xnext_T[j];
        } 
    }


    // Output xnp1
    ofstream xnp1_output, xfull_output, xR_output, xT_output, xsfnew;
    string xnp1FileName = "./write_read/xnp1_Val_" + SaveName + ".dat";
    xnp1_output.open(xnp1FileName);
    xnp1_output << setprecision(16);
    for (i=0; i<N+1; i++)
    {
        for(j=0; j<6;j++)
        {
            xnp1_output << xnp1_save.at(i,j) << " ";
        }
        xnp1_output << endl;
    }
    xnp1_output.close();

    string xsfnewFileName = "./write_read/xsfnew_" + SaveName + ".dat";
    xsfnew.open(xsfnewFileName);
    xsfnew << setprecision(16);
    for(j=0; j<6;j++)
    {
        xsfnew << xs_tf_new[j] << " ";
    }
    xsfnew << endl;
    xsfnew.close();

    // Output xf
    string xfFileName = "./write_read/xfull_Val_" + SaveName + ".dat";
    xfull_output.open(xfFileName);
    xfull_output << setprecision(16);
    for (i=0; i<N+1; i++)
    {
        for(j=0; j<6;j++)
        {
            xfull_output << xnfull_save.at(i,j) << " ";
        }
        xfull_output << endl;
    }
    xfull_output.close();

    // Output R val
    string xRFileName = "./write_read/xfull_ValR_" + SaveName + ".dat";
    xR_output.open(xRFileName);
    xR_output << setprecision(16);
    for (i=0; i<N+1; i++)
    {
        for(j=0; j<6;j++)
        {
            xR_output << xR_save.at(i,j) << " ";
        }
        xR_output << endl;
    }
    xR_output.close();

    // Output T val
    string xTFileName = "./write_read/xfull_ValT_" + SaveName + ".dat";
    xT_output.open(xTFileName);
    xT_output << setprecision(16);
    for (i=0; i<N+1; i++)
    {
        for(j=0; j<6;j++)
        {
            xT_output << xT_save.at(i,j) << " ";
        }
        xT_output << endl;
    }
    xT_output.close();
}





