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
    //////////////////////////////////////////////////////////////// START OF INITIALISATION ///////////////////////////////////////////////////////////////////////////////
    string SaveName = "N100_Julia_Val_TTS";
  
    int i;
    int j;
    int k;
    DA::init( 2, 10 );                                                   // DA(1:6) = xp(tn), DA(7:9) = u(tn),  DA(10) = tCA(tn) 
    double Mu      = 398600;                                            // Gravitational parameter of the Earth     [km^3/s^2]
    double Lsc     = 6378; // km                                                 // Length scale [currently unused]
    double Vsc     = sqrt(Mu/Lsc);
    double Asc     = Vsc*Vsc/Lsc;
    double tsc     = Lsc/Vsc; 
    std::cout.precision(16);                                             // This is to print variables very precisely
    AlgebraicVector<double> xp_t0(6), xs_t0(6);                          // Initialise primary and seconday object state (Cartesian) at t0, which have a close approach at some time tCA
    double tCA_Nom;                                                      // Initialise the time of closest approach for these objects
    AlgebraicVector<double> u_Nom(3);                                    // Initialise nominal thrust (has not been tested for nonzero thrust)
    for (i=0; i<3; i++){                                                 // Nominal thrust is zero
        u_Nom[i]   = 0.0/Asc;
    }
    
    AlgebraicVector<double> xp_tf(6), xs_tf(6);                         // Initialise primary and seconday object states at tCA
    
    double MuEarth = 1;

    double ThrustMagnitude = 1e-7/Asc;                                      // Thrust magntiude                         [km/s^3]
    int Scenario = 2;                                                   // Initial condition scenarios
    // TODO: automate (ideally) -> automated tCA computation, or propagate backwards from tCA
    // Create larger scenario database

    int SaveMode = 1;
    // 1: Only save initial state, time to safety

    int NumPoints = 360;

    AlgebraicVector<double> tf_save_vec(NumPoints);
    AlgebraicMatrix<double> xpf_save_mat(NumPoints,6), xsf_save_mat(NumPoints,6);

    for(k=0;k<NumPoints; k++){
        xp_t0    = InitialXp(Scenario, Mu, Lsc, k, NumPoints);                            // Obtain initial state primary
        xs_t0    = InitialXs(Scenario, Mu, Lsc);                            // Obtain initial state secondary

        tCA_Nom  = Initialtca(Scenario, Mu)/tsc;                           // Obtain initial tCA
        
        for(i=0;i<3;i++){
            xp_t0[i] = xp_t0[i]/Lsc; 
            xs_t0[i] = xs_t0[i]/Lsc; 
            xp_t0[i+3] = xp_t0[i+3]/Vsc; 
            xs_t0[i+3] = xs_t0[i+3]/Vsc; 
        }
        // Everything now dimensionless

        xp_tf = RK78(6, xp_t0, u_Nom, 0.0, tCA_Nom,TBAcc,MuEarth,Lsc);                  // Propagation till nominal tca (primary)
        xs_tf = RK78(6, xs_t0, {0.0, 0.0, 0.0}, 0.0, tCA_Nom,TBAcc,MuEarth,Lsc);      // Propagation till nominal tca (secondary)
        int N;                                                              // Initialise number of thrust arcs in the method 
        N = 100;                                                             // Number should have a relation with the eccentricity of the orbit
        // TODO: Connect N with RTN profile of control; make adaptive

        int DM_Case;                                                                                                  // Choose which distance metric to use for quantifying risk and determining control
        DM_Case = 1;    
        // Case 1: Distance Metric is Euclidean Distance. 
        // Case 2: Distance Metric is Square Mahalanobis Distance (SMD)
        // Case 3: Distance Metric is PoC formula (Serra)

        int tCAHandling;                                                             // Strategy for handling a variable tCA
        tCAHandling = 2; 
        // Case 1: Approximation: no control between tcaNom and tca(u). Effectively: no change in tca
        // As such, a RK78 propagation of DA variables is performed till tcaNom, after which a Kepler propagation 
        // is done around tcaNom with tf a DA variable. Polynomial inversion is done to obtain tca as a 
        // Taylor polynomial of dr, dv and u of the previous segment

        // Case 2: The propagation is performed with a propagator like RK78 to a variable time tca as a DA variable (Picard-Lindelöf)
        // In first approximation, this boils down to using an Euler scheme

        // For iterative loop
        double t_n;
        double StepSizeN = static_cast<double>(1)/N;

        DA DM_tn, tCA_tn; 
        AlgebraicVector<DA> DeltaRB_tn(3);

        AlgebraicVector<double> xp_tn_Vec(6), xs_tn_Vec(6);                          // Initialise primary, secondary state at t(n) as vector
        AlgebraicVector<DA> u_tn(3);                                                 // Initialise control between t(n) and t(n+1) as DA object
        AlgebraicVector<DA> xp_tn_DA(6);                                             // Initialise primary state at t(n) as DA object 
        AlgebraicVector<DA> rp_tn_DA(3), vp_tn_DA(3);                                // Initialise primary position, velocity at t(n) as DA object   

        AlgebraicVector<DA>     xp_tnp1_DA(6);                                       // Initialise primary state at t(n+1) as DA object    
        AlgebraicVector<double> xs_tnp1_Vec(6);

        AlgebraicVector<DA> xs_tCA_DA(6), xp_tCA_DA(6);    

        AlgebraicVector<DA> Evaluated_tCA(10);                                                           

        double  R; // Hard-body radius
        AlgebraicMatrix<double> P(3,3); // Covariance matrix
        P.at(0,0)      = 1;
        P.at(0,1)      = 0.05; 
        P.at(0,2)      = 0.04; 
        P.at(1,0)      = P.at(0,1);
        P.at(1,1)      = 0.9; 
        P.at(1,2)      = 0.045; 
        P.at(2,0)      = P.at(0,2);
        P.at(2,1)      = P.at(1,2); 
        P.at(2,2)      = 0.75; 
        P = P*1e-6; // Small values
        // Symmetric covariance matrix

        R = 0.002; // Hard-body radius

        double tf_save;

        DA DM_NextIt, tCA_NextIt;
        AlgebraicVector<DA> DeltaRB_NextIt;

        double DM_Evaluated_Control;
        DM_Evaluated_Control = InitialDM(Scenario, Lsc);
        double DM_Safety;
        DM_Safety = 1.0/Lsc/Lsc;
        double tCA_Evaluated_Control;
        AlgebraicVector<double> xp_tnp1_Evaluated_Control(6), u_OptFO_tn(3), DeltaRB_Evaluated_Control(3);
        /////////////////////////////////////////////////////////////// END OF INITIALISATION   /////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////// START OF ITERATIVE LOOP /////////////////////////////////////////////////////////////////////////////
        //Discretisation: first consider N timesteps, look at state at t(N-1), thrust arc from t(N-1) to tCA(t(N-1)), see if distance metric is sufficient or not
        int n= N-1;
        while(DM_Evaluated_Control < DM_Safety) {
            t_n = static_cast<double>(n)*1/N*tCA_Nom;                             // Compute the dt for the current thrust arc, based on the chosen value of N and tCA_Nom
            // Obtain the state of the primary, secondary as DA objects at t(n)
            // High-fidelity RK78 propagator
            xp_tn_Vec = RK78(6, xp_t0, u_Nom, 0.0, t_n,TBAcc,MuEarth,Lsc);               // Propagation from t0 to tn (primary)
            xs_tn_Vec = RK78(6, xs_t0, u_Nom, 0.0, t_n,TBAcc,MuEarth,Lsc);               // Propagation from t0 to tn (secondary)
            // You could do this outside of the loop: tabulate it for all time steps
            
            for (i=0; i<3; i++){
                rp_tn_DA[i] = xp_tn_Vec[i]   + DA(i+1);                                  // Make the primary position a DA object with (3) independent DA variables
                vp_tn_DA[i] = xp_tn_Vec[i+3] + DA(i+4);                                  // Make the primary velocity a DA object with (3) independent DA variables
                u_tn[i]     = DA(i+7);                                                   // Make the control acting on the primary a DA object with (3) independent DA variables
            }
            // Propagate rptn, vptn to tf + delta tCA
            
            for (i=0; i<3; i++){
                xp_tn_DA[i]   = rp_tn_DA[i];                                             // Obtain values for primary state from the position of the primary as a DA object
                xp_tn_DA[i+3] = vp_tn_DA[i];                                             // Obtain values for primary state from the velocity of the primary as a DA object
            }
                                                        
            // States at t(n+1)
            xp_tnp1_DA  = RK78(6, xp_tn_DA,  {u_Nom[0] + u_tn[0], u_Nom[1] + u_tn[1], u_Nom[2] + u_tn[2]}, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation from t(n) to t(n+1) of DA object primary
            xs_tnp1_Vec = RK78(6, xs_tn_Vec, {0.0, 0.0, 0.0},                                              0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation from t(n) to t(n+1) of vector object secondary
            
            if (n==N-1){  
                tCA_tn = DA(10);                                                                                                // Initialise tCA as DA object with (1) independent DA variable             
                std::tie(tCA_tn, xp_tCA_DA, xs_tCA_DA) = tcaInversion(tCAHandling, u_Nom, u_tn, xp_tnp1_DA, xs_tnp1_Vec, tCA_tn, MuEarth, Lsc);               // Resolve the dependency of tCA as a DA variable through polynomial inversion 
                // Iterative loop and finding control
                for (i=0;i<6;i++)
                {
                    Evaluated_tCA[i] = DA(i+1);
                }
                for (i=0;i<3;i++)
                {
                    Evaluated_tCA[i+6] = DA(i+7); // *ThrustMagnitude?
                }
                Evaluated_tCA[9] = tCA_tn;
                xp_tCA_DA = xp_tCA_DA.eval(Evaluated_tCA);
                xs_tCA_DA = xs_tCA_DA.eval(Evaluated_tCA);
                std::tie(xp_tnp1_Evaluated_Control, u_OptFO_tn, DeltaRB_Evaluated_Control, DM_Evaluated_Control, tCA_Evaluated_Control, DM_NextIt, tCA_NextIt, DeltaRB_NextIt) =  IterativeDA(n, N, DM_Case, xp_tCA_DA, xs_tCA_DA, ThrustMagnitude, DM_tn, tCA_tn, DeltaRB_tn, P, R); 
            } else {
                std::tie(xp_tnp1_Evaluated_Control, u_OptFO_tn, DeltaRB_Evaluated_Control, DM_Evaluated_Control, tCA_Evaluated_Control, DM_NextIt, tCA_NextIt, DeltaRB_NextIt) =  IterativeDA(n, N, DM_Case, xp_tnp1_DA, xs_tCA_DA, ThrustMagnitude, DM_tn, tCA_tn, DeltaRB_tn, P, R);
            }  
            DM_tn           = DM_NextIt;
            tCA_tn          = tCA_NextIt;
            DeltaRB_tn      = DeltaRB_NextIt;


            switch(SaveMode){
                case 1:
                {
                    // Nothing
                    break;
                } 
                case 2:
                {
                    // Full save - TBD
                    break;
                }
            }
            
            n--;
            if(n==0){
                cout << "Not enough time to avoid collision!" << endl;
                break;
            }
        }
        tf_save = t_n;
        switch(SaveMode){
            case 1:
            {
                tf_save_vec[k] = tf_save*tsc;
                for(i=0;i<3;i++){
                    xpf_save_mat.at(k,i)   = xp_tf[i]*Lsc;
                    xpf_save_mat.at(k,i+3) = xp_tf[i+3]*Vsc;
                    xsf_save_mat.at(k,i)   = xs_tf[i]*Lsc;
                    xsf_save_mat.at(k,i+3) = xs_tf[i+3]*Vsc;             
                }
                break;
            } 
            case 2:
            {
                // Full save - TBD
                break;
            }
        }
        // End of loop
        if (k%10==0){
            cout << static_cast<double>(k)/NumPoints*100 <<" \% done" << endl;
        }
        
    }
    switch(SaveMode){
        case 1:
        {
            FilePrintTTSSimple(SaveName, NumPoints, xpf_save_mat, xsf_save_mat, tf_save_vec);
            break;
        } 
        case 2:
        {
            // Full save - TBD
            break;
        }
    }
}
    


