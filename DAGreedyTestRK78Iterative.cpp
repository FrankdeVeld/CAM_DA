#include <dace/dace.h>
#include "DABackwardFunctionsIterative.h"
#include <cmath>
#include <iomanip>
#include <typeinfo>
using namespace std;
using namespace DACE;

int main( void )
{   
    // DebugMode explicitly prints out some variables to check if something is wrong
    // 1 means on
    int DebugMode = 1;     
    int i;
    DA::init( 1, 10 );                                                   // DA(1:6) = xp(tn), DA(7:9) = u(tn),  DA(10) = tCA(tn) 
    std::cout.precision(10);                                             // This is to print variables very precisely
    AlgebraicVector<double> xp_t0(6), xs_t0(6);                          // Initialise primary and seconday object state (Cartesian) at t0, which have a close approach at some time tCA
    double tCA_Nom;                                                      // Initialise the time of closest approach for these objects
    AlgebraicVector<double> u_Nom(3);                                    // Initialise nominal thrust (has not been tested for nonzero thrust)
    for (i=0; i<3; i++){                                                 // Nominal thrust is zero
        u_Nom[i]   = 0.0;
    }
    
    
    AlgebraicVector<double> xp_tf(6), xs_tf(6);                         // Initialise primary and seconday object states at tCA
    
    double MuEarth = 398600;                                            // Gravitational parameter of the Earth     [km^3/s^2]
    double Lsc     = 1;                         -                       // Length scale [currently unused]
    
    double ThrustMagnitude = 1e-7;                                      // Thrust magntiude                         [km/s^3]
    AlgebraicVector<double> u_Val_R(3);                                 // Validation thrust: pure radial
    AlgebraicVector<double> u_Val_T(3);                                 // Validation thrust: pure tangential

    for (i=0; i<3; i++){                                                // Initialise pure radial, tangential thrust (zero values)
        u_Val_R[i]   = 0.0;
        u_Val_T[i]   = 0.0;
    }
    u_Val_R[0] = ThrustMagnitude;                                       // Initialise pure radial, tangential thrust (non-zero values)
    u_Val_T[1] = ThrustMagnitude;
    int Scenario = 3;                                                   // Initial condition scenarios
    // TODO: automate (ideally) -> automated tCA computation, or propagate backwards from tCA
    // Create larger scenario database

    xp_t0    = InitialXp(Scenario, MuEarth);                            // Obtain initial state primary
    xs_t0    = InitialXs(Scenario, MuEarth);                            // Obtain initial state secondary
    tCA_Nom  = Initialtca(Scenario, MuEarth);                           // Obtain initial tCA

    xp_tf = RK78(6, xp_t0, u_Nom, 0.0, tCA_Nom,TBAcc,MuEarth,Lsc);      // Propagation till nominal tca (primary)
    xs_tf = RK78(6, xs_t0, u_Nom, 0.0, tCA_Nom,TBAcc,MuEarth,Lsc);      // Propagation till nominal tca (secondary)

    int N;                                                              // Initialise number of thrust arcs in the method 
    N = 10;                                                             // Number should have a relation with the eccentricity of the orbit
    // TODO: Connect N with RTN profile of control; make adaptive

    DA EucDis;                                                          // Initialise DA variable Euclidean distance
    DA EucDis_CurrentIt;                                                // Initialise DA variable Euclidean distance (current iteration)
    DA EucDis_PreviousIt;                                               // Initialise DA variable Euclidean distance (previous iteration)
    DA SMD;                                                             // Initialise DA variable Square Mahalanobis Distance

    /////// Discretisation: first consider N timesteps, look at state at t(N-1), thrust arc from t(N-1) to tCA(t(N-1)), see if distance metric is sufficient or not
    for(int n= N-1; n>0; n--) {
        // Initialise 'main' DA object, as well as various parameters 'becoming' a DA object as a result
        AlgebraicVector<DA> rp_tn_DA(3);                                // Initialise primary position at t(n) as DA object
        AlgebraicVector<DA> vp_tn_DA(3);                                // Initialise primary velocity at t(n) as DA object
        AlgebraicVector<DA> u_tn(3);                                    // Initialise control between t(n) and t(n+1) as DA object
        
        
        AlgebraicVector<double> xp_tn_Vec(6), xs_tn_Vec(6);
        double t_n = static_cast<double>(n)*1/N*tCA_Nom;

        // Obtain the state of the primary, secondary as DA objects at t(n)
        // High-fidelity RK78 propagator
        xp_tn_Vec = RK78(6, xp_t0, u_Nom, 0.0, t_n,TBAcc,MuEarth,Lsc); // Propagation from t(N=0) to t(N=9)
        xs_tn_Vec = RK78(6, xs_t0, u_Nom, 0.0, t_n,TBAcc,MuEarth,Lsc);
        // You could do this outside of the loop

        // Make rp, vp, u @ t(N=n) variable DA objects 
        for (i=0; i<3; i++){
            rp_tn_DA[i] = xp_tn_Vec[i]   + DA(i+1);
            vp_tn_DA[i] = xp_tn_Vec[i+3] + DA(i+4);
            u_tn[i]  = DA(i+7);  
        }

        // Propagate rptn, vptn to tf + delta tCA
        AlgebraicVector<DA> xp_tn_DA(6);
        for (i=0; i<3; i++){
            xp_tn_DA[i]   = rp_tn_DA[i];
            xp_tn_DA[i+3] = vp_tn_DA[i];
        }


        int ControlApprox;
        ControlApprox = 1; 
        // Case 1: Approximation: no control between tcaNom and tca(u). Effectively: no change in tca
        // As such, a RK78 propagation of DA variables is performed till tcaNom, after which a Kepler propagation 
        // is done around tcaNom with tf a DA variable. Polynomial inversion is done to obtain tca as a 
        // Taylor polynomial of dr, dv and u of the previous segment

        // ControlApprox = 2;
        // Case 2: The propagation is performed with a propagator like RK78 to a variable time tca as a DA variable (Picard-Lindelöf)

        DA tCA_tn;
        tCA_tn = DA(10);

        tCA_tn = tcaInversion(ControlApprox, xp_tn_DA, xs_tn_Vec, tCA_tn, tCA_Nom, MuEarth, Lsc);
        AlgebraicVector<DA>     xp_tnp1_DA(6); 
        AlgebraicVector<double> xs_tnp1_DA(6);

        AlgebraicVector<double> xp_tnp1_Vec(6); 

        double StepSizeN = static_cast<double>(1)/N;
        xp_tnp1_DA  = RK78(6, xp_tn_DA,  {u_tn[0],u_tn[1],u_tn[2]},      0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation from t(n) to t(n+1) of DA object
        xp_tnp1_Vec = RK78(6, xp_tn_Vec, {u_Nom[0], u_Nom[1], u_Nom[2]}, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation from t(n) to t(n+1) of vector object
        xs_tnp1_DA  = RK78(6, xs_tn_Vec, {u_Nom[0], u_Nom[1], u_Nom[2]}, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation from t(n) to t(n+1) of vector object
        // Perhaps only needed for the last time interval

        AlgebraicVector<double> xp_tnp1_Val_T(6), xp_tnp1_Val_R(6);
        xp_tnp1_Val_R = RK78(6, xp_tn_DA.cons(), u_Val_R, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation till tca
        xp_tnp1_Val_T = RK78(6, xp_tn_DA.cons(), u_Val_T, 0.0, tCA_Nom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation till tca

        int DistanceMetric;
        DistanceMetric = 1;
        // Case 1: Distance Metric is Euclidean Distance. 
        // DistanceMetric = 2;
        // Case 2: Distance Metric is Square Mahalanobis Distance (SMD)
        // DistanceMetric = 3;
        // Case 3: Distance Metric is PoC formula (Serra)
        AlgebraicVector<double> u_OptFO_tn(3);
        AlgebraicVector<double> EvalVector(9);
        AlgebraicVector<DA>     EvalVector_Previous(9);
        AlgebraicVector<DA>     EvalVector_Current(9);
        AlgebraicVector<DA>     DeltaRB_tnp1(3),  DeltaV_tnp1(3);
        AlgebraicVector<double> DeltaRB_tnp1_Val_T(3), DeltaV_tnp1_Val_T(3);
        AlgebraicVector<double> DeltaRB_tnp1_Val_R(3), DeltaV_tnp1_Val_R(3);

        switch(DistanceMetric){
            case 1: 
                {
                for(i=0; i<3; i++){
                    // DA Vectors
                    DeltaRB_tnp1[i] = xp_tnp1_DA[i]   - xs_tnp1_DA[i];
                    DeltaV_tnp1[i]  = xp_tnp1_DA[i+3] - xs_tnp1_DA[i+3];

                    // Vectors
                    DeltaRB_tnp1_Val_T[i] =  xp_tnp1_Val_T[i]   - xs_tnp1_DA[i];
                    DeltaV_tnp1_Val_T[i]  =  xp_tnp1_Val_T[i+3] - xs_tnp1_DA[i+3];

                    // Vectors
                    DeltaRB_tnp1_Val_R[i] = xp_tnp1_Val_R[i]   - xs_tnp1_DA[i];
                    DeltaV_tnp1_Val_R[i]  = xp_tnp1_Val_R[i+3] - xs_tnp1_DA[i+3];
                }

                // Now substitute the DA vectors in the SMD, evaluated at tca (it n!=9)
                if(n==9){
                    EucDis         = DeltaRB_tnp1.vnorm();
                    //EucDis = DeltaRBtnp1;
                    //DA EucDis_no_TCA = EucDis;//SMD.substitute(DeltatCAtN9, findTCA(xrelDATCAProp, 10))
                    // EucDis is now a Taylor polynomial of dr, dv and du at the previous time step

                    double EucDis_Val_T;
                    EucDis_Val_T      = DeltaRB_tnp1_Val_T.vnorm();
                    double EucDis_Val_R;
                    EucDis_Val_R      = DeltaRB_tnp1_Val_R.vnorm();

                    // In first-order approximation, the control can be derived from the partial derivative of the DA object EucDis
                    for (i=0; i<3; i++){
                        u_OptFO_tn[i] = cons(EucDis.deriv(7+i));
                    }
                    u_OptFO_tn = u_OptFO_tn/u_OptFO_tn.vnorm(); // Normalise the control; we set the magnitude independently

                    for (i=0; i<6; i++){
                        EvalVector[i] = 0.0;
                    }
                    for (i=0; i<3; i++){
                        EvalVector[i+6] = ThrustMagnitude*u_OptFO_tn[i];
                    }

                    // For next iteration: First 6 are dx(t=n): free DA variables to be substituted next iteration, last 3 are filled in control
                    for (i=0; i<3; i++){
                        EvalVector_Previous[i]   = rp_tn_DA[i] - xp_tn_Vec[i];
                        EvalVector_Previous[i+3] = vp_tn_DA[i] - xp_tn_Vec[i+3];
                        EvalVector_Previous[i+6] = ThrustMagnitude*u_OptFO_tn[i];//ThrustMagnitude*uOptFOtn[i];
                    }
                    
                    EucDis_PreviousIt = EucDis.eval(EvalVector_Previous);

                    if (DebugMode == 1)
                    {
                        cout << "EucDis raw " << EucDis << endl;
                        cout << "Iteration: " << n << endl << endl;
                        cout << "uRTN opt: " << u_OptFO_tn << endl << endl;         
                        cout << "EucDis no control: "  <<  EucDis.cons() << endl << endl; 
                        cout << "EucDis with control: "  << EucDis.eval(EvalVector) << endl << endl;  
                        cout << "EucDis T control: "  << EucDis_Val_T << endl << endl;  
                        cout << "EucDis R control: "  << EucDis_Val_R << endl << endl;  
                    }
                } else {
                    for (i=0; i<6; i++){ // DA: dependent on xn and un
                        EvalVector_Current[i]   = xp_tnp1_DA[i];
                    }
                    for (i=0; i<3; i++){ //DA
                        EvalVector_Current[i+6] = u_tn[i];
                    }

                    // cout << "EvalVectorCurrent: " << EvalVector_Current << endl << endl;
                    // cout << "EucDisNextIt cons: " << EucDis_NextIt.cons() << endl << endl;

                    // In first-order approximation, the control can be derived from the partial derivative of the DA object EucDis
                    
                    EucDis_CurrentIt = EucDis_PreviousIt.eval(EvalVector_Current);
                    for (i=0; i<3; i++){
                        u_OptFO_tn[i] = cons(EucDis_CurrentIt.deriv(7+i));
                    }
                    u_OptFO_tn = u_OptFO_tn/u_OptFO_tn.vnorm(); // Normalise the control; we set the magnitude independently

                    for (i=0; i<6; i++){
                        EvalVector[i] = 0.0;
                    }
                    for (i=0; i<3; i++){
                        EvalVector[i+6] = ThrustMagnitude*u_OptFO_tn[i];
                    }
                    if (DebugMode == 1){
                        cout << "EucDis_CurrentIt DA: " << EucDis_CurrentIt << endl;
                        cout << "Iteration: " << n << endl << endl;
                        cout << "uRTN opt: " << u_OptFO_tn << endl << endl;          // Weird; control seems correct, SMD not       
                        cout << "EucDis no control: "  <<  EucDis.cons() << endl << endl; 
                        cout << "EucDis previous it: "  << EucDis_PreviousIt.cons() << endl << endl;
                        cout << "EucDis with control: "  << EucDis_CurrentIt.eval(EvalVector) << endl << endl;
                    }
                    // For next iteration: First 6 are dx(t=n): free DA variables to be substituted next iteration, last 3 are filled in control
                    for (i=0; i<3; i++){
                        EvalVector_Previous[i]   = rp_tn_DA[i] - xp_tn_Vec[i];
                        EvalVector_Previous[i+3] = vp_tn_DA[i] - xp_tn_Vec[i+3];
                        EvalVector_Previous[i+6] = ThrustMagnitude*u_OptFO_tn[i];//ThrustMagnitude*uOptFOtn[i];
                    }
                    EucDis_PreviousIt = EucDis.eval(EvalVector_Previous);
                    if (DebugMode == 1){
                    cout << "xp_tnp1_DA no control: "  << xp_tnp1_DA.cons() << endl << endl;
                    cout << "xp_tnp1_DA with control: "  << xp_tnp1_DA.eval(EvalVector) << endl << endl;
                    }
                      
                    // cout << "EucDis T control: "  << EucDisValT << endl << endl;  
                    // cout << "EucDis R control: "  << EucDisValR << endl << endl;  
                }

                


                break;
                }
            //case 2: // TO DO ITERATIVE
            //    {
                // Define SMD
                
                // AlgebraicMatrix<double> P(3,3); 

                // P.at(0,0)      = 1;
                // P.at(0,1)      = 0.05; 
                // P.at(0,2)      = 0.04; 
                // P.at(1,0)      = P.at(0,1);
                // P.at(1,1)      = 0.9; 
                // P.at(1,2)      = 0.045; 
                // P.at(2,0)      = P.at(0,2);
                // P.at(2,1)      = P.at(1,2); 
                // P.at(2,2)      = 0.75; 
                // // Symmetric covariance matrix

                // DeltaRBtnp1[0]  = -(xstnp1N[0] - xptnp1N[0]);
                // DeltaRBtnp1[1]  = -(xstnp1N[1] - xptnp1N[1]);
                // DeltaRBtnp1[2]  = -(xstnp1N[2] - xptnp1N[2]);

                // SMD         = dot(DeltaRBtnp1,P * DeltaRBtnp1);

                // DeltaRBtnp1ValT[0]  = -(xstnp1N[0] - xptnp1ValT[0]);
                // DeltaRBtnp1ValT[1]  = -(xstnp1N[1] - xptnp1ValT[1]);
                // DeltaRBtnp1ValT[2]  = -(xstnp1N[2] - xptnp1ValT[2]);

                // DeltaRBtnp1ValR[0]  = -(xstnp1N[0] - xptnp1ValR[0]);
                // DeltaRBtnp1ValR[1]  = -(xstnp1N[1] - xptnp1ValR[1]);
                // DeltaRBtnp1ValR[2]  = -(xstnp1N[2] - xptnp1ValR[2]);

                // double SMDValT, SMDValR;

                // SMDValT = dot(DeltaRBtnp1ValT,P * DeltaRBtnp1ValT);
                // SMDValR = dot(DeltaRBtnp1ValR,P * DeltaRBtnp1ValR);
    
                // // SMD is now a Taylor polynomial of dr, dv and du at the previous time step

                // DA SMD_no_TCA = SMD;//SMD.substitute(DeltatCAtN9, findTCA(xrelDATCAProp, 10))

                // // In first-order approximation, the control can be derived from the partial derivative of the DA object SMD
                
                // for (i=0; i<3; i++){
                //     uOptFO[i] = cons(SMD_no_TCA.deriv(7+i));
                // }
                // uOptFO = uOptFO/uOptFO.vnorm(); // Normalise the control; we set the magnitude independently

                // for (i=0; i<6; i++){
                //     EvalVector[i] = 0;
                // }
                // for (i=0; i<3; i++){
                //     EvalVector[i+6] = ThrustMagnitude*uOptFO[i];
                // }


                // cout << "SMD no control: "  <<  SMD.cons() << endl << endl; 
                // cout << "DeltaR nominally: "  << DeltaRB.cons() << endl << endl;
                // //cout << "Primary nominally: "  << xptfNew.cons() << endl << endl;
                // //cout << "Secondary nominally: "  << xstfNew.cons() << endl << endl;

                // cout << "SMD with control: "  << SMD.eval(EvalVector) << endl << endl;  
                // cout << "uRTN opt: " << uOptFO << endl << endl;      
                // cout << "DeltaR with control: "  << DeltaRB.eval(EvalVector) << endl << endl;
                // //cout << "Primary with control: "  << xptfNew.eval(EvalVector) << endl << endl;
                // cout << "tca change with control: "  << DeltatCAtN9.eval(EvalVector) << endl << endl;  
                // cout << "SMD R control: "  << SMDValR << endl << endl;  
                // cout << "DeltaR with R control: "  << DeltaRBValR << endl << endl; 
                // //cout << "Primary with R control: "  << xpfValR << endl << endl;
                // cout << "SMD T control: "  << SMDValT << endl << endl;  
                // cout << "DeltaR with T control: "  << DeltaRBValT << endl << endl; 
                // //cout << "Primary with T control: "  << xpfValT << endl << endl; 
                // break;
                // }
            }
    }
}
    

    

// Also do with miss distance just to check

    //double SMDtN9 = SMD.cons(0,0,uOptFO); // Evaluate the SMD as a result of the control u(t(n=9))


    // // print initial and final conditions
    
    // cout << endl << "Initial conditions, primary:" << endl << endl;
    // cout << xp0.cons() << endl << endl;
    
    // cout << endl << "Final conditions, primary:" << endl << endl;
    // cout << xpf.cons() << endl << endl;

    // cout << endl << "Initial conditions, secondary:" << endl << endl;
    // cout << xs0.cons() << endl << endl;
    
    // cout << endl << "Final conditions, secondary:" << endl << endl;
    // cout << std::setw(10) << xsf.cons() << endl << endl;

    // Nominal orbit: conjunction at t=7680.61

    // Subtract constant part to build DirMap
    // AlgebraicVector<DA> DirMapP(6);
    // AlgebraicVector<DA> DirMapS(6);
    
    // DirMapP = xpf-xpf.cons();
    // DirMapS = xsf-xsf.cons();

    // // Invert DirMap to obtain InvMap
    // AlgebraicVector<DA> InvMapP(6);
    // AlgebraicVector<DA> InvMapS(6);
    
    // InvMapP = DirMapP.invert();
    // InvMapS = DirMapS.invert();




