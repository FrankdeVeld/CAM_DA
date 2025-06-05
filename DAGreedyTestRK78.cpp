#include <dace/dace.h>
#include "DABackwardFunctions.h"
#include <cmath>
#include <iomanip>
#include <typeinfo>
using namespace std;
using namespace DACE;

int main( void )
{   
    int i;
    DA::init( 2, 10 ); // DA(1:6) = xp(t(N=9)), DA(7:9) = deltaU,  DA(10) = deltaTCA, DA(11)=SMD 
    std::cout.precision(10);
    AlgebraicVector<double> xp0(6), xpf(6);
    AlgebraicVector<double> uNom(3);
    for (i=0; i<3; i++){
        uNom[i]   = 0.0;
    }
    
    
    AlgebraicVector<double> xs0(6), xsf(6);
    AlgebraicVector<double> xpfValT(6), xpfValR(6);
    double MuEarth = 398600;  // km^3/s^2
    double Lsc     = 1;
    double tcaNom;
    double ThrustMagnitude = 1e-7; // km/s^3
    AlgebraicVector<double> uValR(3);
    AlgebraicVector<double> uValT(3);

    for (i=0; i<3; i++){
        uValR[i]   = 0.0;
        uValT[i]   = 0.0;
    }
    uValR[0] = ThrustMagnitude;
    uValT[1] = ThrustMagnitude;
    int Scenario = 3;

    xp0    = InitialXp(Scenario, MuEarth);
    xs0    = InitialXs(Scenario, MuEarth);
    tcaNom = Initialtca(Scenario, MuEarth);

    xpf = RK78(6, xp0, uNom, 0.0, tcaNom,TBAcc,MuEarth,Lsc); // Propagation till tca
    xsf = RK78(6, xs0, uNom, 0.0, tcaNom,TBAcc,MuEarth,Lsc);

    int N;
    N = 10;

    /////// Discretisation: first consider N=10 timesteps, look at N=9, t(N=9)=6912.549 s (scenario 1)

    // Initialise 'main' DA object, as well as various parameters 'becoming' a DA object as a result
    //AlgebraicVector<DA> DAVar(9);
    AlgebraicVector<DA> rptN9(3); // Primary position at t(N=9)
    AlgebraicVector<DA> vptN9(3); // Primary velocity at t(N=9)
    AlgebraicVector<DA>  utN9(3); // Control between t(N=9) and t(N=10)
    
    // Obtain the state of the primary, secondary as DA objects at t(N=9)
    // High-fidelity RK78 propagator - no DA objects are used yet
    AlgebraicVector<double> xptN9Vec(6), xstN9Vec(6);
    double Nmin1Quot = static_cast<double>(N-1)/N;
    xptN9Vec = RK78(6, xp0, uNom, 0.0, tcaNom*Nmin1Quot,TBAcc,MuEarth,Lsc); // Propagation from t(N=0) to t(N=9)
    xstN9Vec = RK78(6, xs0, uNom, 0.0, tcaNom*Nmin1Quot,TBAcc,MuEarth,Lsc);

    // Make rp, vp, u @ t(N=9) variable DA objects 
    for (i=0; i<3; i++){
        rptN9[i] = xptN9Vec[i] + DA(i+1);
        vptN9[i] = xptN9Vec[i+3] + DA(i+4);
        utN9[i]  = DA(i+7);  
    }
    DA DeltatCAtN9;
    DeltatCAtN9 = DA(10);

    // Propagate rptN9, vptN9 to t(N+10) + delta tCA
    AlgebraicVector<DA> xptN9DA(6);
    for (i=0; i<3; i++){
        xptN9DA[i]   = rptN9[i];
        xptN9DA[i+3] = vptN9[i];
    }

    int ControlApprox;
    ControlApprox = 1; 
    // Case 1: Approximation: no control between tcaNom and tca(u). Effectively: no change in tca
    // As such, a RK78 propagation of DA variables is performed till tcaNom, after which a Kepler propagation 
    // is done around tcaNom with tf a DA variable. Polynomial inversion is done to obtain tca as a 
    // Taylor polynomial of dr, dv and u of the previous segment

    // ControlApprox = 2;
    // Case 2: The propagation is performed with a propagator like RK78 to a variable time tca as a DA variable

    DeltatCAtN9 = tcaInversion(ControlApprox, xptN9DA, xstN9Vec, DeltatCAtN9, tcaNom, MuEarth, Lsc);
    AlgebraicVector<DA>     xptfNew(6); 
    AlgebraicVector<double> xstfNew(6);
    double StepSizeN = static_cast<double>(1)/N;
    xptfNew = RK78(6, xptN9DA, {DA(7),DA(8),DA(9)}, 0.0, tcaNom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation from t(N=9) to t(N=10) of DA object
    xstfNew = RK78(6, xstN9Vec, {0.0, 0.0, 0.0}, 0.0, tcaNom*StepSizeN,TBAcc,MuEarth,Lsc);               // Propagation from t(N=9) to t(N=10) of vector object

    xpfValR = RK78(6, cons(xptN9DA), uValR, 0.0, tcaNom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation till tca
    xpfValT = RK78(6, cons(xptN9DA), uValT, 0.0, tcaNom*StepSizeN,TBAcc,MuEarth,Lsc); // Propagation till tca

    int DistanceMetric;
    DistanceMetric = 1;
    // Case 1: Distance Metric is Euclidean Distance. 
    // DistanceMetric = 2;
    // Case 2: Distance Metric is Square Mahalanobis Distance (SMD)
    // DistanceMetric = 3;
    // Case 3: Distance Metric is PoC formula (Serra)
    AlgebraicVector<double> uOptFO(3);
    AlgebraicVector<double> EvalVector(9);
    AlgebraicVector<DA>     DeltaRB(3);
    AlgebraicVector<double>     DeltaRBValT(3);
    AlgebraicVector<double>     DeltaRBValR(3);
    switch(DistanceMetric){
        case 1: 
            {
            // Define EucDis
            DA EucDis;  

            DeltaRB[0]  = -(xstfNew[0] - xptfNew[0]);
            DeltaRB[1]  = -(xstfNew[1] - xptfNew[1]);
            DeltaRB[2]  = -(xstfNew[2] - xptfNew[2]);

            EucDis         = DeltaRB.vnorm();
            DA EucDis_no_TCA = EucDis;//SMD.substitute(DeltatCAtN9, findTCA(xrelDATCAProp, 10))
            // EucDis is now a Taylor polynomial of dr, dv and du at the previous time step

            DeltaRBValT[0]  = -(xstfNew[0] - xpfValT[0]);
            DeltaRBValT[1]  = -(xstfNew[1] - xpfValT[1]);
            DeltaRBValT[2]  = -(xstfNew[2] - xpfValT[2]);

            DeltaRBValR[0]  = -(xstfNew[0] - xpfValR[0]);
            DeltaRBValR[1]  = -(xstfNew[1] - xpfValR[1]);
            DeltaRBValR[2]  = -(xstfNew[2] - xpfValR[2]);

            double EucDisValT;
            EucDisValT      = DeltaRBValT.vnorm();
            double EucDisValR;
            EucDisValR      = DeltaRBValR.vnorm();

            // In first-order approximation, the control can be derived from the partial derivative of the DA object EucDis
            for (i=0; i<3; i++){
                uOptFO[i] = cons(EucDis_no_TCA.deriv(7+i));
            }
            uOptFO = uOptFO/uOptFO.vnorm(); // Normalise the control; we set the magnitude independently

            for (i=0; i<6; i++){
                EvalVector[i] = 0;
            }
            for (i=0; i<3; i++){
                EvalVector[i+6] = ThrustMagnitude*uOptFO[i];
            }
            cout << "uRTN opt: " << uOptFO << endl << endl;         
            cout << "EucDis no control: "  <<  EucDis.cons() << endl << endl; 
            cout << "EucDis with control: "  << EucDis.eval(EvalVector) << endl << endl;  
            cout << "EucDis T control: "  << EucDisValT << endl << endl;  
            cout << "EucDis R control: "  << EucDisValR << endl << endl;  
            break;
            }
        case 2:
            {
            // Define SMD
            DA SMD;  
            AlgebraicMatrix<double> P(3,3); 

            P.at(0,0)      = 1;
            P.at(0,1)      = 0.05; 
            P.at(0,2)      = 0.04; 
            P.at(1,0)      = P.at(0,1);
            P.at(1,1)      = 0.9; 
            P.at(1,2)      = 0.045; 
            P.at(2,0)      = P.at(0,2);
            P.at(2,1)      = P.at(1,2); 
            P.at(2,2)      = 0.75; 
            // Symmetric covariance matrix

            DeltaRB[0]  = -(xstfNew[0] - xptfNew[0]);
            DeltaRB[1]  = -(xstfNew[1] - xptfNew[1]);
            DeltaRB[2]  = -(xstfNew[2] - xptfNew[2]);

            SMD         = dot(DeltaRB,P * DeltaRB);

            DeltaRBValT[0]  = -(xstfNew[0] - xpfValT[0]);
            DeltaRBValT[1]  = -(xstfNew[1] - xpfValT[1]);
            DeltaRBValT[2]  = -(xstfNew[2] - xpfValT[2]);

            DeltaRBValR[0]  = -(xstfNew[0] - xpfValR[0]);
            DeltaRBValR[1]  = -(xstfNew[1] - xpfValR[1]);
            DeltaRBValR[2]  = -(xstfNew[2] - xpfValR[2]);

            double SMDValT, SMDValR;

            SMDValT = dot(DeltaRBValT,P * DeltaRBValT);
            SMDValR = dot(DeltaRBValR,P * DeltaRBValR);
   
            // SMD is now a Taylor polynomial of dr, dv and du at the previous time step

            DA SMD_no_TCA = SMD;//SMD.substitute(DeltatCAtN9, findTCA(xrelDATCAProp, 10))

            // In first-order approximation, the control can be derived from the partial derivative of the DA object SMD
             
            for (i=0; i<3; i++){
                uOptFO[i] = cons(SMD_no_TCA.deriv(7+i));
            }
            uOptFO = uOptFO/uOptFO.vnorm(); // Normalise the control; we set the magnitude independently

            for (i=0; i<6; i++){
                EvalVector[i] = 0;
            }
            for (i=0; i<3; i++){
                EvalVector[i+6] = ThrustMagnitude*uOptFO[i];
            }


            cout << "SMD no control: "  <<  SMD.cons() << endl << endl; 
            cout << "DeltaR nominally: "  << DeltaRB.cons() << endl << endl;
            //cout << "Primary nominally: "  << xptfNew.cons() << endl << endl;
            //cout << "Secondary nominally: "  << xstfNew.cons() << endl << endl;

            cout << "SMD with control: "  << SMD.eval(EvalVector) << endl << endl;  
            cout << "uRTN opt: " << uOptFO << endl << endl;      
            cout << "DeltaR with control: "  << DeltaRB.eval(EvalVector) << endl << endl;
            //cout << "Primary with control: "  << xptfNew.eval(EvalVector) << endl << endl;
            cout << "tca change with control: "  << DeltatCAtN9.eval(EvalVector) << endl << endl;  
            cout << "SMD R control: "  << SMDValR << endl << endl;  
            cout << "DeltaR with R control: "  << DeltaRBValR << endl << endl; 
            //cout << "Primary with R control: "  << xpfValR << endl << endl;
            cout << "SMD T control: "  << SMDValT << endl << endl;  
            cout << "DeltaR with T control: "  << DeltaRBValT << endl << endl; 
            //cout << "Primary with T control: "  << xpfValT << endl << endl; 
            break;
            }
    }
    

    

// Also do with miss distance just to check!! 

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
    
    
    
    
    
}



