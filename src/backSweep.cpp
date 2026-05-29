#include <dace/dace.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <typeinfo>
#include <nlohmann/json.hpp>
#include "DABackwardFunctionsIterative.h"
#include "prop_utils.h"

using namespace std;
using namespace DACE;
using json = nlohmann::json;

// ---------------------------------------------------------------------------
// Minimal JSON helpers (inline so no extra .cpp is needed)
// ---------------------------------------------------------------------------
static json read_json(const string& path) {
    ifstream f(path);
    if (!f.is_open()) throw runtime_error("Cannot open " + path);
    json j; f >> j; return j;
}

static void write_json(const string& path, const json& j) {
    ofstream f(path);
    if (!f.is_open()) throw runtime_error("Cannot write " + path);
    f << j.dump(2) << "\n";
}

// ---------------------------------------------------------------------------
// Simple wall-clock helper
// ---------------------------------------------------------------------------
#include <chrono>
static long now_ms() {
    return (long)chrono::duration_cast<chrono::milliseconds>(
        chrono::steady_clock::now().time_since_epoch()).count();
}

int main( void )
{
    long t0ms = now_ms();
    bool breakyes;
    //////////////////////////////////////////////////////////////// READ JSON INPUT ///////////////////////////////////////////////////////////////////////////////
    json jin = read_json("./input.json");
    // everything is nondimensionalised before json writing, so all the script works in nondimensional units
    int N                   = jin.at("N").get<int>();                   // Number of nodes
    int DM_Case             = jin.at("metric_case").get<int>();         // 1=MD, 2=SMD
    int tCAHandling         = jin.at("tCAHandling").get<int>();         // 1=no, 2=DA, 3=P-L
    int order               = jin.at("order").get<int>();               // DA order
    int refineLastInterval  = jin.at("refineLastInterval").get<int>();  // 0=no, 1=yes
    double t_back           = jin.at("t_back").get<double>();           // backpropagation time
    double t_start          = jin.at("t_start").get<double>();          // maximum maneuver end time (not implemented)
    double uMax             = jin.at("uMax").get<double>();             // nondimensional thrust magnitude
    double Lsc              = jin.at("Lsc").get<double>();              // [km] position scaling (typically semi-major axis)
    double R                = jin.at("HBR").get<double>();              // hard-body radius
    double lim              = jin.at("lim").get<double>();              // metric limit      
    double smd0             = jin.at("smd0").get<double>();             // ballistic smd at TCA 
    
    // States at TCA, 6-element arrays
    auto xptf_std = jin.at("xp_tCA").get<vector<double>>();
    auto xstf_std = jin.at("xs_tCA").get<vector<double>>();
    auto rb0      = jin.at("rb0").get<vector<double>>();

    // 3x3 covariance in B-plane reference, row-major nested array
    auto P_std = jin.at("P").get<vector<vector<double>>>();
    AlgebraicMatrix<double> DeltaRB_save(N,2), u_save(N,3), P(3,3);
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            P.at(i,j) = P_std[i][j];

    int i, j;

    // Initialize DA framework
    DA::init( order, 10 );
    DA::setEps(1e-15);

    // Fill initial conditions from JSON
    AlgebraicVector<double> xp_tf(6), xs_tf(6), rB0(2), DM_save(N), tCA_save(N), alpha_vec(N), xp_tnp1_Evaluated_Control(6), u_OptFO_tn(3), DeltaRB_Evaluated_Control(2);
    for (i=0;i<6;i++) {
        xp_tf[i] = xptf_std[i];
        xs_tf[i] = xstf_std[i];
    }
    rB0[0] = rb0[0]; rB0[1] = rb0[1];

    // Null control (nominal)
    AlgebraicVector<double> u_Nom = {0.0, 0.0, 0.0};
    
    // constant Delta t between two time nodes 
    double dt = (t_back-t_start)/(N-1);

    // Initialize DA variables
    DA DM_tn, tCA_tn, DM_NextIt, tCA_NextIt;
    AlgebraicVector<DA> u_tn(3), xp_tn_DA(6), xp_tn_Vec(6), xp_tnp1_DA(6), xs_tCA_DA(6), xp_tCA_DA(6), Evaluated_tCA(10), DeltaRB_tn(2), DeltaRB_NextIt(2);

    // Initialize double variables
    double DM_Evaluated_Control, tCA_Evaluated_Control, t_n, alpha_ev;;
    
    // Save position in B-plane and initial tca shift before maneuver execution
    DeltaRB_save.at(N-1,0) = rb0[0];
    DeltaRB_save.at(N-1,1) = rb0[1];
    tCA_save[N-1]          = 0.0;
    if (DM_Case == 1)  {
        DM_save[N-1] = rB0.dot(rB0);
    }
    else {
        DM_save[N-1] = smd0;
    }
    // dummy DA variable for null DA action (needed to use a DA type in the called functions)
    DA da_null = DA(1)*0.0;
    
    // time node shrinking variable (only used for the last time step in the refinement process)
    DA alpha = DA(10);

    // Backward sweep
    for(int n = N-1; n > 0; n--) 
    {
        t_n = dt*(N-n);
        // backpropagation (eq. 34 from the paper)
        xp_tn_Vec = KeplerProp(xp_tf + da_null, - (t_n - alpha), 1.0);
        
        // DA expansion (eq. 35 from the paper)
        for (i=0; i<3; i++){
            xp_tn_DA[i]   = xp_tn_Vec[i]   + DA(i+1);
            xp_tn_DA[i+3] = xp_tn_Vec[i+3]   + DA(i+4);
            u_tn[i]       = DA(i+7);
        }
        
        // DA forward propagation (eq. 36 from the paper)
        xp_tnp1_DA  = RK78(6, xp_tn_DA, {u_Nom[0]+u_tn[0], u_Nom[1]+u_tn[1], u_Nom[2]+u_tn[2]}, da_null, dt - alpha, TBAcc, 1.0, Lsc);

        if (n == N-1){
            tCA_tn      = DA(10);
            // In the last time step, TCA is expanded only if TCAHandling is not 0
            if (tCAHandling > 0){
                // TCA expansion (eq. 27 from the paper)
                tie(tCA_tn, xp_tCA_DA, xs_tCA_DA) = tcaInversion(tCAHandling, u_Nom, u_tn, xp_tnp1_DA, xs_tf, tCA_tn, 1.0, Lsc);
                for (i=0;i<6;i++) Evaluated_tCA[i] = DA(i+1);
                for (i=0;i<3;i++) Evaluated_tCA[i+6] = DA(i+7);
                Evaluated_tCA[9] = tCA_tn;
                xp_tnp1_DA = xp_tCA_DA.eval(Evaluated_tCA);
                xs_tCA_DA = xs_tCA_DA.eval(Evaluated_tCA);
            }
        }
        
        // Computation of new DA maps for next section of the sweep and current DA map to evaluate the metric and check safety 
        tie(xp_tnp1_Evaluated_Control, u_OptFO_tn, DeltaRB_Evaluated_Control, DM_Evaluated_Control, tCA_Evaluated_Control, DM_NextIt, tCA_NextIt, DeltaRB_NextIt, alpha_ev, breakyes) =
                IterativeDA(n, N, DM_Case, xp_tnp1_DA, xs_tCA_DA, uMax, DM_tn, tCA_tn, DeltaRB_tn, P, R, lim, refineLastInterval);
        
        // Save DA maps for the next section of the sweep
        DM_tn      = DM_NextIt;
        tCA_tn     = tCA_NextIt;
        DeltaRB_tn = DeltaRB_NextIt;

        // Save variables at the current section of the sweep
        for(i=0;i<3;i++){
            u_save.at(n-1,i)    = u_OptFO_tn[i];
        }
        DeltaRB_save.at(n-1,0)  = DeltaRB_Evaluated_Control[0];
        DeltaRB_save.at(n-1,1)  = DeltaRB_Evaluated_Control[1];
        DM_save[n-1]            = DM_Evaluated_Control;
        tCA_save[n-1]           = tCA_Evaluated_Control;
        alpha_vec[n-1]          = 0.0;
        // if the section includes the passing of the limit, save the time interval shrinking variable
        if (breakyes == 1) {alpha_vec[n-1] = alpha_ev; break;}
    }
    // Control at TCA (first order hold, just for plotting)
    u_save.at(N-1,0)       = u_save.at(N-2,0); 
    u_save.at(N-1,1)       = u_save.at(N-2,1);  
    u_save.at(N-1,2)       = u_save.at(N-2,2); 
    
    //////////////////////////////////////////////////////////////// WRITE JSON OUTPUT ////////////////////////////////////////////////////////////////////////////
    json jout;
    jout["timeMs"] = now_ms() - t0ms;
    jout["N"]      = N;

    json j_nodes = json::array();
    for (int n = 0; n < N; ++n) 
    {
        json jn;
        jn["i"] = n;
        jn["tNode"] = dt*(N-1-n) - alpha_vec[n];
        jn["alpha"] =  alpha_vec[n];

        jn["controlRTN"] = {
            u_save.at(n,0),
            u_save.at(n,1),
            u_save.at(n,2)
        };

        jn["tcaShift"]    = tCA_save[n];         // delta w.r.t. nominal TCA [s]

        jn["relativePositionBPlane"] = {
            DeltaRB_save.at(n,0),
            DeltaRB_save.at(n,1)
        };

        jn["dangerMetric"] = DM_save[n];
        j_nodes.push_back(jn);
    }
    jout["nodes"] = j_nodes;

    write_json("output.json", jout);
    return 0;
}