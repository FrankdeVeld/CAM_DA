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

    int N            = jin.at("N").get<int>();
    int DM_Case      = jin.at("metric_case").get<int>();        // 1=Euclidean, 2=SMD
    int tCAHandling  = jin.at("tCAHandling").get<int>();
    int order        = jin.at("order").get<int>();        // 1=Euclidean, 2=SMD
    double t_back    = jin.at("t_back").get<double>();   // [s], dimensional
    double uMax      = jin.at("uMax").get<double>();       // [km/s^2], dimensional thrust magnitude
    double Lsc       = jin.at("Lsc").get<double>();        // [km]
    double R         = jin.at("HBR").get<double>();        // hard-body radius [km]
    double lim       = jin.at("lim").get<double>();        
    double smd0      = jin.at("smd0").get<double>();        
    
    // States at TCA in dimensional units [km, km/s], 6-element arrays
    auto xptf_std = jin.at("xp_tCA").get<vector<double>>();
    auto xstf_std = jin.at("xs_tCA").get<vector<double>>();
    auto rb0      = jin.at("rb0").get<vector<double>>();

    // 3x3 covariance in B-plane reference, row-major nested array
    auto P_std = jin.at("P").get<vector<vector<double>>>();
    AlgebraicMatrix<double> P(3,3);
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            P.at(i,j) = P_std[i][j];
    //////////////////////////////////////////////////////////////// SETUP (unchanged from original) //////////////////////////////////////////////////////////////
    int i, j;
    DA::init( order, 10 );
    cout.precision(16);

    // Fill AlgebraicVector states from JSON
    AlgebraicVector<double> xp_tf(6), xs_tf(6), rB0(2);
    for (i=0;i<6;i++) {
        xp_tf[i] = xptf_std[i];
        xs_tf[i] = xstf_std[i];
    }
    rB0[0] = rb0[0]; rB0[1] = rb0[1];

    // Propagate backwards to t=0 to get initial states
    AlgebraicVector<double> u_Nom = {0.0, 0.0, 0.0};
    double dt = t_back/(N-1);

    DA DM_tn, tCA_tn;
    AlgebraicVector<DA> DeltaRB_tn(3);

    AlgebraicVector<DA>     u_tn(3), xp_tn_DA(6), rp_tn_DA(3), vp_tn_DA(3);
    AlgebraicVector<DA>     xp_tn_Vec(6), xp_tnp1_DA(6);
    AlgebraicVector<DA>     xs_tCA_DA(6), xp_tCA_DA(6);
    AlgebraicVector<DA>     Evaluated_tCA(10);

    // Save matrices (same as original, kept for intermediate use)
    AlgebraicMatrix<double> u_save(N,3), DeltaRB_save(N,2);
    AlgebraicVector<double> DM_save(N), tCA_save(N), alpha_vec(N);

    DA DM_NextIt, tCA_NextIt;
    AlgebraicVector<DA> DeltaRB_NextIt;
    double DM_Evaluated_Control, tCA_Evaluated_Control;
    AlgebraicVector<double> xp_tnp1_Evaluated_Control(6), u_OptFO_tn(3), DeltaRB_Evaluated_Control(3);
    
    DeltaRB_save.at(N-1,0) = rb0[0];
    DeltaRB_save.at(N-1,1) = rb0[1];
    u_save.at(N-1,0)       = 0.0; 
    u_save.at(N-1,1)       = 0.0; 
    u_save.at(N-1,2)       = 0.0;
    tCA_save[N-1]          = 0.0;
    if (DM_Case == 1)  {
        DM_save[N-1] = rB0.dot(rB0);
    }
    else {
        DM_save[N-1] = smd0;
    }
    DA da_null = DA(1)*0.0;
    DA alpha = DA(10);
    double t_n, alpha_ev;
    ////////////////////////////////////////////////////////////// START OF ITERATIVE LOOP /////////////////////////////////////////////////////
    for(int n = N-1; n > 0; n--) {
        t_n = dt*(N-n);
        xp_tn_Vec = KeplerProp(xp_tf + da_null, - (t_n - alpha), 1.0);

        for (i=0; i<3; i++){
            rp_tn_DA[i] = xp_tn_Vec[i]   + DA(i+1);
            vp_tn_DA[i] = xp_tn_Vec[i+3] + DA(i+4);
            u_tn[i]     = DA(i+7);
        }
        for (i=0; i<3; i++){
            xp_tn_DA[i]   = rp_tn_DA[i];
            xp_tn_DA[i+3] = vp_tn_DA[i];
        }

        xp_tnp1_DA  = RK78(6, xp_tn_DA, {u_Nom[0]+u_tn[0], u_Nom[1]+u_tn[1], u_Nom[2]+u_tn[2]}, da_null, dt - alpha, TBAcc, 1.0, Lsc);

        if (n == N-1){
            tCA_tn      = DA(10);
            tie(tCA_tn, xp_tCA_DA, xs_tCA_DA) = tcaInversion(tCAHandling, u_Nom, u_tn, xp_tnp1_DA, xs_tf, tCA_tn, 1.0, Lsc);
            for (i=0;i<6;i++) Evaluated_tCA[i] = DA(i+1);
            for (i=0;i<3;i++) Evaluated_tCA[i+6] = DA(i+7);
            Evaluated_tCA[9] = tCA_tn;
            xp_tnp1_DA = xp_tCA_DA.eval(Evaluated_tCA);
            xs_tCA_DA = xs_tCA_DA.eval(Evaluated_tCA);
        } 
        tie(xp_tnp1_Evaluated_Control, u_OptFO_tn, DeltaRB_Evaluated_Control, DM_Evaluated_Control, tCA_Evaluated_Control, DM_NextIt, tCA_NextIt, DeltaRB_NextIt, alpha_ev, breakyes) =
                IterativeDA(n, N, DM_Case, xp_tnp1_DA, xs_tCA_DA, uMax, DM_tn, tCA_tn, DeltaRB_tn, P, R, lim);
        DM_tn      = DM_NextIt;
        tCA_tn     = tCA_NextIt;
        DeltaRB_tn = DeltaRB_NextIt;

        for(i=0;i<3;i++){
            u_save.at(n-1,i)       = u_OptFO_tn[i];
        }
        DeltaRB_save.at(n-1,0) = DeltaRB_Evaluated_Control[0];
        DeltaRB_save.at(n-1,1) = DeltaRB_Evaluated_Control[1];

        DM_save[n-1]  = DM_Evaluated_Control;
        tCA_save[n-1] = tCA_Evaluated_Control;
        alpha_vec[n-1] = 0.0;
        if (breakyes == 1) {alpha_vec[n-1] = alpha_ev; break;}
    }

    //////////////////////////////////////////////////////////////// WRITE JSON OUTPUT ////////////////////////////////////////////////////////////////////////////
    json jout;
    jout["timeMs"] = now_ms() - t0ms;
    jout["N"]      = N;

    json j_nodes = json::array();
    for (int n = 0; n < N; ++n) {
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