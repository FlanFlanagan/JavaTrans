#include <Eigen/Eigen>
#include <tr1/cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <tr1/cmath>

#define MaxLN 9.0
#define bins 10
#define PI 3.14159265359
#define N_A 6.0221415 * pow(10.0, 23)
#define cm2_to_barns pow(10.0, 24)
#define barn_to_cm2 pow(10.0, -24)
#define sec_per_day 24.0 * 3600.0
#define ffactors 6
#define Dimensions 1
#define CONVERGENCE 0.005
#define LEGENDRE 9
#define Ordinate 8
#define SOURCE pow(10.0, 8.0)
#define BENCHMARK 1
#define region_one_size 10.0
#define num_regions 1
#define region_two_size 1
#define Burnup 60.0
#define AVAIL 0.90
#define TIME_STEP 43200
#define ITERATIONS 101
#define Convert 1E-24
#define POWER 1000.;

using namespace std;
using namespace Eigen;

class Isots
{
    public:
        VectorXd CHI;
        RowVectorXd NUFISSION;
        VectorXd TOTAL;
        vector<vector<vector<double> > > SKERNAL;
        int NAME;
        double NumDen;
        int nuclide_type;
        vector<vector<double> > FFACTOR;
        vector<double> New_FFACTOR;
        vector<double> sigma_b;
        vector<double> SCATTER;
        vector<double> ABSORB;
        double CHANGE;
        double ABSORB_C;
        double NUFISSION_C;
    protected:

    private:

};

class mesh
{
    public:
    double x_lower;
    double x_upper;
    double x_int;
    double y_lower;
    double y_upper;
    double z_lower;
    double z_upper;
    vector<vector<vector<vector<double > > > > source;
    vector<vector<vector<vector<vector<double > > > > > flux;
    vector<vector<vector<vector<vector<double > > > > > total_flux;
    vector<vector<vector<vector<vector<double > > > > > old_total_flux;
    vector<double> final_flux_e;
    double final_flux;
};

class Region
{
    public:
        double x_upper;
        double x_lower;
        double y_upper;
        double y_lower;
        double z_upper;
        double z_lower;
        vector<mesh> meshpoints;
        vector<Isots> vIsos;
        double r_type;
        double max_mesh;
        double mesh_size;
        double mesh_points;
        vector<vector<vector<double> > > SKERNAL;
        vector<double> TOTAL;
        vector<double> ABSORB;
        double TOTAL_C;
        double ABSORB_C;
        double NUFISSION_C;
        vector<double> FINAL_FLUX;
        double FINAL_ff;
        double Flux_ff;
        vector<double> SKATTER;
        vector<double> CHI;
        Matrix<double, bins, bins> fission;
        double fission_sum;
        double power_flux;
    protected:
    private:
};


int main()
{
    int total_count = 0;
    double crit1 = 1.0;
    double crit2;
    vector<Region> REGIONS;
    ofstream Isos3;
    ofstream Isos4;
    ofstream Isos5;
    ofstream Isos6;
    Isos3.open("Flux.txt");
    Isos4.open("numbden.txt");
    Isos5.open("crit.txt");
    Isos6.open("ordinateF.txt");
            /// Dimensionality ///
            double EDGES;
            EDGES = 2*(Dimensions)+1;
            /// Ordinates ///
            double mewA2[2] = {0.5774, -0.5774};
            double wmewA2[2] ={1,1};
            double mewA8[8] = {0.9603, 0.7967, 0.5255, 0.1834, -0.1834, -0.5255, -0.7967, -0.9603};
            double wmewA8[8] ={0.1012, 0.2224, 0.3137, 0.3627, 0.3627, 0.3137, 0.2224, 0.1012};
            vector<double>  mew;
            mew.resize(Ordinate);
            vector<double>  Wmew;
            Wmew.resize(Ordinate);
            if(Dimensions == 1){
                if(Ordinate == 2){
                    for(int i = 0; i < 2; i++){
                        mew[i] = mewA2 [i];
                        Wmew[i] = wmewA2[i];
                    }
                }
                else if(Ordinate == 8){
                    for(int i = 0; i < 8; i++){
                        mew[i] = mewA8[i];
                        Wmew[i] = wmewA8[i];
                    }
                }
            }
            /// Regions ///
            if(total_count == 0){
                Region region;
                REGIONS.push_back(region);
                REGIONS[0].x_lower = 0;
                REGIONS[0].x_upper = region_one_size;
                REGIONS[0].r_type = 1;
            /// Setting up the Regions for Benchmark ///
                if(BENCHMARK == 1){
                    if(total_count == 0){
                        for(int iv = 1; iv < num_regions+1; iv ++){
                            REGIONS.push_back(region);
                            REGIONS[iv].x_lower = region_one_size + (region_two_size * (iv-1));
                            REGIONS[iv].x_upper = region_one_size + (region_two_size * iv);
                            REGIONS[iv].r_type = 2;
                        }
                    }
                }
            }
            /// Building the regions ///
            for(unsigned int ri = 0; ri < REGIONS.size(); ri++){
                string region_type;
                string s1 = "./Isotopes2/";
                string s2 = "/Isos.txt";
                stringstream region_num;
                region_num << REGIONS[ri].r_type;
                region_type = s1+region_num.str()+s2;
                ifstream Isos;
                Isos.open(region_type.c_str());
                if(!Isos){
                    cout << "Opps text input broke" << endl;
                }
                Isots nuclide;
                int xx;
                int nuc_type;
                double dens;
                while(Isos >> xx){
                    nuclide.NAME = xx;
                    Isos >> dens;
                    nuclide.NumDen = dens*1E-24;
                    Isos >> nuc_type;
                    nuclide.nuclide_type = nuc_type;
                    REGIONS[ri].vIsos.push_back(nuclide);
                }
                Isos.close();
                /// Set up sizes - Region KERNAL. (HEIGHT x WIDTH)
                REGIONS[ri].SKERNAL.resize(bins);
                for (int i = 0; i < bins; ++i){
                    REGIONS[ri].SKERNAL[i].resize(bins);
                    for (int j = 0; j < bins; ++j){
                        REGIONS[ri].SKERNAL[i][j].resize(MaxLN);
                    }
                }
                /// Set up sizes - Materials KERNAL. (HEIGHT x WIDTH)
                for(unsigned int k = 0; k < REGIONS[ri].vIsos.size(); k++){
                    REGIONS[ri].vIsos[k].SKERNAL.resize(bins);
                    for (int i = 0; i < bins; ++i){
                        REGIONS[ri].vIsos[k].SKERNAL[i].resize(bins);
                        for (int j = 0; j < bins; ++j){
                            REGIONS[ri].vIsos[k].SKERNAL[i][j].resize(MaxLN);
                        }
                    }
                }
                /// Resizing FFRACTOR
                for(int i=0; i < REGIONS[ri].vIsos.size(); i++){
                    REGIONS[ri].vIsos[i].FFACTOR.resize(bins);
                    for(int j=0; j < bins; j++){
                        REGIONS[ri].vIsos[i].FFACTOR[j].resize(ffactors);
                    }
                }

                ///Inputs from cross section files///
                ifstream Isos2;
                for(int i = 0; i < REGIONS[ri].vIsos.size(); i++){
                    string result;
                    string s1 = "./Isotopes2/Isos/";
                    string s2 = ".xs";
                    stringstream isotope_num;
                    isotope_num << REGIONS[ri].vIsos[i].NAME;
                    result = s1+isotope_num.str()+s2;
                    Isos2.open(result.c_str());
                    //cout << result.c_str() << endl;
                    if(!Isos2){
                        cout << "failures!" << endl;
                        cout << result.c_str()<< endl;
                    }
                    string TEXT;
                    REGIONS[ri].vIsos[i].TOTAL.resize(bins);
                    REGIONS[ri].vIsos[i].CHI.resize(bins);
                    REGIONS[ri].vIsos[i].NUFISSION.resize(bins);
                    for(int eng = 0; eng < bins; eng++){
                        REGIONS[ri].vIsos[i].TOTAL[eng] = 0;
                        REGIONS[ri].vIsos[i].CHI[eng] = 0;
                        REGIONS[ri].vIsos[i].NUFISSION[eng] = 0;
                    }

                    double xx;

                    int bin;
                    double LeGrN;
                    int to_grp;
                    int from_grp;
                    while(Isos2){
                        Isos2 >> TEXT;
                        while(Isos2 >> bin){
                            Isos2 >> xx;
                            REGIONS[ri].vIsos[i].TOTAL(10 - bin) = xx;
                        }

                        Isos2.clear();
                        Isos2 >> TEXT;
                       /// FFACTOR INPUT
                       while(Isos2 >> bin){
                            for(int l = 0; l < ffactors; l++){
                                Isos2 >> xx;
                                REGIONS[ri].vIsos[i].FFACTOR[10-bin][l] = xx;
                            }
                        }
                        Isos2.clear();
                        Isos2 >> TEXT;
                        while(Isos2 >> bin){
                            Isos2 >> xx;
                            REGIONS[ri].vIsos[i].CHI(10 - bin) = xx;
                        }
                        Isos2.clear();
                        Isos2 >> TEXT;

                        while(Isos2 >> bin){
                            Isos2 >> xx;
                            REGIONS[ri].vIsos[i].NUFISSION(10 - bin) = xx;
                        }
                        Isos2.clear();
                        Isos2 >> TEXT;

                        while(Isos2 >> to_grp >> from_grp){
                            for(int n = 0; n < MaxLN; n++){
                                Isos2 >> LeGrN;
                                REGIONS[ri].vIsos[i].SKERNAL[10 - from_grp][10 - to_grp][n] = LeGrN;
                            }
                        }
                    }
                    Isos2.clear();
                    Isos2.close();
                }
                /// Input is finished ///
                MatrixXd sigma_total (bins,bins);
                sigma_total = MatrixXd::Zero(bins, bins);
                MatrixXd sigma_s (bins, bins);
                sigma_s = MatrixXd::Zero(bins,bins);
                VectorXd sigma_t_sum (bins);
                sigma_t_sum = VectorXd::Zero(bins);
                REGIONS[ri].TOTAL.resize(bins);
                REGIONS[ri].ABSORB.resize(bins);
                REGIONS[ri].SKATTER.resize(bins);
                REGIONS[ri].CHI.resize(bins);
                /// REGION TOTAL ///
                for(int i = 0;  i < bins; i++){
                    for(int j =0; j < REGIONS[ri].vIsos.size(); j++){
                        REGIONS[ri].TOTAL[i] = REGIONS[ri].TOTAL[i] + REGIONS[ri].vIsos[j].NumDen * REGIONS[ri].vIsos[j].TOTAL(i);
                    }
                }
                /// REGION CHI
                for(int i = 0;  i < bins; i++){
                    REGIONS[ri].CHI[i] = 0;
                    for(int j =0; j < REGIONS[ri].vIsos.size(); j++){
                        REGIONS[ri].CHI[i] = REGIONS[ri].CHI[i] + REGIONS[ri].vIsos[j].CHI(i) / 3.0;
                    }
                }

                /// vIsos SCATTER ///
                for(int k = 0; k < REGIONS[ri].vIsos.size();k++){
                    REGIONS[ri].vIsos[k].SCATTER.resize(bins);
                    REGIONS[ri].vIsos[k].ABSORB.resize(bins);
                }
                for(int i = 0; i < bins; i++){
                    for(int j = 0; j < bins; j++){
                        for(int k=0; k<REGIONS[ri].vIsos.size(); k++){
                            REGIONS[ri].vIsos[k].SCATTER[i] = REGIONS[ri].vIsos[k].SCATTER[i] + REGIONS[ri].vIsos[k].SKERNAL[j][i][0];
                        }
                    }
                }
                /// vIsos ABSORB ///
                for(int i = 0; i < bins; i++){
                    for(int j = 0; j < REGIONS[ri].vIsos.size(); j++){
                        REGIONS[ri].vIsos[j].ABSORB[i] = REGIONS[ri].vIsos[j].TOTAL(i) - REGIONS[ri].vIsos[j].SCATTER[i];
                    }
                }
                /// vIsos Skernal ///
                for(int leg_i = 0; leg_i < LEGENDRE; leg_i++){
                    for(int i = 0; i < bins; i++){
                        for(int j = 0; j < bins; j++){
                            for(int k=0; k<REGIONS[ri].vIsos.size(); k++){
                                REGIONS[ri].vIsos[k].SKERNAL[i][j][leg_i] = REGIONS[ri].vIsos[k].SKERNAL[i][j][leg_i]*REGIONS[ri].vIsos[k].NumDen;
                            }
                        }
                    }
                }
                /// REGION SKERNAL ///
                for(int leg_i = 0; leg_i < LEGENDRE; leg_i++){
                    for(int i = 0; i < bins; i++){
                        for(int j = 0; j < bins; j++){
                            for(int k=0; k<REGIONS[ri].vIsos.size(); k++){
                                REGIONS[ri].SKERNAL[i][j][leg_i] = REGIONS[ri].SKERNAL[i][j][leg_i] + REGIONS[ri].vIsos[k].SKERNAL[i][j][leg_i];
                            }
                        }
                    }
                }
                /// REGION SCATTER ///
                for(int i = 0; i < bins; i++){
                    for(int j = 0; j < bins; j++){
                        for(int k=0; k<REGIONS[ri].vIsos.size(); k++){
                            REGIONS[ri].SKATTER[i] = REGIONS[ri].SKATTER[i] + REGIONS[ri].vIsos[k].SKERNAL[j][i][0];
                        }
                    }
                }
                /// REGION ABSORB ///
                for(int i = 0; i < bins; i++){
                    REGIONS[ri].ABSORB[i] = REGIONS[ri].TOTAL[i] - REGIONS[ri].SKATTER[i];
                }
                /// Setting up the fission calculations ///
                /*for(unsigned int i = 0; i < REGIONS[ri].vIsos.size(); i ++){
                    REGIONS[ri].fission = REGIONS[ri].fission + REGIONS[ri].vIsos[i].CHI * REGIONS[ri].vIsos[i].NUFISSION * REGIONS[ri].vIsos[i].NumDen;
                }*/
                /// Determining MAS MESH ///
                double max_xs;
                max_xs = REGIONS[ri].TOTAL[0];
                for(unsigned int ii = 0; ii < REGIONS[ri].TOTAL.size(); ii++){
                    if(REGIONS[ri].TOTAL[ii] >= max_xs){
                        max_xs = REGIONS[ri].TOTAL[ii];
                    }
                }
                double min_mew;
                min_mew = mew[0];
                for(int ii = 0; ii < Ordinate; ii++){
                    if(abs(mew[ii]) <= min_mew){
                        min_mew = abs(mew[ii]);
                    }
                }
                REGIONS[ri].max_mesh = (2*(1/max_xs)*(min_mew));
                REGIONS[ri].mesh_size = REGIONS[ri].x_upper - REGIONS[ri].x_lower;
                while(REGIONS[ri].mesh_size > REGIONS[ri].max_mesh){
                    REGIONS[ri].mesh_size = REGIONS[ri].mesh_size/2;
                }
                REGIONS[ri].mesh_points = (REGIONS[ri].x_upper - REGIONS[ri].x_lower)/REGIONS[ri].mesh_size;
                REGIONS[ri].meshpoints.resize(REGIONS[ri].mesh_points);
                /// Source and Resizing ///
                int rn = 0;
                if(Dimensions == 1){
                    /// Resizing Flux and Source Vectors///
                    for(int rmesh = 0; rmesh < REGIONS[ri].mesh_points; rmesh++){
                        REGIONS[ri].meshpoints[rmesh].total_flux.resize(3);
                        REGIONS[ri].meshpoints[rmesh].old_total_flux.resize(3);
                        for(int rj = 0; rj < EDGES; rj++){
                            REGIONS[ri].meshpoints[rmesh].total_flux[rj].resize(1);
                            REGIONS[ri].meshpoints[rmesh].total_flux[rj][0].resize(1);
                            REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0].resize(Ordinate);
                            REGIONS[ri].meshpoints[rmesh].old_total_flux[rj].resize(1);
                            REGIONS[ri].meshpoints[rmesh].old_total_flux[rj][0].resize(1);
                            REGIONS[ri].meshpoints[rmesh].old_total_flux[rj][0][0].resize(Ordinate);
                            for(unsigned int rk = 0; rk < REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0].size(); rk++){
                                    REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0][rk].resize(bins);
                            }
                            for(unsigned int rk = 0; rk < REGIONS[ri].meshpoints[rmesh].old_total_flux[rj][0][0].size(); rk++){
                                    REGIONS[ri].meshpoints[rmesh].old_total_flux[rj][0][0][rk].resize(bins);
                            }
                            for(unsigned int rk = 0; rk < REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0].size(); rk++){
                                for(int eng = 0; eng < bins; eng++){
                                    REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0][rk][eng] = 0;
                                }
                            }
                            for(unsigned int rk = 0; rk < REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0].size(); rk++){
                                for(int eng = 0; eng < bins; eng++){
                                    REGIONS[ri].meshpoints[rmesh].old_total_flux[rj][0][0][rk][eng] = 0;
                                }
                            }
                        }
                        REGIONS[ri].meshpoints[rmesh].flux.resize(3);
                        for(int rj = 0; rj < EDGES; rj++){
                            REGIONS[ri].meshpoints[rmesh].flux[rj].resize(1);
                            REGIONS[ri].meshpoints[rmesh].flux[rj][0].resize(1);
                            REGIONS[ri].meshpoints[rmesh].flux[rj][0][0].resize(Ordinate);
                            for(unsigned int rk = 0; rk < REGIONS[ri].meshpoints[rmesh].flux[rj][0][0].size(); rk++){
                                REGIONS[ri].meshpoints[rmesh].flux[rj][0][0][rk].resize(bins);
                                for(int eng = 0; eng < bins; eng++){
                                    REGIONS[ri].meshpoints[rmesh].flux[rj][0][0][rk][eng] =  0;
                                }
                            }
                        }
                        REGIONS[ri].meshpoints[rmesh].source.resize(1);
                        REGIONS[ri].meshpoints[rmesh].source[0].resize(1);
                        REGIONS[ri].meshpoints[rmesh].source[0][0].resize(Ordinate);
                        for(unsigned int rk = 0; rk < REGIONS[ri].meshpoints[rmesh].source[0][0].size(); rk++){
                            REGIONS[ri].meshpoints[rmesh].source[0][0][rk].resize(bins);
                            for(int mk = 0; mk < bins; mk++){
                                REGIONS[ri].meshpoints[rmesh].source[0][0][rk][mk] = 0;
                            }
                        }

                        REGIONS[ri].meshpoints[rmesh].x_lower = REGIONS[ri].x_lower + rn*REGIONS[ri].mesh_size;
                        REGIONS[ri].meshpoints[rmesh].x_int = REGIONS[ri].x_lower + (rn+0.5)*REGIONS[ri].mesh_size;
                        REGIONS[ri].meshpoints[rmesh].x_upper = REGIONS[ri].x_lower + (rn+1)*REGIONS[ri].mesh_size;
                        rn++;
                    }
                    /// Source ///
                    for(unsigned int rmesh = 0; rmesh < REGIONS[ri].meshpoints.size(); rmesh++){
                        for(int mew1 = 0; mew1 < Ordinate; mew1++){
                            for(int to_eng = 0; to_eng < bins; to_eng++){
                                for(int from_e = 0; from_e < bins; from_e++){
                                    REGIONS[ri].meshpoints[rmesh].source[0][0][mew1][to_eng] = SOURCE/((Ordinate/2) * Wmew[mew1]);
                                }
                            }
                        }
                    }
                }
            }

            cout << "The regions work!" << endl;
            /// Regions are finished ///

            /**This is the problem. I need a way to calculate the power based on a burn up of the fuel inside the reactor. Perhaps just use the fissionable material to start?
            double POWER;
            POWER = ((REGIONS[0].vIsos[0].NumDen + REGIONS[0].vIsos[1].NumDen) * 1.0E24)*(REGIONS[0].x_upper - REGIONS[0].x_lower)/(N_A)*238.0 / 1000 * Burnup * (1/(AVAIL* 365.25*6)) * 1.0E6;
            **/
            /// Building Old flux ///
            for(unsigned int fi = 0; fi < REGIONS.size(); fi++){
                for(unsigned fmesh = 0; fmesh < REGIONS[fi].meshpoints.size(); fmesh ++){
                    for(unsigned int mew1 = 0; mew1 < mew.size(); mew1++){
                        for(int eng = 0; eng < bins; eng++){
                            REGIONS[fi].meshpoints[fmesh].old_total_flux[1][0][0][mew1][eng] = 10000;
                        }
                    }
                }
            }
            int t = 0;
            /// Burn up code ///
            while(t < ITERATIONS){
                /// Criticality Calculations ///
                bool tryagain1 = true;
                while(tryagain1 == true){
                    tryagain1 = false;
                    /// Flux Calculations ///
                    int count = 0;
                    bool tryagain = true;
                    while(tryagain == true){
                        tryagain = false;
                        /// SWEEP THE LEG ///
                        int si = 0;
                        while(si < REGIONS.size()){
                            for(int imesh = 0; imesh < REGIONS[si].mesh_points; imesh++){
                                for(unsigned int mew1 = 0; mew1 < mew.size()/2; mew1++){
                                    for(unsigned int eng = 0; eng < bins; eng++){
                                        REGIONS[si].meshpoints[imesh].flux[1][0][0][mew1][eng] = (REGIONS[si].mesh_size*REGIONS[si].meshpoints[imesh].source[0][0][mew1][eng] + 2*mew[mew1]*REGIONS[si].meshpoints[imesh].flux[0][0][0][mew1][eng]) / (2*mew[mew1] + REGIONS[si].mesh_size*REGIONS[si].TOTAL[eng]);
                                        REGIONS[si].meshpoints[imesh].flux[2][0][0][mew1][eng] = 2*(REGIONS[si].meshpoints[imesh].flux[1][0][0][mew1][eng]) - REGIONS[si].meshpoints[imesh].flux[0][0][0][mew1][eng];
                                        /// Calculating Total Flux ///
                                        for(int edge = 0; edge < EDGES; edge++){
                                            if(count > 0){
                                                REGIONS[si].meshpoints[imesh].total_flux[edge][0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].total_flux[edge][0][0][mew1][eng] + REGIONS[si].meshpoints[imesh].flux[edge][0][0][mew1][eng];
                                            }
                                            else if(count == 0){
                                                REGIONS[si].meshpoints[imesh].total_flux[edge][0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].flux[edge][0][0][mew1][eng];
                                            }
                                        }
                                        /// Moving from one region to the next ///
                                        if(imesh == REGIONS[si].mesh_points-1){
                                            if(si < REGIONS.size()-1){
                                                REGIONS[si+1].meshpoints[0].flux[0][0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].flux[2][0][0][mew1][eng];
                                            }
                                            if(si == REGIONS.size()-1 && eng == 9 && mew1 == mew.size()/2){
                                                break;
                                            }
                                        }
                                        else{
                                            REGIONS[si].meshpoints[imesh+1].flux[0][0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].flux[2][0][0][mew1][eng];
                                        }
                                    }
                                }
                            }
                            si=si+1;
                        }
                        /// REVERSE! REVERSE! ///
                        si = REGIONS.size()-1;
                        while(si > -1){
                            for(int imesh = REGIONS[si].mesh_points-1; imesh > -1; imesh --){
                                for(unsigned int mew1 = mew.size()/2; mew1 < mew.size(); mew1++){
                                    for(unsigned int eng = 0; eng < bins; eng++){
                                        REGIONS[si].meshpoints[imesh].flux[1][0][0][mew1][eng] = (REGIONS[si].mesh_size*REGIONS[si].meshpoints[imesh].source[0][0][mew1][eng] - 2*mew[mew1]*REGIONS[si].meshpoints[imesh].flux[2][0][0][mew1][eng]) / (-2*mew[mew1] + REGIONS[si].mesh_size*REGIONS[si].TOTAL[eng]);
                                        REGIONS[si].meshpoints[imesh].flux[0][0][0][mew1][eng] = 2 * REGIONS[si].meshpoints[imesh].flux[1][0][0][mew1][eng] - REGIONS[si].meshpoints[imesh].flux[2][0][0][mew1][eng];
                                        /// Calculating Total Flux ///
                                        for(int edge = 0; edge < EDGES; edge ++){
                                            if(count > 0){
                                                REGIONS[si].meshpoints[imesh].total_flux[edge][0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].total_flux[edge][0][0][mew1][eng] + REGIONS[si].meshpoints[imesh].flux[edge][0][0][mew1][eng];
                                            }
                                            else if(count == 0){
                                                REGIONS[si].meshpoints[imesh].total_flux[edge][0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].flux[edge][0][0][mew1][eng];
                                            }
                                        }
                                        /// Moving from one region to the next ///
                                        if(imesh == 0){
                                            //cout << si << "    "<< imesh << "  "<< mew1 << "   " <<eng << endl;
                                            if(si > 0){
                                                REGIONS[si-1].meshpoints[REGIONS[si-1].mesh_points-1].flux[2][0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].flux[0][0][0][mew1][eng];
                                            }
                                            if(si == 0 && mew1 == mew.size()-1 && eng ==9){
                                                break;
                                            }
                                        }
                                        else{
                                            REGIONS[si].meshpoints[imesh-1].flux[2][0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].flux[0][0][0][mew1][eng];
                                        }
                                    }
                                }
                            }
                            si=si-1;
                        }
                        /// SCATTER!! ///;
                        /// Zeroing Scattering Source ///
                        for(unsigned int si = 0; si < REGIONS.size(); si++){
                            for(int imesh = 0;imesh < REGIONS[si].mesh_points; imesh++){
                                for(unsigned int mew1 = 0; mew1 < mew.size(); mew1++){
                                    for(unsigned int eng = 0; eng < bins; eng++){
                                        REGIONS[si].meshpoints[imesh].source[0][0][mew1][eng] = 0;
                                    }
                                }
                            }
                        }
                        /// Calculating New Scattering Source ///
                        for(unsigned int si = 0; si < REGIONS.size(); si++){
                            for(int imesh = 0;imesh < REGIONS[si].mesh_points; imesh++){
                                for(unsigned int mew1 = 0; mew1 < mew.size(); mew1++){
                                    for(int eng = 0; eng < bins; eng++){
                                        for(int l = 0; l < LEGENDRE; l++){
                                            for(int i = 0; i < bins; i++){
                                                for(unsigned int mew2 = 0; mew2 < mew.size(); mew2++){
                                                    REGIONS[si].meshpoints[imesh].source[0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].source[0][0][mew1][eng] + (0.5)*((2*l+1)*tr1::legendre(l, mew[mew1]))*(REGIONS[si].SKERNAL[eng][i][l])*(Wmew[mew2]*tr1::legendre(l, mew[mew2]) * REGIONS[si].meshpoints[imesh].flux[1][0][0][mew2][i]);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        /// Convergence Test ///
                        if(count > 0){
                            for(unsigned int si = 0; si < REGIONS.size(); si++){
                                for(unsigned int imesh = 0; imesh < REGIONS[si].mesh_points-1; imesh ++){
                                    for(unsigned int mew1 = 0; mew1 < mew.size(); mew1++){
                                        for(unsigned int eng = 0; eng < bins; eng++){
                                            if(REGIONS[si].meshpoints[imesh].flux[1][0][0][mew1][eng] > CONVERGENCE * REGIONS[si].meshpoints[imesh].total_flux[1][0][0][mew1][eng]){
                                                tryagain = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else if(count == 0){
                            tryagain = true;
                        }
                        count = count + 1;
                    }
                    /// nufission ///
                    vector<double> nufission_t;
                    nufission_t.resize(bins);
                    for(unsigned int ci = 0; ci < REGIONS.size(); ci++){
                        for(unsigned int cIsos = 0; cIsos < REGIONS[ci].vIsos.size(); cIsos++){
                            for(int eng = 0; eng < bins; eng++){
                                nufission_t[eng] = nufission_t[eng] + REGIONS[ci].vIsos[cIsos].NumDen * REGIONS[ci].vIsos[cIsos].NUFISSION(eng);
                            }
                        }
                    }
                    /// Old/New nufiss ///
                    double nufiss1 = 0;
                    double nufiss2 = 0;
                    for(unsigned int ci = 0; ci < REGIONS.size(); ci++){
                        for(unsigned int cIsos = 0; cIsos < REGIONS[ci].vIsos.size(); cIsos++){
                            for(int eng = 0; eng < bins; eng++){
                                for(unsigned int cmesh = 0; cmesh < REGIONS[ci].meshpoints.size(); cmesh++){
                                    for(unsigned int mew1 = 0; mew1 < mew.size(); mew1++){
                                        nufiss1 = nufiss1 + nufission_t[eng] * REGIONS[ci].meshpoints[cmesh].total_flux[1][0][0][mew1][eng] * Wmew[mew1]/2;
                                        nufiss2 = nufiss2 + nufission_t[eng] * REGIONS[ci].meshpoints[cmesh].old_total_flux[1][0][0][mew1][eng] * Wmew[mew1]/2;
                                    }
                                }
                            }
                        }
                    }

                    /// Criticality Convergence Test ///
                    crit2 = crit1 * (nufiss1/nufiss2);
                    if(abs(crit2-crit1) > CONVERGENCE*crit1){
                        tryagain1 = true;
                    }
                    /// FLux Convergence Test ///
                    for(unsigned int si = 0; si < REGIONS.size(); si++){
                        for(unsigned int imesh = 0; imesh < REGIONS[si].mesh_points-1; imesh ++){
                            for(unsigned int mew1 = 0; mew1 < mew.size(); mew1++){
                                for(unsigned int eng = 0; eng < bins; eng++){
                                    if(abs(REGIONS[si].meshpoints[imesh].total_flux[1][0][0][mew1][eng] - REGIONS[si].meshpoints[imesh].old_total_flux[1][0][0][mew1][eng]) > (CONVERGENCE * REGIONS[si].meshpoints[imesh].old_total_flux[1][0][0][mew1][eng])){
                                        tryagain1 = true;
                                    }
                                }
                            }
                        }
                    }
                    crit1 = crit2;
                    nufission_t.clear();
                    /// Zeroing Final Flux ///
                    for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                        for(unsigned int omesh = 0; omesh<REGIONS[oi].meshpoints.size(); omesh++){
                            REGIONS[oi].meshpoints[omesh].final_flux_e.resize(bins);
                            for(unsigned int eng = 0; eng < bins; eng ++){
                                REGIONS[oi].meshpoints[omesh].final_flux_e[eng]=0;
                            }
                        }
                    }
                    /// Calculating Final Flux by Energy ///
                    for(unsigned int eng = 0; eng < bins; eng++){
                        for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                            for(unsigned int omesh = 0; omesh<REGIONS[oi].meshpoints.size(); omesh++){
                                for(unsigned int mew1 = 0; mew1 < mew.size(); mew1 ++){
                                    REGIONS[oi].meshpoints[omesh].final_flux_e[eng] = REGIONS[oi].meshpoints[omesh].final_flux_e[eng] + Wmew[mew1]*REGIONS[oi].meshpoints[omesh].total_flux[1][0][0][mew1][eng] / 2.0;
                                }
                            }
                        }
                    }
                    /// Zeroing Final Flux
                    for(unsigned int eng = 0; eng < bins; eng++){
                        for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                            for(unsigned int omesh = 0; omesh<REGIONS[oi].meshpoints.size(); omesh++){
                                    REGIONS[oi].meshpoints[omesh].final_flux = 0;
                            }
                        }
                    }

                    /// Calculating Final Flux
                    for(unsigned int eng = 0; eng < bins; eng++){
                        for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                            for(unsigned int omesh = 0; omesh<REGIONS[oi].meshpoints.size(); omesh++){
                                    REGIONS[oi].meshpoints[omesh].final_flux = REGIONS[oi].meshpoints[omesh].final_flux + REGIONS[oi].meshpoints[omesh].final_flux_e[eng];
                            }
                        }
                    }
                    /// Updating Old Flux
                    for(unsigned int fi = 0; fi < REGIONS.size(); fi++){
                        for(unsigned fmesh = 0; fmesh < REGIONS[fi].meshpoints.size(); fmesh ++){
                            for(unsigned int mew1 = 0; mew1 < mew.size(); mew1++){
                                for(int eng = 0; eng < bins; eng++){
                                    REGIONS[fi].meshpoints[fmesh].old_total_flux[1][0][0][mew1][eng] = REGIONS[fi].meshpoints[fmesh].total_flux[1][0][0][mew1][eng];
                                }
                            }
                        }
                    }
                    /// Zeroing Source
                    for(unsigned int ri = 0; ri < REGIONS.size();ri++){
                        for(unsigned int rmesh = 0; rmesh < REGIONS[ri].meshpoints.size(); rmesh++){
                            for(int mew1 = 0; mew1 < Ordinate; mew1++){
                                for(int to_eng = 0; to_eng < bins; to_eng++){
                                    for(int from_e = 0; from_e < bins; from_e++){
                                        REGIONS[ri].meshpoints[rmesh].source[0][0][mew1][to_eng] = 0;
                                    }
                                }
                            }
                        }
                    }
                    /// Creating new source ///
                    for(unsigned int ri = 0; ri < REGIONS.size();ri++){
                        for(unsigned int rmesh = 0; rmesh < REGIONS[ri].meshpoints.size(); rmesh++){
                            for(int mew1 = 0; mew1 < Ordinate; mew1++){
                                for(int to_e = 0; to_e< bins; to_e++){
                                    for(int from_e = 0; from_e < bins; from_e++){
                                        for(int rIsos = 0; rIsos < REGIONS[ri].vIsos.size(); rIsos++){
                                            REGIONS[ri].meshpoints[rmesh].source[0][0][mew1][to_e] = REGIONS[ri].meshpoints[rmesh].source[0][0][mew1][to_e] + (REGIONS[ri].vIsos[rIsos].CHI(to_e)/crit1)*(REGIONS[ri].meshpoints[rmesh].final_flux_e[from_e]*REGIONS[ri].vIsos[rIsos].NUFISSION[from_e]*REGIONS[ri].vIsos[rIsos].NumDen);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /// Total Flux if convergence fails ///
                    if(tryagain1==true){
                        for(unsigned int ri = 0; ri < REGIONS.size();ri++){
                            for(int rmesh = 0; rmesh < REGIONS[ri].mesh_points; rmesh++){
                                for(int rj = 0; rj < EDGES; rj++){;
                                    for(unsigned int rk = 0; rk < Ordinate; rk++){
                                        for(int eng = 0; eng < bins; eng++){
                                            REGIONS[ri].meshpoints[rmesh].flux[rj][0][0][rk][eng] =  REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0][rk][eng];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                /// BURN UP CODE ///
                /// ONE GROUP COLLAPSE ///
                for(unsigned int ri = 0; ri < REGIONS.size();ri++){
                    REGIONS[ri].FINAL_FLUX.resize(bins);
                    for(int eng = 0 ; eng < bins; eng++){
                        REGIONS[ri].FINAL_FLUX[eng] = 0;
                    }
                    REGIONS[ri].FINAL_ff = 0;
                    for(unsigned int rIsos = 0; rIsos < REGIONS[ri].vIsos.size(); rIsos++){
                        REGIONS[ri].vIsos[rIsos].ABSORB_C = 0;
                        REGIONS[ri].vIsos[rIsos].NUFISSION_C = 0;
                    }
                    for(int eng = 0; eng < bins; eng++){
                        for(int rmesh = 0; rmesh < REGIONS[ri].mesh_points; rmesh++){
                            REGIONS[ri].FINAL_FLUX[eng] =+ REGIONS[ri].meshpoints[rmesh].final_flux_e[eng];
                        }
                        REGIONS[ri].FINAL_ff = REGIONS[ri].FINAL_ff + REGIONS[ri].FINAL_FLUX[eng];
                        for(unsigned int rIsos = 0; rIsos < REGIONS[ri].vIsos.size(); rIsos++){
                            REGIONS[ri].vIsos[rIsos].ABSORB_C = REGIONS[ri].vIsos[rIsos].ABSORB_C + REGIONS[ri].vIsos[rIsos].ABSORB[eng]*REGIONS[ri].FINAL_FLUX[eng];
                            REGIONS[ri].vIsos[rIsos].NUFISSION_C = REGIONS[ri].vIsos[rIsos].NUFISSION_C + REGIONS[ri].vIsos[rIsos].NUFISSION(eng)*REGIONS[ri].FINAL_FLUX[eng];
                        }
                    }
                    for(unsigned int rIsos = 0; rIsos < REGIONS[ri].vIsos.size(); rIsos++){
                        REGIONS[ri].vIsos[rIsos].ABSORB_C = REGIONS[ri].vIsos[rIsos].ABSORB_C/REGIONS[ri].FINAL_ff;
                        REGIONS[ri].vIsos[rIsos].NUFISSION_C = REGIONS[ri].vIsos[rIsos].NUFISSION_C/REGIONS[ri].FINAL_ff;
                    }
                }
                /// Fission_sum for all regions ///
                for(unsigned int ri = 0; ri < REGIONS.size(); ri++){
                    REGIONS[ri].fission_sum = 0;
                    for(unsigned int rIsos = 0; rIsos < REGIONS[ri].vIsos.size(); rIsos++){
                        REGIONS[ri].fission_sum = REGIONS[ri].fission_sum + ((3.2E-11)*(REGIONS[ri].vIsos[rIsos].NUFISSION_C/REGIONS[ri].vIsos[rIsos].ABSORB_C)) * REGIONS[ri].vIsos[rIsos].NumDen * REGIONS[ri].vIsos[rIsos].ABSORB_C;
                    }
                }
                /// Power Rescaling for all regions ///
                for(unsigned int ri = 0; ri < REGIONS.size(); ri++){
                    double fission_sum = 0;
                    for(unsigned int ti = 0; ti < REGIONS.size(); ti++){
                        if(ti != ri){
                            fission_sum =+ REGIONS[ti].fission_sum * REGIONS[ti].FINAL_ff * (REGIONS[ti].x_upper - REGIONS[ti].x_lower);
                            cout << REGIONS[ti].fission_sum << "    " << REGIONS[ti].FINAL_ff << endl;
                        }
                    }
                    double power_gap = POWER -fission_sum;
                    double power = POWER
                    REGIONS[ri].power_flux = power_gap / (REGIONS[ri].fission_sum * (REGIONS[ri].x_upper - REGIONS[ri].x_lower));
                    if( REGIONS.size() == 1){
                        REGIONS[ri].power_flux = power / (REGIONS[ri].fission_sum * (REGIONS[ri].x_upper - REGIONS[ri].x_lower));
                    }
                }
                /// BATEMAN EQUATIONS ///
                for(unsigned int ti = 0; ti < REGIONS.size(); ti++){
                    for(unsigned int tIsos = 0; tIsos < REGIONS[ti].vIsos.size(); tIsos++){
                        REGIONS[ti].vIsos[tIsos].CHANGE = 0;
                    }
                    REGIONS[ti].vIsos[0].CHANGE = REGIONS[ti].vIsos[0].CHANGE - (REGIONS[ti].vIsos[0].NUFISSION_C + REGIONS[ti].vIsos[0].ABSORB_C)*REGIONS[ti].power_flux*REGIONS[ti].vIsos[0].NumDen*Convert;
                    REGIONS[ti].vIsos[1].CHANGE = REGIONS[ti].vIsos[1].CHANGE - (REGIONS[ti].vIsos[1].NUFISSION_C + REGIONS[ti].vIsos[1].ABSORB_C)*REGIONS[ti].power_flux*REGIONS[ti].vIsos[1].NumDen*Convert;
                    REGIONS[ti].vIsos[2].CHANGE = REGIONS[ti].vIsos[2].CHANGE + (REGIONS[ti].vIsos[0].ABSORB_C)*REGIONS[ti].power_flux*REGIONS[ti].vIsos[0].NumDen*Convert;
                    REGIONS[ti].vIsos[3].CHANGE = REGIONS[ti].vIsos[3].CHANGE - ((REGIONS[ti].vIsos[3].NUFISSION_C + REGIONS[ti].vIsos[3].ABSORB_C)*REGIONS[ti].vIsos[3].NumDen - REGIONS[ti].vIsos[1].ABSORB_C * REGIONS[ti].vIsos[1].NumDen)*REGIONS[ti].power_flux*Convert;
                    REGIONS[ti].vIsos[4].CHANGE = REGIONS[ti].vIsos[4].CHANGE + (REGIONS[ti].vIsos[3].ABSORB_C)*REGIONS[ti].power_flux*REGIONS[ti].vIsos[3].NumDen*Convert;
                    REGIONS[ti].vIsos[5].CHANGE = REGIONS[ti].vIsos[5].CHANGE + (1.87*(REGIONS[ti].vIsos[0].NUFISSION_C*REGIONS[ti].vIsos[0].NumDen +REGIONS[ti].vIsos[1].NUFISSION_C*REGIONS[ti].vIsos[1].NumDen+REGIONS[ti].vIsos[3].NUFISSION_C*REGIONS[ti].vIsos[3].NumDen)+REGIONS[ti].vIsos[6].ABSORB_C*REGIONS[ti].vIsos[6].NumDen+REGIONS[ti].vIsos[7].ABSORB_C*REGIONS[ti].vIsos[7].NumDen)*REGIONS[ti].power_flux*Convert;
                    REGIONS[ti].vIsos[6].CHANGE = REGIONS[ti].vIsos[6].CHANGE + (0.065*(REGIONS[ti].vIsos[0].NUFISSION_C * REGIONS[ti].vIsos[0].NumDen + REGIONS[ti].vIsos[1].NUFISSION_C * REGIONS[ti].vIsos[1].NumDen + REGIONS[ti].vIsos[3].NUFISSION_C*REGIONS[ti].vIsos[3].NumDen)-(REGIONS[ti].vIsos[6].ABSORB_C + 3.04E-5)*REGIONS[ti].vIsos[6].NumDen)*REGIONS[ti].power_flux*Convert;
                    REGIONS[ti].vIsos[7].CHANGE = REGIONS[ti].vIsos[7].CHANGE + (0.065*(REGIONS[ti].vIsos[0].NUFISSION_C * REGIONS[ti].vIsos[0].NumDen + REGIONS[ti].vIsos[1].NUFISSION_C * REGIONS[ti].vIsos[1].NumDen + REGIONS[ti].vIsos[3].NUFISSION_C*REGIONS[ti].vIsos[3].NumDen)-REGIONS[ti].vIsos[7].ABSORB_C*REGIONS[ti].vIsos[7].NumDen)*REGIONS[ti].power_flux*Convert;

                    for(unsigned int tIsos = 0; tIsos < REGIONS[ti].vIsos.size(); tIsos++){
                        REGIONS[ti].vIsos[tIsos].NumDen = REGIONS[ti].vIsos[tIsos].NumDen + TIME_STEP*REGIONS[ti].vIsos[tIsos].CHANGE;
                    }
                }
                if(t%10 == 0){
                    for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                        Isos4 << REGIONS[oi].vIsos[3].NumDen << "   ";
                        Isos6 << REGIONS[oi].fission_sum*REGIONS[oi].power_flux*(REGIONS[oi].x_upper - REGIONS[oi].x_lower)<< "     ";
                        }
                    Isos4 << endl;
                    Isos6 << endl;
                    Isos5 << t*TIME_STEP /86400.0 << "   "<<crit1 << endl;
                }
                t=t+1;
            }
            double leak_left = 0;
            double leak_right = 0;
            /// LEAK RIGHT ///
            for(unsigned int mew1 = 0; mew1 < mew.size()/2; mew1++){
                for(unsigned int eng = 0; eng < bins; eng++){
                    leak_right = leak_right + (0.5)*Wmew[mew1]*mew[mew1]*REGIONS[REGIONS.size()-1].meshpoints[REGIONS[REGIONS.size()-1].mesh_points-1].total_flux[2][0][0][mew1][eng];
                }
            }
            /// LEAK LEFT ///
            for(unsigned int mew1 = mew.size()/2; mew1 < mew.size(); mew1++){
                for(unsigned int eng = 0; eng < bins; eng++){
                    leak_left = leak_left + (0.5)* Wmew[mew1]*mew[mew1]*REGIONS[0].meshpoints[0].total_flux[0][0][0][mew1][eng];
                }
            }
            /// Outputting the Flux Distributions ///
            for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                for(unsigned int omesh = 0; omesh < REGIONS[oi].meshpoints.size(); omesh++){
                    for(unsigned int eng =0; eng < bins; eng++){
                        Isos3.width(6);
                        Isos3 << REGIONS[oi].meshpoints[omesh].x_int << "   " << REGIONS[oi].meshpoints[omesh].final_flux_e[eng];
                        Isos3.width(10);
                        Isos3 << "  ";
                    }
                    Isos3 << endl;
                }
            }
            Isos3 << endl;
            Isos4 << endl;
            Isos5 << endl;
            Isos6 << endl;
            /// Zeroing the vectors of the regions ///
            /*for(int ci = 0; ci < REGIONS.size(); ci++){
                REGIONS[ci].meshpoints.clear();
                REGIONS[ci].mesh_points=0;
                REGIONS[ci].mesh_size=0;
                REGIONS[ci].max_mesh=0;
                REGIONS[ci].SKERNAL.clear();
                REGIONS[ci].TOTAL.clear();
                REGIONS[ci].ABSORB.clear();
                REGIONS[ci].SKATTER.clear();
                REGIONS[ci].vIsos.clear();
            }*/


}


