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
#define bins 10.0
#define PI 3.14159265359
#define N_A 6.0221415 * pow(10.0, 23)
#define cm2_to_barns pow(10.0, 24)
#define barn_to_cm2 pow(10.0, -24)
#define sec_per_day 24.0 * 3600.0
#define ffactors 6
#define Dimensions 1
#define CONVERGENCE 0.005
#define LEGENDRE 3
#define REGIONSETS 3
#define Ordinate 8
#define SOURCE pow(10.0, 10.0)
#define REGION_SIZE 2.0
#define MAX_BOUND 20.0
#define BENCHMARK 1

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
    vector<vector<vector<vector<double> > > > source;
    vector<vector<vector<vector<vector<double > > > > > flux;
    vector<vector<vector<vector<vector<double > > > > > total_flux;
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
        vector<double> SKATTER;
    protected:
    private:
};


int main()
{
    int total_count = 0;
    vector<Region> REGIONS;
    ofstream Isos3;
    ofstream Isos4;
    ofstream Isos5;
    ofstream Isos6;
    Isos3.open("Flux.txt");
    Isos4.open("REGIONS.txt");
    Isos5.open("absorb.txt");
    Isos6.open("ordinateF.txt");
    while(total_count <2){

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
            vector<double>  atea;
            atea.resize(Ordinate);
            vector<double>  cacee;
            cacee.resize(Ordinate);
            vector<double>  Wmew;
            Wmew.resize(Ordinate);
            vector<double>  Watea;
            Watea.resize(Ordinate);
            vector<double>  Wcacee;
            Wcacee.resize(Ordinate);
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
            vector<double> nofr_minmax;
            if(total_count == 0){
                srand(time(0));
                double r_num;
                r_num = floor(rand()/double(RAND_MAX)*(MAX_BOUND-REGION_SIZE)*100)/100;
                bool test_region;
                vector<double> nofr;
                nofr.push_back(r_num);
                while(nofr.size() < REGIONSETS){
                    test_region = true;
                    r_num = floor(rand()/double(RAND_MAX)*(MAX_BOUND-REGION_SIZE)*100)/100;
                    for(unsigned int n = 0; n < nofr.size(); n++){
                        if( REGION_SIZE >= abs(nofr[n]-r_num)){
                            test_region = false;
                        }
                    }
                    if(test_region == true){
                        nofr.push_back(r_num);
                    }
                }
                int removed;
                while(nofr.size() > 0){
                    double min_r = nofr[0];
                    for(unsigned int n = 0; n < nofr.size(); n++){
                        if(nofr[n] <= min_r){
                            min_r = nofr[n];
                            removed = n;
                        }
                    }
                    nofr_minmax.push_back(min_r);
                    nofr.erase(nofr.begin()+removed);
                }
                Region region;
            /// No Benchmark ///
                if(BENCHMARK == 0){
                    if(nofr_minmax[0] == 0){
                        for(int iv = 0; iv <6; iv ++){
                            REGIONS.push_back(region);
                        }
                        REGIONS[0].x_lower = 0;
                        REGIONS[0].x_upper = REGION_SIZE;
                        REGIONS[0].r_type = 2;
                        REGIONS[1].x_lower = nofr_minmax[0];
                        REGIONS[1].x_upper = nofr_minmax[1];
                        REGIONS[1].r_type = 1;
                        REGIONS[2].x_lower = nofr_minmax[1];
                        REGIONS[2].x_upper = nofr_minmax[1]+REGION_SIZE;
                        REGIONS[2].r_type = 2;
                        REGIONS[3].x_lower = nofr_minmax[1]+REGION_SIZE;
                        REGIONS[3].x_upper = nofr_minmax[2];
                        REGIONS[3].r_type = 1;
                        REGIONS[4].x_lower = nofr_minmax[2];
                        REGIONS[4].x_upper = nofr_minmax[2]+ REGION_SIZE;
                        REGIONS[4].r_type = 2;
                        REGIONS[5].x_lower = nofr_minmax[2]+REGION_SIZE;
                        REGIONS[5].x_upper = MAX_BOUND;
                        REGIONS[5].r_type = 1;

                    }
                    else if(nofr_minmax[nofr_minmax.size()-1] == MAX_BOUND - REGION_SIZE){
                        for(int iv = 0; iv <6; iv ++){
                            REGIONS.push_back(region);
                        }
                        REGIONS[0].x_lower = 0;
                        REGIONS[0].x_upper = nofr_minmax[0];
                        REGIONS[0].r_type = 1;
                        REGIONS[1].x_lower = nofr_minmax[0];
                        REGIONS[1].x_upper = nofr_minmax[0]+REGION_SIZE;
                        REGIONS[1].r_type = 2;
                        REGIONS[2].x_lower = nofr_minmax[0]+REGION_SIZE;
                        REGIONS[2].x_upper = nofr_minmax[1];
                        REGIONS[2].r_type = 1;
                        REGIONS[3].x_lower = nofr_minmax[1];
                        REGIONS[3].x_upper = nofr_minmax[1]+REGION_SIZE;
                        REGIONS[3].r_type = 2;
                        REGIONS[4].x_lower = nofr_minmax[1]+REGION_SIZE;
                        REGIONS[4].x_upper = nofr_minmax[2];
                        REGIONS[4].r_type = 1;
                        REGIONS[5].x_lower = nofr_minmax[2];
                        REGIONS[5].x_upper = MAX_BOUND;
                        REGIONS[5].r_type = 2;
                    }
                    else{
                        for(int iv = 0; iv <7; iv ++){
                            REGIONS.push_back(region);
                        }
                        REGIONS[0].x_lower = 0;
                        REGIONS[0].x_upper = nofr_minmax[0];
                        REGIONS[0].r_type = 1;
                        REGIONS[1].x_lower = nofr_minmax[0];
                        REGIONS[1].x_upper = nofr_minmax[0]+REGION_SIZE;
                        REGIONS[1].r_type = 2;
                        REGIONS[2].x_lower = nofr_minmax[0]+REGION_SIZE;
                        REGIONS[2].x_upper = nofr_minmax[1];
                        REGIONS[2].r_type = 1;
                        REGIONS[3].x_lower = nofr_minmax[1];
                        REGIONS[3].x_upper = nofr_minmax[1]+REGION_SIZE;
                        REGIONS[3].r_type = 2;
                        REGIONS[4].x_lower = nofr_minmax[1]+REGION_SIZE;
                        REGIONS[4].x_upper = nofr_minmax[2];
                        REGIONS[4].r_type = 1;
                        REGIONS[5].x_lower = nofr_minmax[2];
                        REGIONS[5].x_upper = nofr_minmax[2]+REGION_SIZE;
                        REGIONS[5].r_type = 2;
                        REGIONS[6].x_lower = nofr_minmax[2]+REGION_SIZE;
                        REGIONS[6].x_upper = MAX_BOUND;
                        REGIONS[6].r_type = 1;
                    }
                }
            /// Setting up the Regions for Benchmark ///
                if(BENCHMARK == 1){
                    if(total_count == 0){
                        for(int iv = 0; iv <7; iv ++){
                            REGIONS.push_back(region);
                        }
                    }
                    REGIONS[0].x_lower = 0;
                    REGIONS[0].x_upper = 3.5;
                    REGIONS[0].r_type = 1;
                    REGIONS[1].x_lower = 3.5;
                    REGIONS[1].x_upper = 5.5;
                    REGIONS[1].r_type = 2;
                    REGIONS[2].x_lower = 5.5;
                    REGIONS[2].x_upper = 9;
                    REGIONS[2].r_type = 1;
                    REGIONS[3].x_lower = 9;
                    REGIONS[3].x_upper = 11;
                    REGIONS[3].r_type = 2;
                    REGIONS[4].x_lower = 11;
                    REGIONS[4].x_upper = 14.5;
                    REGIONS[4].r_type = 1;
                    REGIONS[5].x_lower = 14.5;
                    REGIONS[5].x_upper = 16.5;
                    REGIONS[5].r_type = 2;
                    REGIONS[6].x_lower = 16.5;
                    REGIONS[6].x_upper = 20;
                    REGIONS[6].r_type = 1;
                    }
            }
            if(total_count == 1){
                for(int k =0; k < REGIONS.size(); k++){
                    if(REGIONS[k].r_type == 1){
                        REGIONS[k].r_type = 0;
                    }
                    else if(REGIONS[k].r_type == 2){
                        REGIONS[k].r_type = 10;
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
                double dens;
                while(Isos >> xx){
                    nuclide.NAME = xx;
                    Isos >> dens;
                    nuclide.NumDen = dens*1E-24;
                    REGIONS[ri].vIsos.push_back(nuclide);
                }
                Isos.close();
                // Set up sizes - Region KERNAL. (HEIGHT x WIDTH)
                REGIONS[ri].SKERNAL.resize(bins);
                for (int i = 0; i < bins; ++i){
                    REGIONS[ri].SKERNAL[i].resize(bins);
                    for (int j = 0; j < bins; ++j){
                        REGIONS[ri].SKERNAL[i][j].resize(MaxLN);
                    }
                }
                // Set up sizes - Materials KERNAL. (HEIGHT x WIDTH)
                for(unsigned int k = 0; k < REGIONS[ri].vIsos.size(); k++){
                    REGIONS[ri].vIsos[k].SKERNAL.resize(bins);
                    for (int i = 0; i < bins; ++i){
                        REGIONS[ri].vIsos[k].SKERNAL[i].resize(bins);
                        for (int j = 0; j < bins; ++j){
                            REGIONS[ri].vIsos[k].SKERNAL[i][j].resize(MaxLN);
                        }
                    }
                }

                ///Inputs from cross section files///
                int vIsos_size;
                vIsos_size = REGIONS[ri].vIsos.size();
                ifstream Isos2;
                for(int i = 0; i < vIsos_size; i++){
                    string result;
                    string s1 = "./";
                    string s2 = ".xs";
                    stringstream isotope_num;
                    isotope_num << REGIONS[ri].vIsos[i].NAME;
                    result = s1+isotope_num.str()+s2;
                    Isos2.open(result.c_str());
                    if(!Isos2){
                        cout << "failures!" << endl;
                        cout << result.c_str()<< endl;
                    }
                    string TEXT;
                    REGIONS[ri].vIsos[i].TOTAL.resize(bins);
                    REGIONS[ri].vIsos[i].CHI.resize(bins);
                    REGIONS[ri].vIsos[i].NUFISSION.resize(bins);
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
                       /* FFACTOR INPUT
                       while(Isos2 >> bin){
                            for(int l = 0; l < ffactors; l++){
                                Isos2 >> xx;
                                REGIONS[ri].vIsos[i].FFACTOR[10-bin][l] = xx;
                            }
                        }
                        Isos2.clear();
                        Isos2 >> TEXT;*/
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
                /// REGION TOTAL ///
                for(int i = 0;  i < bins; i++){
                    for(int j =0; j < vIsos_size; j++){
                        sigma_t_sum[i] = sigma_t_sum(i) + REGIONS[ri].vIsos[j].NumDen * REGIONS[ri].vIsos[j].TOTAL(i);
                        REGIONS[ri].TOTAL[i] = REGIONS[ri].TOTAL[i] + REGIONS[ri].vIsos[j].NumDen * REGIONS[ri].vIsos[j].TOTAL(i);
                    }
                }
                /// vIsos SKERNAL ///
                for(int leg_i = 0; leg_i < LEGENDRE; leg_i++){
                    for(int i = 0; i < bins; i++){
                        for(int j = 0; j < bins; j++){
                            for(int k=0; k<vIsos_size; k++){
                                REGIONS[ri].vIsos[k].SKERNAL[i][j][leg_i] = REGIONS[ri].vIsos[k].SKERNAL[i][j][leg_i]*REGIONS[ri].vIsos[k].NumDen;
                            }
                        }
                    }
                }
                /// REGION SKERNAL ///
                for(int leg_i = 0; leg_i < LEGENDRE; leg_i++){
                    for(int i = 0; i < bins; i++){
                        for(int j = 0; j < bins; j++){
                            for(int k=0; k<vIsos_size; k++){
                                REGIONS[ri].SKERNAL[i][j][leg_i] = REGIONS[ri].SKERNAL[i][j][leg_i] + REGIONS[ri].vIsos[k].SKERNAL[i][j][leg_i];
                            }
                        }
                    }
                }
                /// REGION SCATTER ///
                for(int i = 0; i < bins; i++){
                    for(int j = 0; j < bins; j++){
                        for(int k=0; k<vIsos_size; k++){
                            REGIONS[ri].SKATTER[i] = REGIONS[ri].SKATTER[i] + REGIONS[ri].vIsos[k].SKERNAL[j][i][0];
                        }
                    }
                }
                /// REGION ABSORB ///
                for(int i = 0; i < bins; i++){
                    REGIONS[ri].ABSORB[i] = REGIONS[ri].TOTAL[i] - REGIONS[ri].SKATTER[i];
                }

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
                cout << REGIONS[ri].max_mesh << endl;
                REGIONS[ri].mesh_size = REGIONS[ri].x_upper - REGIONS[ri].x_lower;
                while(REGIONS[ri].mesh_size > REGIONS[ri].max_mesh){
                    REGIONS[ri].mesh_size = REGIONS[ri].mesh_size/2;
                }
                REGIONS[ri].mesh_points = (REGIONS[ri].x_upper - REGIONS[ri].x_lower)/REGIONS[ri].mesh_size;
                REGIONS[ri].meshpoints.resize(REGIONS[ri].mesh_points);
                /// Source and Resizing ///
                int rn = 0;
                if(Dimensions == 1){
                    /// Resizing ///
                    for(int rmesh = 0; rmesh < REGIONS[ri].mesh_points; rmesh++){
                        REGIONS[ri].meshpoints[rmesh].total_flux.resize(3);
                        for(int rj = 0; rj < EDGES; rj++){
                            REGIONS[ri].meshpoints[rmesh].total_flux[rj].resize(1);
                            REGIONS[ri].meshpoints[rmesh].total_flux[rj][0].resize(1);
                            REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0].resize(Ordinate);
                            for(unsigned int rk = 0; rk < REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0].size(); rk++){
                                    REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0][rk].resize(bins);
                            }
                            for(unsigned int rk = 0; rk < REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0].size(); rk++){
                                for(int eng = 0; eng < bins; eng++){
                                    REGIONS[ri].meshpoints[rmesh].total_flux[rj][0][0][rk][eng] = 0;
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
                                    REGIONS[ri].meshpoints[rmesh].flux[rj][0][0][rk][eng] = 0;
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
                    if(REGIONS[ri].r_type == 2){
                        for(unsigned int rmesh = 0; rmesh < REGIONS[ri].meshpoints.size(); rmesh++){
                            for(int mew1 = 0; mew1 < Ordinate; mew1++){
                                REGIONS[ri].meshpoints[rmesh].source[0][0][mew1][0] = SOURCE/(1.)*Wmew[mew1];
                            }
                        }
                    }
                    if(REGIONS[ri].r_type == 10){
                        for(unsigned int rmesh = 0; rmesh < REGIONS[ri].meshpoints.size(); rmesh++){
                            for(int mew1 = 0; mew1 < Ordinate; mew1++){
                                REGIONS[ri].meshpoints[rmesh].source[0][0][mew1][0] = SOURCE/(1.)*Wmew[mew1];
                            }
                        }
                    }
                }
            }
            cout << "The regions work!" << endl;

            /// Regions are finished ///

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
                for(int si = 0; si < REGIONS.size(); si++){
                    for(int imesh = 0;imesh < REGIONS[si].mesh_points; imesh++){
                        for(unsigned int mew1 = 0; mew1 < mew.size(); mew1++){
                            for(unsigned int eng = 0; eng < bins; eng++){
                                REGIONS[si].meshpoints[imesh].source[0][0][mew1][eng] = 0;
                            }
                        }
                    }
                }
                /// Calculating New Scattering Source ///
                for(int si = 0; si < REGIONS.size(); si++){
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
                    for(int si = 0; si < REGIONS.size(); si++){
                        for(unsigned int imesh = 0; imesh < REGIONS[si].mesh_points-1; imesh ++){
                            for(unsigned int mew1 = 0; mew1 < mew.size(); mew1++){
                                for(unsigned int eng = 0; eng < bins; eng++){
                                    if(REGIONS[si].meshpoints[imesh].flux[1][0][0][mew1][eng] > CONVERGENCE * REGIONS[si].meshpoints[imesh].total_flux[1][0][0][mew1][eng]){
                                        tryagain = true;
                                        break;
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
                cout << count << endl;
            };
            cout << count << endl;
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
            cout << -leak_left/(SOURCE*3*REGION_SIZE) << "  " << leak_right/(SOURCE*3*REGION_SIZE) << endl;
            /// Calculating the Final flux by energy group ///
            for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                for(unsigned int omesh =0; omesh<REGIONS[oi].meshpoints.size(); omesh++){
                    for(unsigned int eng = 0; eng < bins; eng ++){
                        REGIONS[oi].meshpoints[omesh].final_flux_e.resize(bins);
                    }
                }
            }
            for(unsigned int eng = 0; eng < bins; eng++){
                for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                    for(unsigned int omesh =0; omesh<REGIONS[oi].meshpoints.size(); omesh++){
                        for(unsigned int mew1 = 0; mew1 < mew.size(); mew1 ++){
                            REGIONS[oi].meshpoints[omesh].final_flux_e[eng] = REGIONS[oi].meshpoints[omesh].final_flux_e[eng] + Wmew[mew1]*REGIONS[oi].meshpoints[omesh].total_flux[1][0][0][mew1][eng] / 2.0;
                        }
                    }
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
            /// Outputting Region Lower Bounds
            for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                Isos4 << REGIONS[oi].x_lower << "   " << REGIONS[oi].r_type << endl;
            }
            /// Outputting Region Absorbtion Rates

            for(unsigned int oi = 0; oi < REGIONS.size(); oi++){
                for(unsigned int omesh =0; omesh<REGIONS[oi].meshpoints.size(); omesh++){
                    for(unsigned int eng = 0; eng < bins; eng ++){
                        REGIONS[oi].meshpoints[omesh].final_flux = REGIONS[oi].meshpoints[omesh].final_flux + REGIONS[oi].meshpoints[omesh].final_flux_e[eng]*REGIONS[oi].ABSORB[eng];
                    }
                    Isos5 << REGIONS[oi].meshpoints[omesh].x_int << "   "<< REGIONS[oi].meshpoints[omesh].final_flux << endl;
                }
            }

            for(int mew1 = 0; mew1 < Ordinate; mew1++){
                for(int eng = 0; eng < bins; eng ++){
                    Isos6 << REGIONS[1].meshpoints[REGIONS[1].meshpoints.size()/2].flux[1][0][0][mew1][eng] * Wmew[mew1] << "    ";
                }
                Isos6 << endl;
            }
            Isos3 << endl;
            Isos4 << -leak_left/(SOURCE*3*REGION_SIZE) << "  " << leak_right/(SOURCE*3*REGION_SIZE) << endl << count;
            Isos4 << endl;
            Isos5 << endl;
            Isos6 << endl;

            for(int ci = 0; ci < REGIONS.size(); ci++){
                REGIONS[ci].meshpoints.clear();
                REGIONS[ci].mesh_points=0;
                REGIONS[ci].mesh_size=0;
                REGIONS[ci].max_mesh=0;
                REGIONS[ci].SKERNAL.clear();
                REGIONS[ci].TOTAL.clear();
                REGIONS[ci].ABSORB.clear();
                REGIONS[ci].SKATTER.clear();
                REGIONS[ci].vIsos.clear();
            }

            total_count = total_count +1;
    }
    return 0;
}

