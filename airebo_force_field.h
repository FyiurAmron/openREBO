// written by Szymon Winczewski
// based on LAMMPS implementation of AIREBO force field

#ifndef airebo_force_field_h
#define airebo_force_field_h

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>


const double AIREBO_PI = 4.0 * atan(1.0);

#define TOL 1.0e-9

// #define AIREBO_DEBUG


typedef struct
{
    double x, y, z;
    double r, r_sq;
} vec3d;


class AIREBOForceField
{
public:
    AIREBOForceField(std::string para_file_name, double para_cutlj,
                     bool para_ljflag, bool para_torflag, int para_max_number_of_REBO_neighbours);
    ~AIREBOForceField();

    double getCutoffRadius();
    double compute(int para_number_of_atoms, int *para_type,
                   int *para_neighbours_num, int **para_neighbours_list,
                   vec3d **para_neighbours_bonds);

    double getTotalEnergy();
    double getREBOEnergy();
    double getLJEnergy();
    double getTORSIONEnergy();


private:
// ogolne parametry obliczen
    int natoms;
    int *type;
    bool allocated;

// sumy funkcji wagowych (czlon REBO)
    double *nC, *nH;

// sasiedzi REBO
    int *REBO_neighbours_num;
    int **REBO_neighbours_list;
    vec3d **REBO_neighbours_bonds;

// wszyscy sasiedzi
    int *neighbours_num;
    int **neighbours_list;
    vec3d **neighbours_bonds;

// promien odciecia dla czlonu LJ, wyrazony w wielokrotnosciach parametru sigma potencjalu LJ
    double cutlj;
    double **cutljsq;
    double **lj1, **lj2, **lj3, **lj4;
    double cut3rebo, cutljrebo, cutljrebosq;
    double cutmax;

// flagi: czy uwzgledniac czlon lj oraz czlon torsyjny
    bool ljflag, torflag;
    int max_number_of_REBO_neighbours;

    double total_energy, energy_rebo, energy_lj, energy_torsion;

// parametry potencjalu AIREBO
// a) czlon REBO
    double rcmin_CC, rcmin_CH, rcmin_HH;
    double rcmax_CC, rcmax_CH, rcmax_HH;
    double rcmin[2][2], rcmax[2][2];
    double rcmaxsq[2][2];
    double rcmaxp_CC, rcmaxp_CH, rcmaxp_HH;
    double rcmaxp[2][2];

    double smin;
    double Nmin, Nmax;
    double NCmin, NCmax;
    double Q_CC, Q_CH, Q_HH;
    double Q[2][2];
    double alpha_CC, alpha_CH, alpha_HH;
    double alpha[2][2];
    double A_CC, A_CH, A_HH;
    double A[2][2];
    double BIJc_CC1, BIJc_CC2, BIJc_CC3;
    double BIJc_CH1, BIJc_CH2, BIJc_CH3;
    double BIJc_HH1, BIJc_HH2, BIJc_HH3;
    double BIJc[2][2][3];
    double Beta_CC1, Beta_CC2, Beta_CC3;
    double Beta_CH1, Beta_CH2, Beta_CH3;
    double Beta_HH1, Beta_HH2, Beta_HH3;
    double Beta[2][2][3];
    double rho_CC, rho_CH, rho_HH;
    double rho[2][2];

// b) czlon LJ
    double rcLJmin_CC, rcLJmin_CH, rcLJmin_HH;
    double rcLJmin[2][2];
    double rcLJmax_CC, rcLJmax_CH, rcLJmax_HH;
    double rcLJmax[2][2];
    double rcLJmaxsq[2][2];
    double bLJmin_CC, bLJmin_CH, bLJmin_HH;
    double bLJmin[2][2];
    double bLJmax_CC, bLJmax_CH, bLJmax_HH;
    double bLJmax[2][2];
    double epsilon_CC, epsilon_CH, epsilon_HH;
    double epsilon[2][2];
    double sigma_CC, sigma_CH, sigma_HH;
    double sigma[2][2];

// c) czlon torsyjny
    double thmin, thmax;
    double epsilonT_CCCC, epsilonT_CCCH, epsilonT_HCCH;
    double epsilonT[2][2];

// d) splajny
    double gCdom[5], gC1[4][6], gC2[4][6];
    double gHdom[4], gH[3][6];
    double pCCdom[2][2], pCC[4][4][16];
    double pCHdom[2][2], pCH[4][4][16];
    double piCCdom[3][2], piCC[4][4][9][64];
    double piCHdom[3][2], piCH[4][4][9][64];
    double piHHdom[3][2], piHH[4][4][9][64];
    double Tijdom[3][2], Tijc[4][4][9][64];

// splajny (polozenia wezlow)
    double PCCf[5][5], PCCdfdx[5][5], PCCdfdy[5][5];
    double PCHf[5][5], PCHdfdx[5][5], PCHdfdy[5][5];
    double piCCf[5][5][11], piCCdfdx[5][5][11], piCCdfdy[5][5][11], piCCdfdz[5][5][11];
    double piCHf[5][5][11], piCHdfdx[5][5][11], piCHdfdy[5][5][11], piCHdfdz[5][5][11];
    double piHHf[5][5][11], piHHdfdx[5][5][11], piHHdfdy[5][5][11], piHHdfdz[5][5][11];
    double Tf[5][5][10], Tdfdx[5][5][10], Tdfdy[5][5][10], Tdfdz[5][5][10];




    void readParameters(std::string file_name);
    void allocateMemory();
    void deallocateMemory();
    void initialize_constants();
    void initialize_splines();

// odczyt parametrow potencjalu
    void getLine(std::ifstream &file);
    double getLineAndConvertToDouble(std::ifstream &file);
    int getLineAndConvertToInt(std::ifstream &file);

// funkcje matematyczne
    double min(double val1, double val2);
    double max(double val1, double val2);
    double square(double arg);
    double cube(double arg);
    double pow4(double arg);
    double pow5(double arg);

// delta Kroneckera
    double kronecker(const int a, const int b) const;

// obliczanie funkcji odciecia
    double Sp(double Xij, double Xmin, double Xmax) const;
    double Sp2(double Xij, double Xmin, double Xmax) const;

// ewaluacja splajnow
    double gSpline(double costh, double Nij, int typei);
    double PijSpline(double NijC, double NijH, int typei, int typej);
    double piRCSpline(double Nij, double Nji, double Nijconj, int typei, int typej);
    double TijSpline(double Nij, double Nji, double Nijconj);

    double Sp5th(double x, double *coeffs);
    double Spbicubic(double x, double y, double *coeffs);
    double Sptricubic(double x, double y, double z, double *coeffs);

    double bondorder(int i, int j, double *rji_vec, double rji);
    double bondorderLJ(int i, int j, double *rji_vec, double rji, double rji0);

    void REBO_neighbours();
    void E_REBO();
    void E_LJ();
    void E_TORSION();
};


inline void AIREBOForceField::getLine(std::ifstream &file)
{
    std::string line;

    getline(file, line);
    if ( file.eof() == 1 )
    {
        std::cout << "error in AIREBOForceField::(): end of file reached!" << std::endl;
        exit(0);
    }
}


inline double AIREBOForceField ::getLineAndConvertToDouble(std::ifstream &file)
{
    std::string line;
    std::istringstream iss;
    double value;

    getline(file, line);
    if ( file.eof() == 1 )
    {
        std::cout << "error in AIREBOForceField::getLineAndConvertToDouble(): end of file reached!" << std::endl;
        exit(0);
    }

    iss.str(line);
    iss >> value;
    if ( value != value )
    {
        std::cout << "error in AIREBOForceField::getLineAndConvertToDouble(): conversion error!" << std::endl;
        exit(0);
    }

    return value;
}


inline int AIREBOForceField ::getLineAndConvertToInt(std::ifstream &file)
{
    std::string line;
    std::istringstream iss;
    int value;

    getline(file, line);
    if ( file.eof() == 1 )
    {
        std::cout << "error in AIREBOForceField::getLineAndConvertToInt(): end of file reached!" << std::endl;
        exit(0);
    }

    iss.str(line);
    iss >> value;
    if ( value != value )
    {
        std::cout << "error in AIREBOForceField::getLineAndConvertToInt(): conversion error!" << std::endl;
        exit(0);
    }

    return value;
}


inline double AIREBOForceField::min(double val1, double val2)
{
    if ( val1 < val2 )
        return val1;
    else
        return val2;
}


inline double AIREBOForceField::max(double val1, double val2)
{
    if ( val1 > val2 )
        return val1;
    else
        return val2;
}


inline double AIREBOForceField::square(double arg)
{
    return arg * arg;
}


inline double AIREBOForceField::cube(double arg)
{
    return arg * arg * arg;
}


inline double AIREBOForceField::pow4(double arg)
{
    return arg * arg * arg * arg;
}


inline double AIREBOForceField::pow5(double arg)
{
    return arg * arg * arg * arg * arg;
}


inline double AIREBOForceField::kronecker(const int a, const int b) const
{
    return (a == b) ? 1.0 : 0.0;
}


inline double AIREBOForceField::Sp(double Xij, double Xmin, double Xmax) const
{
// funkcja odciecia Sp
    double cutoff;
    double t = ( Xij - Xmin ) / ( Xmax - Xmin );
// !!!nieoptymalnie, konieczne odraczanie dzielenia
    if ( t <= 0.0 )
        cutoff = 1.0;
    else if ( t >= 1.0 )
        cutoff = 0.0;
    else
        cutoff = 0.5 * ( 1.0 + cos( t * AIREBO_PI ) );
    return cutoff;
}


inline double AIREBOForceField::Sp2(double Xij, double Xmin, double Xmax) const
{
// funkcja odciecia Sp2
    double cutoff;
    double t = ( Xij - Xmin ) / ( Xmax - Xmin );
// !!!nieoptymalnie, konieczne odraczanie dzielenia
    if ( t <= 0.0 )
        cutoff = 1.0;
    else if ( t >= 1.0 )
        cutoff = 0.0;
    else
        cutoff = ( 1.0 - ( t * t * ( 3.0 - 2.0 * t ) ) );
    return cutoff;
}

#endif
