// written by Szymon Winczewski
// based on LAMMPS implementation of AIREBO force field

#include "airebo_force_field.h"

using namespace std;


AIREBOForceField::AIREBOForceField(std::string para_file_name, double para_cutlj,
                                   bool para_ljflag, bool para_torflag, int para_max_number_of_REBO_neighbours)
{
    natoms = 0;
    type = NULL;
    allocated = 0;
    nC = NULL;
    nH = NULL;
    cutljsq = NULL;
    lj1 = NULL;
    lj2 = NULL;
    lj3 = NULL;
    lj4 = NULL;

    total_energy   = 0.0;
    energy_rebo    = 0.0;
    energy_lj      = 0.0;
    energy_torsion = 0.0;

    readParameters(para_file_name);
    if ( para_cutlj <= 0.0 )
    {
        cout << "error in AIREBOForceField::AIREBOForceField(): incorrect parameter (para_cutlj)!" << endl;
        exit(0);
    }
    cutlj = para_cutlj;
    ljflag = para_ljflag;
    torflag = para_torflag;

    if ( para_max_number_of_REBO_neighbours < 1 )
    {
        cout << "error in AIREBOForceField::AIREBOForceField(): incorrect parameter (para_max_number_of_REBO_neighbours)!" << endl;
        exit(0);
    }
    max_number_of_REBO_neighbours = para_max_number_of_REBO_neighbours;

    allocateMemory();
    initialize_constants();
    initialize_splines();
}


AIREBOForceField::~AIREBOForceField()
{
    deallocateMemory();
}


double AIREBOForceField::getCutoffRadius()
{
    return cutmax;
}


double AIREBOForceField::compute(int para_number_of_atoms, int *para_type,
                                 int *para_neighbours_num, int **para_neighbours_list,
                                 vec3d **para_neighbours_bonds)
{
    #ifdef AIREBO_DEBUG
    cout << "AIREBOForceField::compute() started!" << endl;
    cout << "   para_number_of_atoms    " << para_number_of_atoms << endl;
    cout << "   para_type               " << para_type << endl;
    cout << "   para_neighbours_num     " << para_neighbours_num << endl;
    cout << "   para_neighbours_list    " << para_neighbours_list << endl;
    cout << "   para_neighbours_bonds   " << para_neighbours_bonds << endl;
    cout << endl;
    #endif

    int i;

    if ( para_number_of_atoms <= 0 )
    {
        cout << "error in AIREBOForceField::compute(): incorrect parameter (para_number_of_atoms)" << endl;
        exit(0);
    }

// zaalokowano juz wczesniej pamiec
    if ( allocated == 1 )
    {
// potrzebna jest pamiec o innym rozmiarze
        if ( natoms != para_number_of_atoms )
        {
            delete [] nC;
            delete [] nH;
            for (i = 0; i < natoms; i++)
            {
                delete [] REBO_neighbours_list[i];
                delete [] REBO_neighbours_bonds[i];
            }
            delete [] REBO_neighbours_num;
            delete [] REBO_neighbours_list;
            delete [] REBO_neighbours_bonds;
            allocated = 0;

            natoms = para_number_of_atoms;
            nC = new double [natoms];
            nH = new double [natoms];
            REBO_neighbours_num = new int [natoms];
            REBO_neighbours_list = new int *[natoms];
            REBO_neighbours_bonds = new vec3d *[natoms];
            for (i = 0; i < natoms; i++)
            {
                REBO_neighbours_list[i] = new int [max_number_of_REBO_neighbours];
                REBO_neighbours_bonds[i] = new vec3d [max_number_of_REBO_neighbours];
            }
            allocated = 1;
        }
    }
    else
    {
// nie zaalokowano wczesniej pamieci
        natoms = para_number_of_atoms;
        nC = new double [natoms];
        nH = new double [natoms];
        REBO_neighbours_num = new int [natoms];
        REBO_neighbours_list = new int *[natoms];
        REBO_neighbours_bonds = new vec3d *[natoms];
        for (i = 0; i < natoms; i++)
        {
            REBO_neighbours_list[i] = new int [max_number_of_REBO_neighbours];
            REBO_neighbours_bonds[i] = new vec3d [max_number_of_REBO_neighbours];
        }
        allocated = 1;
    }

    if ( para_type == NULL )
    {
        cout << "error in AIREBOForceField::compute(): incorrect (NULL) parameter (para_type)" << endl;
        exit(0);
    }
    type = para_type;

    if ( para_neighbours_num == NULL )
    {
        cout << "error in AIREBOForceField::compute(): incorrect (NULL) parameter (para_neighbours_num)" << endl;
        exit(0);
    }
    neighbours_num = para_neighbours_num;

    if ( para_neighbours_list == NULL )
    {
        cout << "error in AIREBOForceField::compute(): incorrect (NULL) parameter (para_neighbours_list)" << endl;
        exit(0);
    }
    neighbours_list = para_neighbours_list;

    if ( para_neighbours_bonds == NULL )
    {
        cout << "error in AIREBOForceField::compute(): incorrect (NULL) parameter (para_neighbours_bonds)" << endl;
        exit(0);
    }
    neighbours_bonds = para_neighbours_bonds;

    total_energy   = 0.0;
    energy_rebo    = 0.0;
    energy_lj      = 0.0;
    energy_torsion = 0.0;

    #ifdef AIREBO_DEBUG
    cout << "starting REBO_neighbours()!" << endl;
    #endif

    REBO_neighbours();

    #ifdef AIREBO_DEBUG
    cout << "REBO_neighbours() finished!" << endl;
    #endif

    #ifdef AIREBO_DEBUG
    cout << "starting E_REBO()!" << endl;
    #endif

    E_REBO();

    #ifdef AIREBO_DEBUG
    cout << "E_REBO() finished!" << endl;
    cout << "energy_rebo = " << energy_rebo << endl;
    cout << endl;
    #endif

    if ( ljflag == 1 )
    {
        #ifdef AIREBO_DEBUG
        cout << "starting E_LJ()!" << endl;
        #endif

        E_LJ();

        #ifdef AIREBO_DEBUG
        cout << "E_LJ() finished!" << endl;
        cout << "energy_lj = " << energy_lj << endl;
        cout << endl;
        #endif
    }

    if ( torflag == 1 )
    {
        #ifdef AIREBO_DEBUG
        cout << "starting E_TORSION()!" << endl;
        #endif

        E_TORSION();

        #ifdef AIREBO_DEBUG
        cout << "E_TORSION() finished!" << endl;
        cout << "energy_torsion = " << energy_torsion << endl;
        cout << endl;
        #endif
    }

    total_energy = energy_rebo + energy_lj + energy_torsion;

    return total_energy;
}


double AIREBOForceField::getTotalEnergy()
{
    return total_energy;
}


double AIREBOForceField::getREBOEnergy()
{
    return energy_rebo;
}


double AIREBOForceField::getLJEnergy()
{
    return energy_lj;
}


double AIREBOForceField::getTORSIONEnergy()
{
    return energy_torsion;
}


void AIREBOForceField::readParameters(string file_name)
{
// wczytuje i parsuje plik z parametrami potencjalu AIREBO
    ifstream file;
    string line;
    int i, j, k, l;
    int number_of_domains;

// otwieramy plik do odczytu
    file.open(file_name.c_str());
    if ( file.is_open() == 0 )
    {
        cout << "error in AIREBOForceField::readParameters(): could not locate AIREBO potential file \"" << file_name << "\"!" << endl;
        exit(0);
    }

// pomijamy poczatkowe linie stanowiace komentarz
    while (1)
    {
        getline(file, line);
        if ( line.length() == 0 )
            break;
        else if ( line.at(0) != '#' )
            break;
    }

// a) czlon REBO
    rcmin_CC = getLineAndConvertToDouble(file);
    rcmin_CH = getLineAndConvertToDouble(file);
    rcmin_HH = getLineAndConvertToDouble(file);
    rcmax_CC = getLineAndConvertToDouble(file);
    rcmax_CH = getLineAndConvertToDouble(file);
    rcmax_HH = getLineAndConvertToDouble(file);
    rcmaxp_CC = getLineAndConvertToDouble(file);
    rcmaxp_CH = getLineAndConvertToDouble(file);
    rcmaxp_HH = getLineAndConvertToDouble(file);

    smin = getLineAndConvertToDouble(file);
    Nmin = getLineAndConvertToDouble(file);
    Nmax = getLineAndConvertToDouble(file);
    NCmin = getLineAndConvertToDouble(file);
    NCmax = getLineAndConvertToDouble(file);

    Q_CC = getLineAndConvertToDouble(file);
    Q_CH = getLineAndConvertToDouble(file);
    Q_HH = getLineAndConvertToDouble(file);
    alpha_CC = getLineAndConvertToDouble(file);
    alpha_CH = getLineAndConvertToDouble(file);
    alpha_HH = getLineAndConvertToDouble(file);
    A_CC = getLineAndConvertToDouble(file);
    A_CH = getLineAndConvertToDouble(file);
    A_HH = getLineAndConvertToDouble(file);

    BIJc_CC1 = getLineAndConvertToDouble(file);
    BIJc_CC2 = getLineAndConvertToDouble(file);
    BIJc_CC3 = getLineAndConvertToDouble(file);
    BIJc_CH1 = getLineAndConvertToDouble(file);
    BIJc_CH2 = getLineAndConvertToDouble(file);
    BIJc_CH3 = getLineAndConvertToDouble(file);
    BIJc_HH1 = getLineAndConvertToDouble(file);
    BIJc_HH2 = getLineAndConvertToDouble(file);
    BIJc_HH3 = getLineAndConvertToDouble(file);

    Beta_CC1 = getLineAndConvertToDouble(file);
    Beta_CC2 = getLineAndConvertToDouble(file);
    Beta_CC3 = getLineAndConvertToDouble(file);
    Beta_CH1 = getLineAndConvertToDouble(file);
    Beta_CH2 = getLineAndConvertToDouble(file);
    Beta_CH3 = getLineAndConvertToDouble(file);
    Beta_HH1 = getLineAndConvertToDouble(file);
    Beta_HH2 = getLineAndConvertToDouble(file);
    Beta_HH3 = getLineAndConvertToDouble(file);

    rho_CC = getLineAndConvertToDouble(file);
    rho_CH = getLineAndConvertToDouble(file);
    rho_HH = getLineAndConvertToDouble(file);

// b) czlon LJ
    rcLJmin_CC = getLineAndConvertToDouble(file);
    rcLJmin_CH = getLineAndConvertToDouble(file);
    rcLJmin_HH = getLineAndConvertToDouble(file);
    rcLJmax_CC = getLineAndConvertToDouble(file);
    rcLJmax_CH = getLineAndConvertToDouble(file);
    rcLJmax_HH = getLineAndConvertToDouble(file);

    bLJmin_CC = getLineAndConvertToDouble(file);
    bLJmin_CH = getLineAndConvertToDouble(file);
    bLJmin_HH = getLineAndConvertToDouble(file);
    bLJmax_CC = getLineAndConvertToDouble(file);
    bLJmax_CH = getLineAndConvertToDouble(file);
    bLJmax_HH = getLineAndConvertToDouble(file);

    epsilon_CC = getLineAndConvertToDouble(file);
    epsilon_CH = getLineAndConvertToDouble(file);
    epsilon_HH = getLineAndConvertToDouble(file);
    sigma_CC = getLineAndConvertToDouble(file);
    sigma_CH = getLineAndConvertToDouble(file);
    sigma_HH = getLineAndConvertToDouble(file);

// c) czlon REBO
    epsilonT_CCCC = getLineAndConvertToDouble(file);
    epsilonT_CCCH = getLineAndConvertToDouble(file);
    epsilonT_HCCH = getLineAndConvertToDouble(file);

// d) splajny
// splajny gC1 i gC2
    getLine(file); getLine(file); getLine(file);
// liczba wezlow
    number_of_domains = getLineAndConvertToInt(file);
// gCdom
    for (i = 0; i < number_of_domains; i++)
        gCdom[i] = getLineAndConvertToDouble(file);
    getLine(file);
// gC1
    for (i = 0; i < number_of_domains - 1; i++)
        for (j = 0; j < 6; j++)
            gC1[i][j] = getLineAndConvertToDouble(file);
    getLine(file);
// gC2
    for (i = 0; i < number_of_domains - 1; i++)
        for (j = 0; j < 6; j++)
            gC2[i][j] = getLineAndConvertToDouble(file);

// splajn gH
    getLine(file); getLine(file); getLine(file);
// liczba wezlow
    number_of_domains = getLineAndConvertToInt(file);
// gHdom
    for (i = 0; i < number_of_domains; i++)
        gHdom[i] = getLineAndConvertToDouble(file);
    getLine(file);
// gH
    for (i = 0; i < number_of_domains - 1; i++)
        for (j = 0; j < 6; j++)
            gH[i][j] = getLineAndConvertToDouble(file);

// splajn pCC
    getLine(file); getLine(file); getLine(file);
// liczba wezlow
    number_of_domains = getLineAndConvertToInt(file);
// pCCdom
    for (i = 0; i < number_of_domains / 2; i++)
        for (j = 0; j < number_of_domains / 2; j++)
            pCCdom[i][j] = getLineAndConvertToDouble(file);
    getLine(file);
// pCC
    for (i = 0; i < (int) pCCdom[0][1]; i++)
        for (j = 0; j < (int) pCCdom[1][1]; j++)
            for (k = 0; k < 16; k++)
                pCC[i][j][k] = getLineAndConvertToDouble(file);

// splajn pCH
    getLine(file); getLine(file); getLine(file);
// liczba wezlow
    number_of_domains = getLineAndConvertToInt(file);
// pCHdom
    for (i = 0; i < number_of_domains / 2; i++)
        for (j = 0; j < number_of_domains / 2; j++)
            pCHdom[i][j] = getLineAndConvertToDouble(file);
    getLine(file);
// pCH
    for (i = 0; i < (int) pCHdom[0][1]; i++)
        for (j = 0; j < (int) pCHdom[1][1]; j++)
            for (k = 0; k < 16; k++)
                pCH[i][j][k] = getLineAndConvertToDouble(file);

// splajn piCC
    getLine(file); getLine(file); getLine(file);
// liczba wezlow
    number_of_domains = getLineAndConvertToInt(file);
// piCCdom
    for (i = 0; i < number_of_domains / 2; i++)
        for (j = 0; j < number_of_domains / 3; j++)
            piCCdom[i][j] = getLineAndConvertToDouble(file);
    getLine(file);
// piCC
    for (i = 0; i < (int) piCCdom[0][1]; i++)
        for (j = 0; j < (int) piCCdom[1][1]; j++)
            for (k = 0; k < (int) piCCdom[2][1]; k++)
                for (l = 0; l < 64; l++)
                    piCC[i][j][k][l] = getLineAndConvertToDouble(file);

// splajn piCH
    getLine(file); getLine(file); getLine(file);
    number_of_domains = getLineAndConvertToInt(file);
// piCHdom
    for (i = 0; i < number_of_domains / 2; i++)
        for (j = 0; j < number_of_domains / 3; j++)
            piCHdom[i][j] = getLineAndConvertToDouble(file);
    getLine(file);
// piCH
    for (i = 0; i < (int) piCHdom[0][1]; i++)
        for (j = 0; j < (int) piCHdom[1][1]; j++)
            for (k = 0; k < (int) piCHdom[2][1]; k++)
                for (l = 0; l < 64; l++)
                    piCH[i][j][k][l] = getLineAndConvertToDouble(file);

// splajn piHH
    getLine(file); getLine(file); getLine(file);
    number_of_domains = getLineAndConvertToInt(file);
// piHHdom
    for (i = 0; i < number_of_domains / 2; i++)
        for (j = 0; j < number_of_domains / 3; j++)
            piHHdom[i][j] = getLineAndConvertToDouble(file);
    getLine(file);
// piHH
    for (i = 0; i < (int) piHHdom[0][1]; i++)
        for (j = 0; j < (int) piHHdom[1][1]; j++)
            for (k = 0; k < (int) piHHdom[2][1]; k++)
                for (l = 0; l < 64; l++)
                    piHH[i][j][k][l] = getLineAndConvertToDouble(file);

// splajn Tij
    getLine(file); getLine(file); getLine(file);
    number_of_domains = getLineAndConvertToInt(file);
// Tijdom
    for (i = 0; i < number_of_domains / 2; i++)
        for (j = 0; j < number_of_domains / 3; j++)
            Tijdom[i][j] = getLineAndConvertToDouble(file);
    getLine(file);
// Tijc
    for (i = 0; i < (int) Tijdom[0][1]; i++)
        for (j = 0; j < (int) Tijdom[1][1]; j++)
            for (k = 0; k < (int) Tijdom[2][1]; k++)
                for (l = 0; l < 64; l++)
                    Tijc[i][j][k][l] = getLineAndConvertToDouble(file);
// zamykamy plik
    file.close();

// przepisujemy wczytane parametry do macierzy
// a) czlon REBO
    rcmin[0][0] = rcmin_CC;
    rcmin[0][1] = rcmin_CH;
    rcmin[1][0] = rcmin_CH;
    rcmin[1][1] = rcmin_HH;

    rcmax[0][0] = rcmax_CC;
    rcmax[0][1] = rcmax_CH;
    rcmax[1][0] = rcmax_CH;
    rcmax[1][1] = rcmax_HH;

    rcmaxsq[0][0] = rcmax[0][0] * rcmax[0][0];
    rcmaxsq[1][0] = rcmax[1][0] * rcmax[1][0];
    rcmaxsq[0][1] = rcmax[0][1] * rcmax[0][1];
    rcmaxsq[1][1] = rcmax[1][1] * rcmax[1][1];

    rcmaxp[0][0] = rcmaxp_CC;
    rcmaxp[0][1] = rcmaxp_CH;
    rcmaxp[1][0] = rcmaxp_CH;
    rcmaxp[1][1] = rcmaxp_HH;

    Q[0][0] = Q_CC;
    Q[0][1] = Q_CH;
    Q[1][0] = Q_CH;
    Q[1][1] = Q_HH;

    alpha[0][0] = alpha_CC;
    alpha[0][1] = alpha_CH;
    alpha[1][0] = alpha_CH;
    alpha[1][1] = alpha_HH;

    A[0][0] = A_CC;
    A[0][1] = A_CH;
    A[1][0] = A_CH;
    A[1][1] = A_HH;

    BIJc[0][0][0] = BIJc_CC1;
    BIJc[0][0][1] = BIJc_CC2;
    BIJc[0][0][2] = BIJc_CC3;
    BIJc[0][1][0] = BIJc_CH1;
    BIJc[0][1][1] = BIJc_CH2;
    BIJc[0][1][2] = BIJc_CH3;
    BIJc[1][0][0] = BIJc_CH1;
    BIJc[1][0][1] = BIJc_CH2;
    BIJc[1][0][2] = BIJc_CH3;
    BIJc[1][1][0] = BIJc_HH1;
    BIJc[1][1][1] = BIJc_HH2;
    BIJc[1][1][2] = BIJc_HH3;

    Beta[0][0][0] = Beta_CC1;
    Beta[0][0][1] = Beta_CC2;
    Beta[0][0][2] = Beta_CC3;
    Beta[0][1][0] = Beta_CH1;
    Beta[0][1][1] = Beta_CH2;
    Beta[0][1][2] = Beta_CH3;
    Beta[1][0][0] = Beta_CH1;
    Beta[1][0][1] = Beta_CH2;
    Beta[1][0][2] = Beta_CH3;
    Beta[1][1][0] = Beta_HH1;
    Beta[1][1][1] = Beta_HH2;
    Beta[1][1][2] = Beta_HH3;

    rho[0][0] = rho_CC;
    rho[0][1] = rho_CH;
    rho[1][0] = rho_CH;
    rho[1][1] = rho_HH;

// b) czlon LJ
    rcLJmin[0][0] = rcLJmin_CC;
    rcLJmin[0][1] = rcLJmin_CH;
    rcLJmin[1][0] = rcLJmin_CH;
    rcLJmin[1][1] = rcLJmin_HH;

    rcLJmax[0][0] = rcLJmax_CC;
    rcLJmax[0][1] = rcLJmax_CH;
    rcLJmax[1][0] = rcLJmax_CH;
    rcLJmax[1][1] = rcLJmax_HH;

    rcLJmaxsq[0][0] = rcLJmax[0][0] * rcLJmax[0][0];
    rcLJmaxsq[1][0] = rcLJmax[1][0] * rcLJmax[1][0];
    rcLJmaxsq[0][1] = rcLJmax[0][1] * rcLJmax[0][1];
    rcLJmaxsq[1][1] = rcLJmax[1][1] * rcLJmax[1][1];

    bLJmin[0][0] = bLJmin_CC;
    bLJmin[0][1] = bLJmin_CH;
    bLJmin[1][0] = bLJmin_CH;
    bLJmin[1][1] = bLJmin_HH;

    bLJmax[0][0] = bLJmax_CC;
    bLJmax[0][1] = bLJmax_CH;
    bLJmax[1][0] = bLJmax_CH;
    bLJmax[1][1] = bLJmax_HH;

    epsilon[0][0] = epsilon_CC;
    epsilon[0][1] = epsilon_CH;
    epsilon[1][0] = epsilon_CH;
    epsilon[1][1] = epsilon_HH;

    sigma[0][0] = sigma_CC;
    sigma[0][1] = sigma_CH;
    sigma[1][0] = sigma_CH;
    sigma[1][1] = sigma_HH;

// d) czlon torsyjny
    thmin = -1.0;
    thmax = -0.995;
    epsilonT[0][0] = epsilonT_CCCC;
    epsilonT[0][1] = epsilonT_CCCH;
    epsilonT[1][0] = epsilonT_CCCH;
    epsilonT[1][1] = epsilonT_HCCH;
}


void AIREBOForceField::allocateMemory()
{
    cutljsq = new double *[2];
    cutljsq[0] = new double [2];
    cutljsq[1] = new double [2];

    lj1 = new double *[2];
    lj1[0] = new double [2];
    lj1[1] = new double [2];

    lj2 = new double *[2];
    lj2[0] = new double [2];
    lj2[1] = new double [2];

    lj3 = new double *[2];
    lj3[0] = new double [2];
    lj3[1] = new double [2];

    lj4 = new double *[2];
    lj4[0] = new double [2];
    lj4[1] = new double [2];
}


void AIREBOForceField::deallocateMemory()
{
    delete [] cutljsq[0];
    delete [] cutljsq[1];
    delete [] cutljsq;

    delete [] lj1[0];
    delete [] lj1[1];
    delete [] lj1;

    delete [] lj2[0];
    delete [] lj2[1];
    delete [] lj2;

    delete [] lj3[0];
    delete [] lj3[1];
    delete [] lj3;

    delete [] lj4[0];
    delete [] lj4[1];
    delete [] lj4;

    if ( allocated == 1 )
    {
        int i;
        delete [] nC;
        delete [] nH;
        for (i = 0; i < natoms; i++)
        {
            delete [] REBO_neighbours_list[i];
            delete [] REBO_neighbours_bonds[i];
        }
        delete [] REBO_neighbours_num;
        delete [] REBO_neighbours_list;
        delete [] REBO_neighbours_bonds;
    }
}


void AIREBOForceField::initialize_constants()
{
    int i, j;
    double tmp_cutmax;

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
            lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
            lj3[i][j] =  4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
            lj4[i][j] =  4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
            cutljsq[i][j] = cutlj * sigma[i][j] * cutlj * sigma[i][j];
        }
    }

// promien odciecia dla czlonu REBO zadany przez promien odciecia dla C,
// gdyz jest on najwiekszy
    cut3rebo = 3.0 * rcmax[0][0];
    cutljrebo = rcLJmax[0][0] + rcmax[0][0];
    cutljrebosq = cutljrebo * cutljrebo;

    cutmax = cut3rebo;

    if ( ljflag == 1 )
    {
        tmp_cutmax = rcLJmax[0][0] + 2.0 * rcmax[0][0];
        if ( tmp_cutmax > cutmax )
            cutmax = tmp_cutmax;
        tmp_cutmax = cutlj * sigma[0][0];
        if ( tmp_cutmax > cutmax )
            cutmax = tmp_cutmax;
    }
}


void AIREBOForceField::initialize_splines()
{
    int i, j, k;

    for (i = 0; i < 5; i++)
    {
        for (j = 0; j < 5; j++)
        {
            PCCf[i][j]    = 0.0;
            PCCdfdx[i][j] = 0.0;
            PCCdfdy[i][j] = 0.0;
            PCHf[i][j]    = 0.0;
            PCHdfdx[i][j] = 0.0;
            PCHdfdy[i][j] = 0.0;
        }
    }

    PCCf[0][2] = -0.00050;
    PCCf[0][3] =  0.0161253646;
    PCCf[1][1] = -0.010960;
    PCCf[1][2] =  0.00632624824;
    PCCf[2][0] = -0.0276030;
    PCCf[2][1] =  0.00317953083;

    PCHf[0][1] =  0.209336733;
    PCHf[0][2] = -0.0644496154;
    PCHf[0][3] = -0.303927546;
    PCHf[1][0] =  0.010;
    PCHf[1][1] = -0.125123401;
    PCHf[1][2] = -0.298905246;
    PCHf[2][0] = -0.122042146;
    PCHf[2][1] = -0.300529172;
    PCHf[3][0] = -0.307584705;

    for (i = 0; i < 5; i++)
    {
        for (j = 0; j < 5; j++)
        {
            for (k = 0; k < 10; k++)
            {
                piCCf[i][j][k]    = 0.0;
                piCCdfdx[i][j][k] = 0.0;
                piCCdfdy[i][j][k] = 0.0;
                piCCdfdz[i][j][k] = 0.0;
                piCHf[i][j][k]    = 0.0;
                piCHdfdx[i][j][k] = 0.0;
                piCHdfdy[i][j][k] = 0.0;
                piCHdfdz[i][j][k] = 0.0;
                piHHf[i][j][k]    = 0.0;
                piHHdfdx[i][j][k] = 0.0;
                piHHdfdy[i][j][k] = 0.0;
                piHHdfdz[i][j][k] = 0.0;
                Tf[i][j][k]       = 0.0;
                Tdfdx[i][j][k]    = 0.0;
                Tdfdy[i][j][k]    = 0.0;
                Tdfdz[i][j][k]    = 0.0;
            }
        }
    }

    for (i = 3; i < 10; i++)
        piCCf[0][0][i] = 0.0049586079;
    piCCf[1][0][1] = 0.021693495;
    piCCf[0][1][1] = 0.021693495;
    for (i = 2; i < 10; i++)
        piCCf[1][0][i] = 0.0049586079;
    for (i = 2; i < 10; i++)
        piCCf[0][1][i] = 0.0049586079;
    piCCf[1][1][1] = 0.05250;
    piCCf[1][1][2] = -0.002088750;
    for (i = 3; i < 10; i++)
        piCCf[1][1][i] = -0.00804280;
    piCCf[2][0][1] = 0.024698831850;
    piCCf[0][2][1] = 0.024698831850;
    piCCf[2][0][2] = -0.00597133450;
    piCCf[0][2][2] = -0.00597133450;
    for (i = 3; i < 10; i++)
        piCCf[2][0][i] = 0.0049586079;
    for (i = 3; i < 10; i++)
        piCCf[0][2][i] = 0.0049586079;
    piCCf[2][1][1] = 0.00482478490;
    piCCf[1][2][1] = 0.00482478490;
    piCCf[2][1][2] = 0.0150;
    piCCf[1][2][2] = 0.0150;
    piCCf[2][1][3] = -0.010;
    piCCf[1][2][3] = -0.010;
    piCCf[2][1][4] = -0.01168893870;
    piCCf[1][2][4] = -0.01168893870;
    piCCf[2][1][5] = -0.013377877400;
    piCCf[1][2][5] = -0.013377877400;
    piCCf[2][1][6] = -0.015066816000;
    piCCf[1][2][6] = -0.015066816000;
    for (i = 7; i < 10; i++)
        piCCf[2][1][i] = -0.015066816000;
    for (i = 7; i < 10; i++)
        piCCf[1][2][i] = -0.015066816000;
    piCCf[2][2][1] = 0.0472247850;
    piCCf[2][2][2] = 0.0110;
    piCCf[2][2][3] = 0.0198529350;
    piCCf[2][2][4] = 0.01654411250;
    piCCf[2][2][5] = 0.013235290;
    piCCf[2][2][6] = 0.00992646749999 ;
    piCCf[2][2][7] = 0.006617644999;
    piCCf[2][2][8] = 0.00330882250;
    piCCf[3][0][1] = -0.05989946750;
    piCCf[0][3][1] = -0.05989946750;
    piCCf[3][0][2] = -0.05989946750;
    piCCf[0][3][2] = -0.05989946750;
    for (i = 3; i < 10; i++)
        piCCf[3][0][i] = 0.0049586079;
    for (i = 3; i < 10; i++)
        piCCf[0][3][i] = 0.0049586079;
    piCCf[3][1][2] = -0.0624183760;
    piCCf[1][3][2] = -0.0624183760;
    for (i = 3; i < 10; i++)
        piCCf[3][1][i] = -0.0624183760;
    for (i = 3; i < 10; i++)
        piCCf[1][3][i] = -0.0624183760;
    piCCf[3][2][1] = -0.02235469150;
    piCCf[2][3][1] = -0.02235469150;
    for (i = 2; i < 10; i++)
        piCCf[3][2][i] = -0.02235469150;
    for (i = 2; i < 10; i++)
        piCCf[2][3][i] = -0.02235469150;

    piCCdfdx[2][1][1] = -0.026250;
    piCCdfdx[2][1][5] = -0.0271880;
    piCCdfdx[2][1][6] = -0.0271880;
    for (i = 7; i < 10; i++)
        piCCdfdx[2][1][i] = -0.0271880;
    piCCdfdx[1][3][2] = 0.0187723882;
    for (i = 2; i < 10; i++)
        piCCdfdx[2][3][i] = 0.031209;

    piCCdfdy[1][2][1] = -0.026250;
    piCCdfdy[1][2][5] = -0.0271880;
    piCCdfdy[1][2][6] = -0.0271880;
    for (i = 7; i < 10; i++)
        piCCdfdy[1][2][i] = -0.0271880;
    piCCdfdy[3][1][2] = 0.0187723882;
    for (i = 2; i < 10; i++)
        piCCdfdy[3][2][i] = 0.031209;

    piCCdfdz[1][1][2] = -0.0302715;
    piCCdfdz[2][1][4] = -0.0100220;
    piCCdfdz[1][2][4] = -0.0100220;
    piCCdfdz[2][1][5] = -0.0100220;
    piCCdfdz[1][2][5] = -0.0100220;
    for (i = 4; i < 9; i++)
        piCCdfdz[2][2][i] = -0.0033090;

// !!!komentarz: make top end of piCC flat instead of zero
//               also enforces some symmetry
    i = 4;
    for (j = 0; j < 4; j++)
        for (k = 1; k < 11; k++)
            piCCf[i][j][k] = piCCf[i-1][j][k];
    for (i = 0; i < 4; i++)
        for (j = i+1; j < 5; j++)
            for (k = 1; k < 11; k++)
                piCCf[i][j][k] = piCCf[j][i][k];
    for (k = 1; k < 11; k++)
        piCCf[4][4][k] = piCCf[3][4][k];
    k = 10;
    for (i = 0; i < 5; i++)
        for (j = 0; j < 5; j++)
            piCCf[i][j][k] = piCCf[i][j][k-1];

    piCHf[1][1][1] = -0.050;
    piCHf[1][1][2] = -0.050;
    piCHf[1][1][3] = -0.30;
    for (i = 4; i < 10; i++)
        piCHf[1][1][i] = -0.050;
    for (i = 5; i < 10; i++)
        piCHf[2][0][i] = -0.004523893758064;
    for (i = 5; i < 10; i++)
        piCHf[0][2][i] = -0.004523893758064;
    piCHf[2][1][2] = -0.250;
    piCHf[1][2][2] = -0.250;
    piCHf[2][1][3] = -0.250;
    piCHf[1][2][3] = -0.250;
    piCHf[3][1][1] = -0.10;
    piCHf[1][3][1] = -0.10;
    piCHf[3][1][2] = -0.125;
    piCHf[1][3][2] = -0.125;
    piCHf[3][1][3] = -0.125;
    piCHf[1][3][3] = -0.125;
    for (i = 4; i < 10; i++)
        piCHf[3][1][i] = -0.10;
    for (i = 4; i < 10; i++)
        piCHf[1][3][i] = -0.10;

// !!!komentarz: make top end of piCH flat instead of zero
//               also enforces some symmetry
    i = 4;
    for (j = 0; j < 4; j++)
        for (k = 1; k < 11; k++)
            piCHf[i][j][k] = piCHf[i-1][j][k];
    for (i = 0; i < 4; i++)
        for (j = i+1; j < 5; j++)
            for (k = 1; k < 11; k++)
                piCHf[i][j][k] = piCHf[j][i][k];
    for (k = 1; k < 11; k++)
        piCHf[4][4][k] = piCHf[3][4][k];
    k = 10;
    for (i = 0; i < 5; i++)
        for (j = 0; j < 5; j++)
            piCHf[i][j][k] = piCHf[i][j][k-1];

    piHHf[1][1][1] = 0.124915958;

    Tf[2][2][1] = -0.035140;
    for (i = 2; i < 10; i++)
        Tf[2][2][i] = -0.0040480;
}
