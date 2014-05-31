// written by Szymon Winczewski
// based on LAMMPS implementation of AIREBO force field

#include "airebo_force_field.h"

using namespace std;


void AIREBOForceField::E_LJ()
{
    double rljmin, rljmax;
    double sigcut, sigmin, sigwid;
    int i, j, k, m;
    int jj, kk, mm;
    int itype, jtype, ktype, mtype;

    int neighbours_num_i;
    int *neighbours_list_i;
    vec3d *neighbours_bonds_i;

    int REBO_neighbours_num_i;
    int *REBO_neighbours_list_i;
    vec3d *REBO_neighbours_bonds_i;

    int REBO_neighbours_num_k;
    int *REBO_neighbours_list_k;
    vec3d *REBO_neighbours_bonds_k;

    double best;
    int testpath, done;

    double rji, rji_sq, rji_vec[3];
    double rki, rki_sq, rki_vec[3];
    double rkj, rkj_sq, rkj_vec[3];
    double rmk, rmk_sq, rmk_vec[3];
    double rmj, rmj_sq, rmj_vec[3];

    double wki, wkj, wmk, wmj;

    double cij;
    double slw;
    double drji;
    double swidth, tee, tee2;
    double r2inv, r6inv;
    double vdw, VLJ, Str, VA, Stb;
    double scale, delscale[3];

    rljmin = 0.0;
    rljmax = 0.0;
    sigcut = 0.0;
    sigmin = 0.0;
    sigwid = 0.0;

    for (i = 0; i < natoms; i++)
    {
        itype = type[i];

        neighbours_num_i = neighbours_num[i];
        neighbours_list_i = neighbours_list[i];
        neighbours_bonds_i = neighbours_bonds[i];

        for (jj = 0; jj < neighbours_num_i; jj++)
        {
            j = neighbours_list_i[jj];

            if ( i < j )
                continue;

            jtype = type[j];

            rji_sq = neighbours_bonds_i[jj].r_sq;
            if ( rji_sq >= cutljsq[itype][jtype] )
                continue;

            rji = neighbours_bonds_i[jj].r;
            rji_vec[0] = neighbours_bonds_i[jj].x;
            rji_vec[1] = neighbours_bonds_i[jj].y;
            rji_vec[2] = neighbours_bonds_i[jj].z;

            if ( rji >= cut3rebo )
            {
                best = 0.0;
                testpath = 0;
            }
            else if ( rji >= rcMax[itype][jtype] )
            {
                best = 0.0;
                testpath = 1;
            }
            else
            {
                best = Sp(rji, rcMin[itype][jtype], rcMax[itype][jtype]);
                if ( best < 1.0 )
                    testpath = 1;
                else
                    testpath = 0;
            }

            done = 0;
            if ( testpath )
            {
                REBO_neighbours_num_i = REBO_neighbours_num[i];
                REBO_neighbours_list_i = REBO_neighbours_list[i];
                REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];

                for (kk = 0; ( ( kk < REBO_neighbours_num_i ) && ( done==0 ) ); kk++)
                {
                    k = REBO_neighbours_list_i[kk];

                    if ( k == j )
                        continue;

                    ktype = type[k];

                    rki_vec[0] = REBO_neighbours_bonds_i[kk].x;
                    rki_vec[1] = REBO_neighbours_bonds_i[kk].y;
                    rki_vec[2] = REBO_neighbours_bonds_i[kk].z;
                    rki = REBO_neighbours_bonds_i[kk].r;
                    rki_sq = REBO_neighbours_bonds_i[kk].r_sq;

                    if ( rki_sq < rcMaxSq[itype][ktype] )
                        wki = Sp(rki, rcMin[itype][ktype], rcMax[itype][ktype]);
                    else
                        wki = 0.0;

                    if ( wki > best )
                    {
                        rkj_vec[0] = rki_vec[0] - rji_vec[0];
                        rkj_vec[1] = rki_vec[1] - rji_vec[1];
                        rkj_vec[2] = rki_vec[2] - rji_vec[2];
                        rkj_sq = rkj_vec[0] * rkj_vec[0] + rkj_vec[1] * rkj_vec[1] + rkj_vec[2] * rkj_vec[2];

                        if ( rkj_sq < rcMaxSq[ktype][jtype] )
                        {
                            rkj = sqrt(rkj_sq);
                            wkj = Sp(rkj, rcMin[ktype][jtype], rcMax[ktype][jtype]);
                            if ( wki * wkj > best )
                            {
                                best = wki * wkj;
                                if ( best == 1.0 )
                                {
                                    done = 1;
                                    break;
                                }
                            }
                        }

                        REBO_neighbours_num_k = REBO_neighbours_num[k];
                        REBO_neighbours_list_k = REBO_neighbours_list[k];
                        REBO_neighbours_bonds_k = REBO_neighbours_bonds[k];
                        for (mm = 0; ( ( mm < REBO_neighbours_num_k ) && ( done == 0 ) ); mm++)
                        {
                            m = REBO_neighbours_list_k[mm];

                            if ( ( m == i ) || ( m == j ) )
                                continue;

                            mtype = type[m];

                            rmk_vec[0] = REBO_neighbours_bonds_k[mm].x;
                            rmk_vec[1] = REBO_neighbours_bonds_k[mm].y;
                            rmk_vec[2] = REBO_neighbours_bonds_k[mm].z;
                            rmk = REBO_neighbours_bonds_k[mm].r;
                            rmk_sq = REBO_neighbours_bonds_k[mm].r_sq;

                            if ( rmk_sq < rcMaxSq[ktype][mtype] )
                                wmk = Sp(rmk, rcMin[ktype][mtype], rcMax[ktype][mtype]);
                            else
                                wmk = 0.0;

                            if ( wki * wmk > best )
                            {
                                rmj_vec[0] = rmk_vec[0] + rkj_vec[0];
                                rmj_vec[1] = rmk_vec[1] + rkj_vec[1];
                                rmj_vec[2] = rmk_vec[2] + rkj_vec[2];
                                rmj_sq = rmj_vec[0] * rmj_vec[0] + rmj_vec[1] * rmj_vec[1] + rmj_vec[2] * rmj_vec[2];

                                if ( rmj_sq < rcMaxSq[mtype][jtype] )
                                {
                                    rmj = sqrt(rmj_sq);
                                    wmj = Sp(rmj, rcMin[mtype][jtype], rcMax[mtype][jtype]);
                                    if ( wki * wmk * wmj > best )
                                    {
                                        best = wki * wmk * wmj;
                                        if ( best == 1.0 )
                                        {
                                            done = 1;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            cij = 1.0 - best;
            if ( cij == 0.0 )
                continue;

            sigwid = 0.84;
            sigcut = 3.0;
            sigmin = sigcut - sigwid;

            rljmin = sigma[itype][jtype];
            rljmax = sigcut * rljmin;
            rljmin = sigmin * rljmin;

            if ( rji > rljmax )
                slw = 0.0;
            else if ( rji > rljmin )
            {
                drji = rji - rljmin;
                swidth = rljmax - rljmin;
                tee = drji / swidth;
                tee2 = tee * tee;
                slw = 1.0 - tee2 * ( 3.0 - 2.0 * tee );
            }
            else
                slw = 1.0;

            r2inv = 1.0 / rji_sq;
            r6inv = r2inv * r2inv * r2inv;

            vdw = r6inv * ( lj3[itype][jtype] * r6inv - lj4[itype][jtype] );
            VLJ = vdw * slw;

            Str = Sp2(rji, rcLJmin[itype][jtype], rcLJmax[itype][jtype]);
            VA = Str * cij * VLJ;
            if ( Str > 0.0 )
            {
                scale = rcMin[itype][jtype] / rji;
                delscale[0] = scale * rji_vec[0];
                delscale[1] = scale * rji_vec[1];
                delscale[2] = scale * rji_vec[2];
                Stb = bondorderLJ(i, j, delscale, rcMin[itype][jtype], rji);
            }
            else
                Stb = 0.0;

            energy_lj += VA * Stb + ( 1.0 - Str ) * cij * VLJ;
        }
    }
}


double AIREBOForceField::bondorderLJ(int i, int j, double rji_vec[3], double rji, double rji0)
{
    int atomi, itype;
    int atomj, jtype;
    int k, atomk, ktype;
    int l, atoml, ltype;
    double wij;
    double NijC, NijH, NjiC, NjiH;
    double bij;
    double NconjtmpI, NconjtmpJ;
    double Etmp;
    double Stb;

    int REBO_neighbours_num_i;
    int *REBO_neighbours_list_i;
    vec3d *REBO_neighbours_bonds_i;
    int REBO_neighbours_num_j;
    int *REBO_neighbours_list_j;
    vec3d *REBO_neighbours_bonds_j;

    double rji_sq;
    double rki_vec[3], rki, rki_sq;
    double rlj_vec[3], rlj, rlj_sq;
    double rkj_vec[3], rkj_sq;
    double rli_vec[3], rli_sq;

    double wki, wlj;
    double lamdajik, lamdaijl;
    double Nki, Nlj;
    double cosjik, cosijl;
    double g;

    double PijS, PjiS;
    double pij, pji;

    double Nijconj;
    double piRC, Tij;

    double cos321, cos234;
    double costmp;
    double tspjik, tspijl;

    double crosskij[3], crosskijmag;
    double crossijl[3], crossijlmag;
    double omkijl;

    atomi = i;
    atomj = j;
    itype = type[atomi];
    jtype = type[atomj];

    wij = Sp(rji0, rcMin[itype][jtype], rcMax[itype][jtype]);
    NijC = nC[atomi] - ( wij * kronecker(jtype, 0) );
    NijH = nH[atomi] - ( wij * kronecker(jtype, 1) );
    NjiC = nC[atomj] - ( wij * kronecker(itype, 0) );
    NjiH = nH[atomj] - ( wij * kronecker(itype, 1) );

    bij = 0.0;
    NconjtmpI = 0.0;
    NconjtmpJ = 0.0;
    Etmp = 0.0;
    Stb = 0.0;

    REBO_neighbours_num_i = REBO_neighbours_num[i];
    REBO_neighbours_list_i = REBO_neighbours_list[i];
    REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];
    for (k = 0; k < REBO_neighbours_num_i; k++)
    {
        atomk = REBO_neighbours_list_i[k];
        if ( atomk != atomj )
        {
            ktype = type[atomk];

            rki_vec[0] = REBO_neighbours_bonds_i[k].x;
            rki_vec[1] = REBO_neighbours_bonds_i[k].y;
            rki_vec[2] = REBO_neighbours_bonds_i[k].z;
            rki = REBO_neighbours_bonds_i[k].r;

            lamdajik = 4.0 * kronecker(itype, 1) *
                       ( ( rho[ktype][1] - rki ) - ( rho[jtype][1] - rji ) );
            wki = Sp(rki, rcMin[itype][ktype], rcMax[itype][ktype]);
            Nki = nC[atomk] - ( wki * kronecker(itype, 0) ) +
                  nH[atomk] - ( wki * kronecker(itype, 1) );
            cosjik = ( ( rji_vec[0] * rki_vec[0] ) + ( rji_vec[1] * rki_vec[1] ) + ( rji_vec[2] * rki_vec[2] ) ) / ( rji * rki );
            cosjik = min(cosjik, 1.0);
            cosjik = max(cosjik, -1.0);

            g = gSpline(cosjik, ( NijC + NijH ), itype);
            Etmp += ( wki * g * exp(lamdajik) );
            NconjtmpI = NconjtmpI + ( kronecker(ktype, 0) * wki * Sp(Nki, Nmin, Nmax) );
        }
    }

    PijS = 0.0;
    PijS = PijSpline(NijC, NijH, itype, jtype);
    pij = pow( 1.0 + Etmp + PijS, -0.5);

    Etmp = 0.0;

    REBO_neighbours_num_j = REBO_neighbours_num[j];
    REBO_neighbours_list_j = REBO_neighbours_list[j];
    REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];
    for (l = 0; l < REBO_neighbours_num_j; l++)
    {
        atoml = REBO_neighbours_list_j[l];
        if ( atoml != atomi )
        {
            ltype = type[atoml];
            rlj_vec[0] = REBO_neighbours_bonds_j[l].x;
            rlj_vec[1] = REBO_neighbours_bonds_j[l].y;
            rlj_vec[2] = REBO_neighbours_bonds_j[l].z;
            rlj = REBO_neighbours_bonds_j[l].r;

            lamdaijl = 4.0 * kronecker(jtype, 1) *
                       ( ( rho[ltype][1] - rlj ) - ( rho[itype][1] - rji ) );
            wlj = Sp(rlj, rcMin[jtype][ltype], rcMax[jtype][ltype]);
            Nlj = nC[atoml] - ( wlj * kronecker(jtype, 0) ) +
                  nH[atoml] - ( wlj * kronecker(jtype, 1) );
            cosijl = -1.0 * ( ( rji_vec[0] * rlj_vec[0] ) + ( rji_vec[1] * rlj_vec[1] ) + ( rji_vec[2] * rlj_vec[2] ) ) / ( rji * rlj );
            cosijl = min(cosijl, 1.0);
            cosijl = max(cosijl, -1.0);

            g = gSpline(cosijl, NjiC + NjiH, jtype);
            Etmp += (wlj * g * exp(lamdaijl) );
            NconjtmpJ = NconjtmpJ + ( kronecker(ltype, 0) * wlj * Sp(Nlj, Nmin, Nmax) );
        }
    }

    PjiS = 0.0;
    PjiS = PijSpline(NjiC, NjiH, jtype, itype);
    pji = pow( 1.0 + Etmp + PjiS, -0.5);

    Nijconj = 1.0 + ( NconjtmpI * NconjtmpI ) + ( NconjtmpJ * NconjtmpJ );
    piRC = piRCSpline(NijC + NijH, NjiC + NjiH, Nijconj, itype, jtype);
    Tij = 0.0;
    if ( ( itype == 0 ) && ( jtype == 0) )
        Tij = TijSpline(NijC + NijH, NjiC + NjiH, Nijconj);

    Etmp = 0.0;
    if ( fabs(Tij) > TOL )
    {
        REBO_neighbours_num_i = REBO_neighbours_num[i];
        REBO_neighbours_list_i = REBO_neighbours_list[i];
        REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];
        for (k = 0; k < REBO_neighbours_num_i; k++)
        {
            atomk = REBO_neighbours_list_i[k];
            ktype = type[atomk];
            if ( atomk != atomj )
            {
                rki_vec[0] = REBO_neighbours_bonds_i[k].x;
                rki_vec[1] = REBO_neighbours_bonds_i[k].y;
                rki_vec[2] = REBO_neighbours_bonds_i[k].z;
                rki = REBO_neighbours_bonds_i[k].r;
                cos321 = ( ( rji_vec[0] * rki_vec[0] ) + ( rji_vec[1] * rki_vec[1] ) + ( rji_vec[2] * rki_vec[2] ) ) / ( rji * rki );
                cos321 = min(cos321, 1.0);
                cos321 = max(cos321, -1.0);

                rkj_vec[0] = rki_vec[0] - rji_vec[0];
                rkj_vec[1] = rki_vec[1] - rji_vec[1];
                rkj_vec[2] = rki_vec[2] - rji_vec[2];
                rkj_sq = ( rkj_vec[0] * rkj_vec[0] ) + ( rkj_vec[1] * rkj_vec[1] ) + ( rkj_vec[2] * rkj_vec[2] );

                rji_sq = rji * rji;
                rki_sq = rki * rki;
                costmp = 0.5 * ( rji_sq + rki_sq - rkj_sq ) / rji / rki;
                tspjik = Sp2(costmp, thmin, thmax);

                if ( sqrt(1.0 - cos321 * cos321) > sqrt(TOL) )
                {
                    wki = Sp(rki, rcMin[itype][ktype], rcMaxP[itype][ktype]);

                    REBO_neighbours_num_j = REBO_neighbours_num[j];
                    REBO_neighbours_list_j = REBO_neighbours_list[j];
                    REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];
                    for (l = 0; l < REBO_neighbours_num_j; l++)
                    {
                        atoml = REBO_neighbours_list_j[l];
                        ltype = type[atoml];
                        if ( !( (atoml == atomi) || ( atoml == atomk ) ) )
                        {
                            rlj_vec[0] = REBO_neighbours_bonds_j[l].x;
                            rlj_vec[1] = REBO_neighbours_bonds_j[l].y;
                            rlj_vec[2] = REBO_neighbours_bonds_j[l].z;
                            rlj = REBO_neighbours_bonds_j[l].r;

                            cos234 = - ( ( rji_vec[0] * rlj_vec[0] ) + ( rji_vec[1] * rlj_vec[1] ) + ( rji_vec[2] * rlj_vec[2] ) ) / ( rji * rlj );
                            cos234 = min(cos234, 1.0);
                            cos234 = max(cos234, -1.0);

                            rli_vec[0] = rji_vec[0] + rlj_vec[0];
                            rli_vec[1] = rji_vec[1] + rlj_vec[1];
                            rli_vec[2] = rji_vec[2] + rlj_vec[2];

                            rli_sq = ( rli_vec[0] * rli_vec[0] ) + ( rli_vec[1] * rli_vec[1] ) + ( rli_vec[2] * rli_vec[2] );
                            rlj_sq = rlj * rlj;
                            costmp = 0.5 * ( rji_sq + rlj_sq - rli_sq ) / rji / rlj;
                            tspijl = Sp2(costmp, thmin, thmax);

                            if ( sqrt(1.0 - cos234 * cos234) > sqrt(TOL) )
                            {
                                wlj = Sp(rlj, rcMin[jtype][ltype], rcMaxP[jtype][ltype]);
                                crosskij[0] = ( rji_vec[1] * rki_vec[2] - rji_vec[2] * rki_vec[1] );
                                crosskij[1] = ( rji_vec[2] * rki_vec[0] - rji_vec[0] * rki_vec[2] );
                                crosskij[2] = ( rji_vec[0] * rki_vec[1] - rji_vec[1] * rki_vec[0] );
                                crosskijmag = sqrt( crosskij[0] * crosskij[0] + crosskij[1] * crosskij[1] + crosskij[2] * crosskij[2] );
                                crossijl[0] = ( rji_vec[1] * rlj_vec[2] - rji_vec[2] * rlj_vec[1] );
                                crossijl[1] = ( rji_vec[2] * rlj_vec[0] - rji_vec[0] * rlj_vec[2] );
                                crossijl[2] = ( rji_vec[0] * rlj_vec[1] - rji_vec[1] * rlj_vec[0] );
                                crossijlmag = sqrt( crossijl[0] * crossijl[0] + crossijl[1] * crossijl[1] + crossijl[2] * crossijl[2] );
                                omkijl = -1.0 * ( ( ( crosskij[0] * crossijl[0] ) + ( crosskij[1] * crossijl[1] ) + ( crosskij[2] * crossijl[2] ) ) / ( crosskijmag * crossijlmag) );
                                Etmp += ( ( 1.0 - square(omkijl) ) * wki * wlj ) * ( 1.0 - tspjik ) * ( 1.0 - tspijl );
                            }
                        }
                    }
                }
            }
        }
    }

    bij = ( 0.5 * ( pij + pji ) ) + piRC + ( Tij * Etmp );
    Stb = Sp2(bij, bLJmin[itype][jtype], bLJmax[itype][jtype]);

    return Stb;
}
