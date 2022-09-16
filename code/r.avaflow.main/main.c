/*
##############################################################################
#
# MODULE:       r.avaflow.main.c
#
# AUTHORS:      Martin Mergili and Shiva P. Pudasaini
# CONTRIBUTORS: Massimiliano Alvioli, Wolfgang Fellin, Jan-Thomas Fischer,
#               Sigridur S. Gylfadottir, Markus Metz, Markus Neteler, 
#               Alexander Ostermann, Matthias Rauter
#
# PURPOSE:      The simulation model for avalanche and debris flows
#               Flow routing script
#
# COPYRIGHT:    (c) 2008 - 2021 by the authors
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2013 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 2014 - 2021 by the University of Bonn
#               (c) 2000 - 2021 by the GRASS Development Team
#
# VERSION:      20210728 (28 July 2021)
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
# OPEN ISSUES   !!! Results seem to be symmetric around each axis, but for 
#               highly viscous flows, behaviour is slightly different in 
#               diagonal directions compared to the directions along the axes.
#
#               !!! Models for enhanced gravity, dispersion,  
#               dynamic adaptation of friction parameters, and diffusion  
#               control are in an experimental stage and do not necessarily
#               yield the expected results.
#
##############################################################################
*/


#define WITHGRASS // use of GRASS GIS (REMOVE THIS LINE IF GRASS IS NOT USED)


#include <fcntl.h> // libraries
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>


#ifdef WITHGRASS


    #include <grass/gmath.h>
    #include <grass/gis.h>
    #include <grass/gprojects.h>
    #include <grass/glocale.h>
    #include <grass/raster.h>
    #include <grass/segment.h>


#endif


struct ico { // global structure for input constants and control variables

    char *MAINMAPSET; int MULT; int LIMITER; int MODEL; int M; int N; int IMAX; int ENTR; int ENTR2; int ENTR3; int ZONES; int CENTR; int CVSHEAR;
    int RELM; int RELM2; int RELM3; int RELV; int RELV2; int RELV3; int PHI; int PHI2; int PHI3; int DELTAB; int TUFRI; int DELTA; int DELTA2; int DELTA3; int NYSS; int NYFS; 
    int NYFF; int AMBDRAG; int FLUFRI; int TRANSSSFS; int TRANSSSFF; int TRANSFSFF; int TRELEASE; int TRELSTOP; int STOPTIME; int TSLIDE; float PI; float UNDEF; float HFLOWMIN; 
    float CSZ; float BDWEST; float BDNORTH; float BDSOUTH; float GRAVITY; float CFL[2]; float IMPTHR[3]; int CORRHEIGHT; int CURVCTRL; int NOOSC;
    int ENTRAINMENT; int STOPPING; int BETACTRL; int NVECTMIN; int PMAX; int PHASES[3]; float EKINCOEF; float FRICMIN; float FRI_exp; float DYNFRIC; int MESH; float SLOMO;
    int HYDADD; int ORIGINAL; int SEPFLUX; int COLLAPSE; int AFLAG; float SLIDERAD; float SLIDEEXP; float SLIDEDEF;


    #ifdef WITHGRASS


        float NSEGRC; float NSEGS;


    #endif


};

struct flow { // global structure for flow parameters

    float PHI[3]; float PHI0[3]; float RHO1; float RHO2; float DELTA[3]; float DELTA0[3]; float REP; int DRAG_je; float NY0[3]; float TAUY[3]; float AMBDRAG0; float AMBDRAG; 
    float FLUFRI0; float FLUFRI; float RHO3; float J; float NY[3]; float TUFRI0; float TUFRI; float TRANSSSFS0; float TRANSSSFF0; 
    float TRANSFSFF0; float TRANSSSFS; float TRANSSSFF; float TRANSFSFF; float CENTR0; float CENTR; float CVSHEAR0; float CVSHEAR; float DELTAB0; float DELTAB; float THETAS; float CSTOP; 
    float VM_n0; float VM_l; float VM_n; float VM_fact; float DRAG_k; float DRAG_m; int DRAG_n; float DRAG_rep; float DRAG_vterm; 
    float VIS_chi[3]; float VIS_xi[3]; float VIS_a[3]; float VIS_ry; float NY_exp; float KPMAX[3]; float RHOB1; float RHOB2; float RHOB3; float sep_SHEAR; float disp_MULT; float CCONST;
};


// -- START -- Functions ----------------------------------------------------------------------------------------


float ffmax ( float gval1, float gval2 ) { // function for maximum value (float)

    float gmax;
    if ( gval2 > gval1 ) gmax = gval2;
    else gmax = gval1;
    return gmax;
}

float ffmin ( float gval1, float gval2 ) { // function for minimum value (float)

    float gmin;
    if ( gval2 < gval1 ) gmin = gval2;
    else gmin = gval1;
    return gmin;
}

float fabsmin ( float gval1, float gval2 ) { // function for minimum absolute value (float)

    float gmin;
    if ( fabs(gval2) < fabs(gval1) ) gmin = gval2;
    else gmin = gval1;
    return gmin;
}

float fsign ( float gval ) { // function for sign

    float gsign;
    if ( gval < 0 ) gsign = -1;
    else if ( gval == 0 ) gsign = 0;
    else gsign = 1;
    return gsign;
}

int fvalid ( float ghflow, float ghflow2, float ghflow3, struct ico sico ) { // function for minimum flow height

    int gvalid;

    gvalid = 0;
    if ( ghflow > sico.HFLOWMIN ) gvalid = 1;
    if ( sico.MODEL == 7 && ghflow2 > sico.HFLOWMIN ) gvalid = 1;
    if ( sico.MODEL == 7 && ghflow3 > sico.HFLOWMIN ) gvalid = 1;

    return gvalid;
}

int **alloc_imatrix(int mrows, int mcols) { // function for allocating 2D integer arrays
    int **mm;
    int mi;
    mm = (int **) calloc(mrows, sizeof(int *));
    mm[0] = (int *) calloc(mrows * mcols, sizeof(int));
    for (mi = 1; mi < mrows; mi++)
    mm[mi] = mm[mi - 1] + mcols;
    return mm;
}

float **alloc_dmatrix (int mrows, int mcols) { // function for allocating 2D float arrays
    float **mm;
    int mi;
    mm = (float **) calloc(mrows, sizeof(float *));
    mm[0] = (float *) calloc(mrows * mcols, sizeof(float));
    for (mi = 1; mi < mrows; mi++)
    mm[mi] = mm[mi - 1] + mcols;
    return mm;
}

char **alloc_cmatrix (int mrows, int mcols) { // function for allocating 2D char arrays
    char **mm;
    int mi;
    mm = (char **) calloc(mrows, sizeof(char *));
    mm[0] = (char *) calloc(mrows * mcols, sizeof(char));
    for (mi = 1; mi < mrows; mi++)
    mm[mi] = mm[mi - 1] + mcols;
    return mm;
}

float ***alloc_dmatrix3 (int mrows, int mcols, int mdepths) { // function for allocating 3D float arrays
    float ***mm;
    int mi, mj;
    mm = (float***)calloc(mrows, sizeof(float**));

    for (mi=0; mi<mrows; mi++) {
        mm[mi] = (float**)calloc(mcols, sizeof(float*));

        for (mj=0; mj<mcols; mj++) {
            mm[mi][mj] = (float*)calloc(mdepths, sizeof(float));
        }
    }
    return mm;
}

void free_dmatrix3 (float ***mm,int mi,int mj) { // function for freeing 3D float arrays
    int i,j;

    for(i=0;i<mi;i++)
    {
        for(j=0;j<mj;j++)
        {
            free(mm[i][j]);
        }
        free(mm[i]);
    }
    free(mm);
}

int fiparam ( FILE *gfparam ) { // function for reading integer parameters

    int gparam;
    char gsparam[50];
    char **strtodhlp=NULL;
    if ( fgets(gsparam, 50, gfparam) == NULL ) return -1;
    gparam=strtod(gsparam, strtodhlp);
    return gparam;
}

float fdparam ( FILE *gfparam ) { // function for reading float parameters

    float gparam;
    char gsparam[50];
    char **strtodhlp=NULL;
    if ( fgets(gsparam, 50, gfparam) == NULL ) return -1;
    gparam=strtod(gsparam, strtodhlp);
    return gparam;
}

char *fcparam ( FILE *gfparam ) { // function for reading character parameters

    char *gsparam = (char*) calloc ( 1000, sizeof(char));
    if ( fgets(gsparam, 1000, gfparam) == NULL ) return "-1";
    gsparam[strnlen(gsparam, 1010) - 1] = '\0';
    return gsparam;
}

char *flparam ( FILE *gfparam ) { // function for reading ascii raster line

    int asclinemax = 100000;
    char *gsparam2 = (char*) calloc (( asclinemax + 20 ), sizeof(char));
    if ( fgets(gsparam2, asclinemax, gfparam) == NULL ) return "-1";
    gsparam2[strnlen(gsparam2, asclinemax+10) - 1] = '\0';
    return gsparam2;
}

int *fin ( int gi, float *gpelev, struct ico sico ) { // function for identifying cell environment

    int *gin;

    gin = (int*) calloc( 9, sizeof(int));

    gin[0] = gi;
    if ( gi%sico.N == 0 || (gi+1)%sico.N == 0 || gi < sico.N || gi > sico.IMAX - sico.N || gpelev[gi] == sico.UNDEF ) { // edge cells (to be excluded from computation)
    
        gin[1] = sico.IMAX;
        gin[2] = sico.IMAX;
        gin[3] = sico.IMAX;
        gin[4] = -1;
        gin[5] = -1;
        gin[6] = -1;
        gin[7] = sico.IMAX;
        gin[8] = -1;
    }
    else { // cells to be included in computation
        gin[1] = gi + 1;
        gin[2] = gi + sico.N;
        gin[3] = gi + sico.N + 1;
        gin[4] = gi - 1;
        gin[5] = gi - sico.N;
        gin[6] = gi - sico.N - 1;
        gin[7] = gi + sico.N - 1;
        gin[8] = gi - sico.N + 1;
    }

    return gin;
}

float fdiv ( float gabove, float gbelow, float gcriterion ) { // function to avoid division by zero

    float gdiv;
    if ( gbelow > gcriterion ) gdiv = gabove / gbelow; else gdiv = 0;
    return gdiv;
}

float fdynfric ( float gfric, float gfricmin, float ggekin, float ggalpha, struct ico sico ) { // function for dynamically varying the frictions

    float gdynfric;
    gdynfric = ( gfricmin + ( gfric - gfricmin ) * exp( -ggekin * sico.EKINCOEF )) * pow( ggalpha, sico.FRI_exp );
    return gdynfric;
}

float fbeta ( float gelevm, float gelevp, float gspan, struct ico sico ) { // function for topographic slopes in x and y direction

    float gbeta;
    gbeta = atan(( gelevm - gelevp ) / ( gspan * sico.CSZ ));
    return gbeta;
}

float fbetaxy ( float gbetax, float gbetay ) { // function for topographic slope

    float gbetaxy;
    if ( gbetax != 0 || gbetay != 0 ) gbetaxy = atan( pow( tan ( gbetax) * tan ( gbetax) + tan( gbetay ) * tan( gbetay ), 0.5 ));
    else gbetaxy = 0;
    return gbetaxy;
}

float falpha ( float gbetax, float gbetay, struct ico sico ) { // function for topographic aspect

    float galpha = 0;

    if ( gbetax == 0 && gbetay == 0 ) galpha = 0;
    else if ( gbetax == 0 && gbetay > 0 ) galpha = sico.PI * 0.5;
    else if ( gbetax == 0 && gbetay < 0 ) galpha = 3 * sico.PI * 0.5;
    else if ( gbetax >= 0 ) galpha = atan( tan( gbetay ) / tan( gbetax ));
    else if ( gbetax < 0 ) galpha = atan( tan( gbetay ) / tan( gbetax )) + sico.PI;

    if ( galpha < 0 ) galpha += 2 * sico.PI;

    return galpha;
}

float falphav ( float gmomx, float gmomy, struct ico sico ) { // function for direction of movement

    float galphav = 0;

    if ( gmomx == 0 && gmomy == 0 ) galphav = 0;
    else if ( gmomx == 0 && gmomy > 0 ) galphav = sico.PI * 0.5;
    else if ( gmomx == 0 && gmomy < 0 ) galphav = 3 * sico.PI * 0.5;
    else if ( gmomx >= 0 ) galphav = atan( gmomy / gmomx );
    else if ( gmomx < 0 ) galphav = atan( gmomy / gmomx ) + sico.PI;

    if ( galphav < 0 ) galphav += 2 * sico.PI;

    return galphav;
}

float **finhyd ( char *gname, int ghydtmax, int ghydtmaxx ) { // function for input of hydrograph

    FILE *gfhyd;
    char gpath[200], *gtesth, *gtesth2;
    int gi, gj;
    char *tok;
    char gdelim[] = " \t\n";
    char **strtodhlp=NULL;

    float **ghyd = alloc_dmatrix( ghydtmaxx+1, 7 );

    sprintf(gpath, "%s", gname );
    gfhyd=fopen(gpath, "r");

    gtesth = fcparam ( gfhyd );
    free( gtesth );
    for ( gi=0; gi<=ghydtmax; gi++ ) {
        gtesth = fcparam ( gfhyd );
        tok=strtok(gtesth, gdelim);
        gj = 0;
        while ( tok != NULL ) {
            gtesth2=tok;
            tok=strtok(NULL, gdelim);
            ghyd[gi][gj] = strtod (gtesth2, strtodhlp);
            gj += 1;
        }
        free(gtesth);
    }

    fclose(gfhyd);
    return ghyd;
}

int *fhydp ( int ghydx, int ghydy, float ghydl, int *gpx, int *gpy, float galpha, float ghydlmax, struct ico sico ) { // function for input hydrograph profile

    int gdistx, gdisty, gi, gj, gx0, gy0, gdprex, gdprey, gctrl, *gx, *gy;
    float gd, gdxtest, gdytest;

    int *ghydp = (int*) calloc( (int)(2*ghydlmax/sico.CSZ+2), sizeof(int) );

    gdistx = (int)( ghydl * fabs( cos( galpha )) / ( 2 * sico.CSZ ) + 0.5 );
    gdisty = (int)( ghydl * fabs( sin( galpha )) / ( 2 * sico.CSZ ) + 0.5 );

    gdxtest = fabs((float)( gdistx ) / cos( galpha ));
    gdytest = fabs((float)( gdisty ) / sin( galpha ));

    gx = (int*) calloc( 4 * (( gdistx + 1 ) * ( gdisty + 1 )), sizeof(int));
    gy = (int*) calloc( 4 * (( gdistx + 1 ) * ( gdisty + 1 )), sizeof(int));

    gctrl = 1;

    if ( galpha == 0 || galpha == sico.PI * 0.5 || gdytest == 0 ) {

        gx0 = ghydx - gdistx;
        gy0 = ghydy;

        gx[gctrl] = gx0;
        gy[gctrl] = gy0;

        while ( gx[gctrl] <= ghydx + gdistx ) {

            gctrl += 1;
            gx[gctrl] = gx0 + gctrl;
            gy[gctrl] = gy0;
        }
        
    } else if ( galpha == sico.PI || galpha == 3 * sico.PI * 0.5 || gdxtest == 0 ) {

        gx0 = ghydx;
        gy0 = ghydy - gdisty;

        gx[gctrl] = gx0;
        gy[gctrl] = gy0;

        while ( gy[gctrl] <= ghydy + gdisty ) {

            gctrl += 1;
            gx[gctrl] = gx0;
            gy[gctrl] = gy0 + gctrl;
        }

    } else {

        if ( gdxtest > gdytest ) {

            gx0 = ghydx - gdistx;

            if (( galpha > 0 && galpha < sico.PI * 0.5 ) || ( galpha > sico.PI && galpha < 3 * sico.PI * 0.5 ))
                gy0 = ghydy - gdisty;
            else gy0 = ghydy + gdisty;

            gdprex = 1;
            if (( galpha > 0 && galpha < sico.PI * 0.5 ) || ( galpha > sico.PI && galpha < 3 * sico.PI * 0.5 )) gdprey = 1;
            else gdprey = -1;
            
        } else {

            if (( galpha > 0 && galpha < sico.PI * 0.5 ) || ( galpha > sico.PI && galpha < 3 * sico.PI * 0.5 ))
                gx0 = ghydx - gdistx;
            else gx0 = ghydx + gdistx;
            gy0 = ghydy - gdisty;
            if (( galpha > 0 && galpha < sico.PI * 0.5 ) || ( galpha > sico.PI && galpha < 3 * sico.PI * 0.5 )) gdprex = 1;
            else gdprex = -1;
            gdprey = 1;
        }
        
        if ( fabs( 1 / cos( galpha )) < fabs(1 / sin ( galpha ))) gd = fabs( 1 / cos( galpha )); else gd = fabs(1 / sin ( galpha ));

        gx[gctrl] = gx0;
        gy[gctrl] = gy0;

        while ( gx[gctrl] <= ghydx + gdistx && gy[gctrl] <= ghydy + gdisty
            && gx[gctrl] >= ghydx - gdistx && gy[gctrl] >= ghydy - gdisty ) {

            gctrl += 1;
            gx[gctrl] = gx0 + gdprex * (int)( gctrl * gd * fabs(cos( galpha )) + 0.5 );
            gy[gctrl] = gy0 + gdprey * (int)( gctrl * gd * fabs(sin( galpha )) + 0.5 );
        }
    }
    
    gctrl -= 1;
    ghydp[0] = gctrl;

    for ( gi=1; gi<sico.IMAX; gi++ ) {
        for ( gj=1; gj<=gctrl; gj++ ) {
            if ( gpx[gi] == gx[gj] && gpy[gi] == gy[gj] ) {
                ghydp[gj] = gi;
            }
        }
    }

    free(gx); free(gy);

    return ghydp;
}

float fconvin ( float ghflow, float gbetaxy, struct ico sico ) { // function for converting heights into depths

    float gdflow;
    if ( sico.CORRHEIGHT == 0 ) gdflow = ghflow; 
    else gdflow = ghflow * cos ( gbetaxy );
    return gdflow;
}

float fk ( float gh, float ghflow, float gvflowx, float gvflowy, float gdvflow, float gekin, int gp, struct ico sico, struct flow sflow ) { // function for earth pressure coefficients

    float ggk, gka, gkp, galpha, gdelta, gphi;

    if ( ghflow > sico.HFLOWMIN ) galpha = ghflow / gh; else galpha = 0;
    
    if ( gp == 9 ) {
        
        if ( sico.DYNFRIC == 1 ) {
        
            gdelta = fdynfric( gvflowx * sflow.DELTA[0] + gvflowy * sflow.DELTA[1], sico.FRICMIN, gekin, galpha, sico );
            gphi = fdynfric( gvflowx * sflow.PHI[0] + gvflowy * sflow.PHI[1], sico.FRICMIN, gekin, galpha, sico );
            
        } else {
        
            gdelta = gvflowx * sflow.DELTA[0] + gvflowy * sflow.DELTA[1];
            gphi = gvflowx * sflow.PHI[0] + gvflowy * sflow.PHI[1];         
        }
        gp = 0;      
        
    } else {
    
        if ( sico.DYNFRIC == 1 ) gdelta = fdynfric( sflow.DELTA[gp], sico.FRICMIN, gekin, galpha, sico ); else gdelta = sflow.DELTA[gp];
        if ( sico.DYNFRIC == 1 ) gphi = fdynfric( sflow.PHI[gp], sico.FRICMIN, gekin, galpha, sico );else gphi = sflow.PHI[gp];
    }

    if ( gphi < gdelta ) gphi = gdelta;
    if ( gphi == sico.UNDEF ) gphi = gdelta; // friction angles

    gka = 2 * ( 1 - pow( 1 - pow( cos( gphi ) / cos( gdelta ), 2 ), 0.5 )) / pow( cos( gphi ), 2 ) - 1; // active earth pressure coefficient
    gkp = 2 * ( 1 + pow( 1 - pow( cos( gphi ) / cos( gdelta ), 2 ), 0.5 )) / pow( cos( gphi ), 2 ) - 1; // passive eath pressure coefficient

    if ( gdvflow >= 0 ) ggk=gka;
    else ggk = ffmin(sflow.KPMAX[gp], gkp ); // earth pressure coefficient

    return ggk;
}

float *fcurv ( float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, float gvflowy2, float gvflowy3, float *gelev, float *ggrav, struct ico sico ) { // function for curvature

    int gl;
    float gu[3], gv[3], gkappax, gkappaxy, gkappay, *gkappau;
    gkappau = (float*) calloc( sico.PMAX, sizeof(float));

    gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
    gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

    gkappax = 8 * sico.CSZ * ( gelev[2] + gelev[5] - 2 * gelev[0] )
        / (( 4 * pow ( sico.CSZ, 2 ) + pow( gelev[2] - gelev[5], 2 ))
        * (pow( 4 * pow ( sico.CSZ, 2 ) + pow( gelev[2] - gelev[5] , 2 ) + pow( gelev[1] - gelev[4] , 2 ), 0.5 )));  

    gkappaxy = 2 * sico.CSZ * ( gelev[3] + gelev[6] - gelev[8] - gelev[7] )
	/ ( pow( 4 * pow ( sico.CSZ, 2 ) + pow( gelev[2] - gelev[5], 2 ), 0.5 ) * pow( 4 * pow ( sico.CSZ, 2 ) 
        + pow( gelev[1] - gelev[4] , 2 ), 0.5 ) * pow( 4 * sico.CSZ * sico.CSZ + pow( gelev[2] - gelev[5] , 2 ) 
        + pow( gelev[1] - gelev[4] , 2 ), 0.5 ));

    gkappay = 8 * sico.CSZ * ( gelev[1] + gelev[4] - 2 * gelev[0] )
	/ (( 4 * pow ( sico.CSZ, 2 ) + pow( gelev[1] - gelev[4], 2 )) 
        * (pow( 4 * pow ( sico.CSZ, 2 ) + pow( gelev[2] - gelev[5] , 2 ) + pow( gelev[1] - gelev[4] , 2 ), 0.5 ))); 

    for ( gl = 0; gl < sico.PMAX; gl++ ) {

        if ( gu[gl] != 0 || gv[gl] != 0 ) gkappau[gl] = ( pow( gu[gl], 2 ) * gkappax + 2 * gu[gl] * gv[gl] * gkappaxy + pow( gv[gl], 2 ) * gkappay );
        else gkappau[gl] = 0;
    }

    return gkappau;
}

float *fvm ( float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, float gvflowy2, float gvflowy3,
    struct ico sico, struct flow sflow ) { // function for virtual mass

    float *gvm, gh, galpha[3], ggamma[3], gu[3], gv[3];
    gvm = (float*) calloc( 21, sizeof(float));

    if ( sico.MODEL == 7 ) {

        gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
        gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

        gh = ghflow + ghflow2 + ghflow3;
        if ( gh > sico.HFLOWMIN ) {

            if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
            if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
            if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

        if ( galpha[0] > 0  && galpha[2] > 0 ) ggamma[0] = sflow.RHO3 / sflow.RHO1; else ggamma[0] = 0;
        if ( galpha[0] > 0  && galpha[1] > 0 ) ggamma[1] = sflow.RHO2 / sflow.RHO1; else ggamma[1] = 0;
        if ( galpha[1] > 0  && galpha[2] > 0 ) ggamma[2] = sflow.RHO3 / sflow.RHO2; else ggamma[2] = 0;

        if ( galpha[2] > 0 && galpha[0] > 0 )
            gvm[9] = ( sflow.VM_n0 * ( sflow.VM_l + pow( galpha[0], sflow.VM_n )) - 1 ) / ( galpha[0] / galpha[2] + ggamma[0] ); else gvm[9] = 0;
        if ( galpha[1] > 0 && galpha[0] > 0 )
            gvm[10] = ( sflow.VM_n0 * ( sflow.VM_l + pow( galpha[0], sflow.VM_n )) - 1 ) / ( galpha[0] / galpha[1] + ggamma[1] ); else gvm[10] = 0;
        if ( galpha[2] > 0 && galpha[1] > 0 )
            gvm[11] = ( sflow.VM_n0 * ( sflow.VM_l + pow( galpha[1], sflow.VM_n )) - 1 ) / ( galpha[1] / galpha[2] + ggamma[2] ); else gvm[11] = 0;

        if ( galpha[0] > 0 ) {

            gvm[0] = ggamma[0] * gvm[9] * ( gu[2] - gu[0] ) + ggamma[1] * gvm[10] * ( gu[1] - gu[0] );
            gvm[1] = ggamma[0] * gvm[9] * ( gu[2] * gu[2] - gu[0] * gu[0]) + ggamma[1] * gvm[10] * ( gu[1] * gu[1] - gu[0] * gu[0]);
            gvm[2] = ggamma[0] * gvm[9] * ( gu[2] * gv[2] - gu[0] * gv[0] ) + ggamma[1] * gvm[10] * ( gu[1] * gv[1] - gu[0] * gv[0] );

            if ( fabs( gvm[1] ) > sflow.VM_fact * fabs( gu[2] - gu[0] ) && fabs( gvm[1] ) > sflow.VM_fact * fabs( gu[1] - gu[0] ))
                gvm[1] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[0] ), fabs( gu[1] - gu[0] )));
            if ( fabs( gvm[2] ) > sflow.VM_fact * fabs( gv[2] - gv[0] ) && fabs( gvm[2] ) > sflow.VM_fact * fabs( gv[1] - gv[0] ))
                gvm[2] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[0] ), fabs( gv[1] - gv[0] )));

            gvm[12] = ggamma[0] * gvm[9] * ( gv[2] - gv[0] ) + ggamma[1] * gvm[10] * ( gv[1] - gv[0] );
            gvm[13] = ggamma[0] * gvm[9] * ( gv[2] * gv[2] - gv[0] * gv[0]) + ggamma[1] * gvm[10] * ( gv[1] * gv[1] - gv[0] * gv[0]);
            gvm[14] = ggamma[0] * gvm[9] * ( gu[2] * gv[2] - gu[0] * gv[0] ) + ggamma[1] * gvm[10] * ( gu[1] * gv[1] - gu[0] * gv[0] );

            if ( fabs( gvm[13] ) > sflow.VM_fact * fabs( gv[2] - gv[0] ) && fabs( gvm[13] ) > sflow.VM_fact * fabs( gv[1] - gv[0] ))
                gvm[13] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[0] ), fabs( gv[1] - gv[0] )));
            if ( fabs( gvm[14] ) > sflow.VM_fact * fabs( gu[2] - gu[0] ) && fabs( gvm[14] ) > sflow.VM_fact * fabs( gu[1] - gu[0] ))
                gvm[14] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[0] ), fabs( gu[1] - gu[0] )));

        }
        else { gvm[0] = 0; gvm[1] = 0; gvm[2] = 0; gvm[12] = 0; gvm[13] = 0; gvm[14] = 0; }

        if ( galpha[1] > 0 ) {

            gvm[3] = ggamma[2] * gvm[11] * ( gu[2] - gu[1] ) - galpha[0] / galpha[1] * gvm[10] * ( gu[1] - gu[0] );
            gvm[4] = ggamma[2] * gvm[11] * ( gu[2] * gu[2] - gu[1] * gu[1]) - galpha[0] / galpha[1] * gvm[10] * ( gu[1] * gu[1] - gu[0] * gu[0]);
            gvm[5] = ggamma[2] * gvm[11] * ( gu[2] * gv[2] - gu[1] * gv[1] ) - galpha[0] / galpha[1] * gvm[10] * ( gu[1] * gv[1] - gu[0] * gv[0] );

            if ( fabs( gvm[4] ) > sflow.VM_fact * fabs( gu[2] - gu[1] ) && fabs( gvm[4] ) > sflow.VM_fact * fabs( gu[0] - gu[1] ))
                gvm[4] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[1] ), fabs( gu[0] - gu[1] )));
            if ( fabs( gvm[5] ) > sflow.VM_fact * fabs( gv[2] - gv[1] ) && fabs( gvm[5] ) > sflow.VM_fact * fabs( gv[0] - gv[1] ))
                gvm[5] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[1] ), fabs( gv[0] - gv[1] )));

            gvm[15] = ggamma[2] * gvm[11] * ( gv[2] - gv[1] ) - galpha[0] / galpha[1] * gvm[10] * ( gv[1] - gv[0] );
            gvm[16] = ggamma[2] * gvm[11] * ( gv[2] * gv[2] - gv[1] * gv[1]) - galpha[0] / galpha[1] * gvm[10] * ( gv[1] * gv[1] - gv[0] * gv[0]);
            gvm[17] = ggamma[2] * gvm[11] * ( gu[2] * gv[2] - gu[1] * gv[1] ) - galpha[0] / galpha[1] * gvm[10] * ( gu[1] * gv[1] - gu[0] * gv[0] );

            if ( fabs( gvm[16] ) > sflow.VM_fact * fabs( gv[2] - gv[1] ) && fabs( gvm[16] ) > sflow.VM_fact * fabs( gv[0] - gv[1] ))
                gvm[16] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[1] ), fabs( gv[0] - gv[1] )));
            if ( fabs( gvm[17] ) > sflow.VM_fact * fabs( gu[2] - gu[1] ) && fabs( gvm[17] ) > sflow.VM_fact * fabs( gu[0] - gu[1] ))
                gvm[17] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[1] ), fabs( gu[0] - gu[1] )));

        }
        else { gvm[3] = 0; gvm[4] = 0; gvm[5] = 0; gvm[15] = 0; gvm[16] = 0; gvm[17] = 0; }

        if ( galpha[2] > 0 ) {

            gvm[6] = galpha[0] / galpha[2] * gvm[9] * ( gu[2] - gu[0] ) + galpha[1] / galpha[2] * gvm[11] * ( gu[2] - gu[1] );
            gvm[7] = galpha[0] / galpha[2] * gvm[9] * ( gu[2] * gu[2] - gu[0] * gu[0]) + galpha[1] / galpha[2] * gvm[11] * ( gu[2] * gu[2] - gu[1] * gu[1]);
            gvm[8] = galpha[0] / galpha[2] * gvm[9] * ( gu[2] * gv[2] - gu[0] * gv[0] ) + galpha[1] / galpha[2] * gvm[11] * ( gu[2] * gv[2] - gu[1] * gv[1] );

            if ( fabs( gvm[7] ) > sflow.VM_fact * fabs( gu[2] - gu[0] ) && fabs( gvm[7] ) > sflow.VM_fact * fabs( gu[2] - gu[1] ))
                gvm[7] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[0] ), fabs( gu[2] - gu[1] )));
            if ( fabs( gvm[8] ) > sflow.VM_fact * fabs( gv[2] - gv[0] ) && fabs( gvm[8] ) > sflow.VM_fact * fabs( gv[2] - gv[1] ))
                gvm[8] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[0] ), fabs( gv[2] - gv[1] )));

            gvm[18] = galpha[0] / galpha[2] * gvm[9] * ( gv[2] - gv[0] ) + galpha[1] / galpha[2] * gvm[11] * ( gv[2] - gv[1] );
            gvm[19] = galpha[0] / galpha[2] * gvm[9] * ( gv[2] * gv[2] - gv[0] * gv[0]) + galpha[1] / galpha[2] * gvm[11] * ( gv[2] * gv[2] - gv[1] * gv[1]);
            gvm[20] = galpha[0] / galpha[2] * gvm[9] * ( gu[2] * gv[2] - gu[0] * gv[0] ) + galpha[1] / galpha[2] * gvm[11] * ( gu[2] * gv[2] - gu[1] * gv[1] );

            if ( fabs( gvm[19] ) > sflow.VM_fact * fabs( gv[2] - gv[0] ) && fabs( gvm[19] ) > sflow.VM_fact * fabs( gv[2] - gv[1] ))
                gvm[19] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[0] ), fabs( gv[2] - gv[1] )));
            if ( fabs( gvm[20] ) > sflow.VM_fact * fabs( gu[2] - gu[0] ) && fabs( gvm[20] ) > sflow.VM_fact * fabs( gu[2] - gu[1] ))
                gvm[20] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[0] ), fabs( gu[2] - gu[1] )));
        }
        else { gvm[6] = 0; gvm[7] = 0; gvm[8] = 0; gvm[18] = 0; gvm[19] = 0; gvm[20] = 0; }
    }

    return gvm;
}

float *fdrag ( float ghflow, float ghflow2, float ghflow3, struct ico sico, struct flow sflow ) { // function for drag

    float *gcdrag, gh, galpha[3], galphac[3], ggamma[3], gf, gg, gp, gsp;
    gcdrag = (float*) calloc( 3, sizeof(float));

    if ( sico.MODEL == 7 ) {

        gh = ghflow + ghflow2 + ghflow3;
        if ( gh > sico.HFLOWMIN ) {

            if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
            if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
            if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

        if ( galpha[0] > 0 && galpha[2] > 0 ) ggamma[0] = sflow.RHO3 / sflow.RHO1; else ggamma[0] = 0;
        if ( galpha[0] > 0 && galpha[1] > 0 ) ggamma[1] = sflow.RHO2 / sflow.RHO1; else ggamma[1] = 0;
        if ( galpha[1] > 0 && galpha[2] > 0 ) ggamma[2] = sflow.RHO3 / sflow.RHO2; else ggamma[2] = 0;

        galphac[0] = galpha[0] * galpha[2];
        galphac[1] = galpha[0] * galpha[1];
        galphac[2] = galpha[1] * galpha[2];

        if ( galpha[0] > 0 && galpha[2] > 0 ) {

            gf = ggamma[0] / 180 * pow( galpha[2] / galpha[0], 3 ) * sflow.DRAG_rep;
            gg = pow( galpha[2], sflow.DRAG_m - 1 );
            gp = pow( galpha[0], sflow.DRAG_n );
            gsp = ( gp / galpha[0] + ( 1 - gp ) / galpha[2] ) * sflow.DRAG_k;
            gcdrag[0] = pow( 1, sflow.DRAG_je )
                * galphac[0] * ( 1 - ggamma[0] ) * sico.GRAVITY / pow( sflow.DRAG_vterm * ( gp * gf + ( 1 - gp ) * gg ) + gsp, sflow.DRAG_je );

        } else gcdrag[0] = 0;

        if ( galpha[0] > 0 && galpha[1] > 0 ) {

            gf = ggamma[0] / 180 * pow( galpha[1] / galpha[0], 3 ) * sflow.DRAG_rep;
            gg = pow( galpha[1], sflow.DRAG_m - 1 );
            gp = pow( galpha[0], sflow.DRAG_n );
            gsp = ( gp / galpha[0] + ( 1 - gp ) / galpha[1] ) * sflow.DRAG_k;
            gcdrag[1] = pow( 1, sflow.DRAG_je )
                * galphac[1] * ( 1 - ggamma[1] ) * sico.GRAVITY / pow( sflow.DRAG_vterm * ( gp * gf + ( 1 - gp ) * gg ) + gsp, sflow.DRAG_je );

        } else gcdrag[1] = 0;

        if ( galpha[1] > 0 && galpha[2] > 0 ) {

            gf = ggamma[1] / 180 * pow( galpha[2] / galpha[1], 3 ) * sflow.DRAG_rep;
            gg = pow( galpha[2], sflow.DRAG_m - 1 );
            gp = pow( galpha[1], sflow.DRAG_n );
            gsp = ( gp / galpha[1] + ( 1 - gp ) / galpha[2] ) * sflow.DRAG_k;
            gcdrag[2] = pow( 1, sflow.DRAG_je )
                * galphac[2] * ( 1 - ggamma[2] ) * sico.GRAVITY / pow( sflow.DRAG_vterm * ( gp * gf + ( 1 - gp ) * gg ) + gsp, sflow.DRAG_je );

        } else gcdrag[2] = 0;

    } else { gcdrag[0] = 0; gcdrag[1] = 0; gcdrag[2] = 0; }

    return gcdrag;
}

float *fgze ( float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, float gvflowy2, float gvflowy3, 
    float *gnh, float gghx, float gghy, float *ggux, float *ggvy, float *gguy, float *ggvx, float *gnuss, float *gnvss, float *gnufs, float *gnvfs, float *gnuff, float *gnvff, 
    float *ggalpha10, float *ggalpha20, float *ggalpha30, float gdx, float gdy, float *ggrav, float gbetax, float gbetay, float *gcdrag, float gtsum, 
    float *gkappau, struct ico sico, struct flow sflow ) { // function for enhanced gravity

    int gb, gj, gl, gm, go;
    float *ggze, gh, gu[3], gv[3], gnu[3][9], gnv[3][9], ggalpha[3][9], gw[3], gww[3], guvw[3], galpha[3], ggamma[3], ggamma2[3], gdrag[3], gght, ggtlength[3], gtest, gtestx[3], gtesty[3];
    ggze = (float*) calloc( 12, sizeof(float));

    for ( gl = 0; gl < sico.PMAX; gl++ ) {

        if ( sico.CURVCTRL == 2 || sico.CURVCTRL == 5 ) gtest = ffmax( sico.GRAVITY, ffmin( sico.GRAVITY + gkappau[gl], sico.GRAVITY * 10 )); // curvature effects
        else gtest = sico.GRAVITY;
        
        gtestx[gl] = gtest * cos( gbetax ); gtesty[gl] = gtest * cos( gbetay );
    }

    if ( sico.MODEL == 7 ) { // multi-phase model (gravity enhanced for buoyancy and curvature)
    
        gh = ghflow + ghflow2 + ghflow3; // total flow height

        if ( gh > sico.HFLOWMIN ) {

            if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0; // phase fractions
            if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
            if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

        if ( galpha[0] > 0 && galpha[2] > 0 ) ggamma[0] = sflow.RHO3 / sflow.RHO1; else ggamma[0] = 0; // density ratios
        if ( galpha[0] > 0 && galpha[1] > 0 ) ggamma[1] = sflow.RHO2 / sflow.RHO1; else ggamma[1] = 0;
        if ( galpha[1] > 0 && galpha[2] > 0 ) ggamma[2] = sflow.RHO3 / sflow.RHO2; else ggamma[2] = 0;

        if ( galpha[0] > 0 && galpha[2] > 0 ) ggamma2[0] = sflow.RHO3 / sflow.RHO1; else ggamma2[0] = 1;
        if ( galpha[0] > 0 && galpha[1] > 0 ) ggamma2[1] = sflow.RHO2 / sflow.RHO1; else ggamma2[1] = 1;
        if ( galpha[1] > 0 && galpha[2] > 0 ) ggamma2[2] = sflow.RHO3 / sflow.RHO2; else ggamma2[2] = 1;

        ggze[0] = gtestx[0] * ( 1 - ggamma[0] ); // gz in x direction, solid PHASE 1
        ggze[3] = gtestx[0] * ggamma2[0]; // gz in x direction, solid PHASE 1 (for gh/gx)
        ggze[4] = gtesty[0] * ( 1 - ggamma[0] ); // gz in y direction, solid PHASE 1
        ggze[7] = gtesty[0] * ggamma2[0]; // gz in y direction, solid PHASE 1 (for gh/gy)
     
        ggze[8] = gtestx[1] * ( 1 - ggamma[2] ); // gz in x direction, solid PHASE 2
        ggze[9] = gtestx[1] * ggamma2[2]; // gz in x direction, solid PHASE 2 (for gh/gx)
        ggze[10] = gtesty[1] * ( 1 - ggamma[2] ); // gz in y direction, solid PHASE 2
        ggze[11] = gtesty[1] * ggamma2[2]; // gz in y direction, solid PHASE 2 (for gh/gy)

        ggze[1] = gtestx[1] * ggamma2[2]; // gz in x direction, fine solid PHASE 2
        ggze[5] = gtesty[1] * ggamma2[2]; // gz in y direction, fine solid PHASE 2

        ggze[2] = gtestx[2]; // gz in x direction, fluid PHASE 3
        ggze[6] = gtesty[2]; // gz in y direction, fluid PHASE 3
        
        if ( sico.NOOSC == 2 || sico.NOOSC == 4 || sico.NOOSC == 5 || sico.NOOSC == 7 ) { // advanced enhanced gravity model
        
            for( gj=0; gj<9; gj++ ) { 

                gnu[0][gj] = gnuss[gj]; gnu[1][gj] = gnufs[gj]; gnu[2][gj] = gnuff[gj]; gnv[0][gj] = gnvss[gj]; gnv[1][gj] = gnvfs[gj]; gnv[2][gj] = gnvff[gj];
                ggalpha[0][gj] = ggalpha10[gj]; ggalpha[1][gj] = ggalpha20[gj]; ggalpha[2][gj] = ggalpha30[gj];
            }

            gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
            gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

            gght = -0.5 / gdx * ( gnh[2] * ( ggalpha[0][2] * gnu[0][2] + ggalpha[1][2] * gnu[1][2] + ggalpha[2][2] * gnu[2][2] ) - gnh[5] * ( ggalpha[0][5] * gnu[0][5] 
                + ggalpha[1][5] * gnu[1][5] + ggalpha[2][5] * gnu[2][5] )) - 0.5 / gdy * ( gnh[1] * ( ggalpha[0][1] * gnv[0][1] + ggalpha[1][1] * gnv[1][1] 
                + ggalpha[2][1] * gnv[2][1] ) - gnh[4] * ( ggalpha[0][4] * gnv[0][4] + ggalpha[1][4] * gnv[1][4] + ggalpha[2][4] * gnv[2][4] ));

            for ( go = 0; go < 3; go++ ) {

                ggtlength[go] = -0.5 * ( -gh * ( pow(ggux[go], 2 ) + 2 * ggvx[go] * gguy[go] + pow( ggvy[go], 2 )) 
                    + gght * ( ggux[go] + ggvy[go] ) + ( gu[go] * gghx + gv[go] * gghy ) * ( ggux[go] + ggvy[go] ));

                gw[go] = gu[go] * tan(gbetax) + gv[go] * tan(gbetay);
                gww[go] = (gu[go] * tan(gbetax) + gv[go] * tan( gbetay )) - ( ggux[go] + ggvy[go] ) * gh * 0.5;
                guvw[go] = pow( pow( gu[go], 2 ) + pow( gv[go], 2 ) + pow( gw[go], 2 ), 0.5 );
            }

            gdrag[0] = gcdrag[0] * pow( fabs( guvw[2] - guvw[0] ), sflow.DRAG_je - 1 );
            gdrag[1] = gcdrag[1] * pow( fabs( guvw[1] - guvw[0] ), sflow.DRAG_je - 1 );
            gdrag[2] = gcdrag[2] * pow( fabs( guvw[2] - guvw[1] ), sflow.DRAG_je - 1 );
            
            for ( gm=0; gm<2; gm++ ) {

                gb = 4 * gm;

                if ( galpha[0] > 0 ) ggze[gb+0] = ggze[gb+0] + ggtlength[0] - 1 / galpha[0] * ( gdrag[0] * ( gww[2] - gww[0] ) + gdrag[1] * ( gww[1] - gww[0] ));
                if ( galpha[1] > 0 ) ggze[gb+1] = ggze[gb+1] + ggtlength[1] - 1 / galpha[1] * ( -1 / ggamma2[1] * gdrag[1] * ( gww[1] - gww[0] ) + gdrag[2] * ( gww[2] - gww[1] ));
                if ( galpha[2] > 0 ) ggze[gb+2] = ggze[gb+2] + ggtlength[2] + 1 / galpha[2] * ( 1 / ggamma2[0] * gdrag[0] * ( gww[2] - gww[0] ) + 1 / ggamma2[2] * gdrag[2] * ( gww[2] - gww[1] ));
            }
            
            if ( galpha[1] > 0 ) ggze[8] = ggze[8] + ggtlength[1] - 1 / galpha[1] * ( gdrag[1] * ( gww[1] - gww[0] ) + gdrag[2] * ( gww[2] - gww[1] ));
            if ( galpha[1] > 0 ) ggze[10] = ggze[10] + ggtlength[1] - 1 / galpha[1] * ( gdrag[1] * ( gww[1] - gww[0] ) + gdrag[2] * ( gww[2] - gww[1] ));           
        }
        
    } else { // one-phase models (gravity only enhanced for curvature)

        ggze[0] = gtestx[0];
        ggze[3] = gtestx[0];
        ggze[4] = gtesty[0];
        ggze[7] = gtesty[0];
     
        ggze[8] = gtestx[0];
        ggze[9] = gtestx[0];
        ggze[10] = gtesty[0];
        ggze[11] = gtesty[0];

        ggze[1] = gtestx[0];
        ggze[5] = gtesty[0];

        ggze[2] = gtestx[0];
        ggze[6] = gtesty[0];
    }

    return ggze;
}

float *fdisp ( float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, float gvflowy2, float gvflowy3, 
    float *ggux, float *ggvy, float *gguy, float *ggvx, float gbetax, float gbetay, float *gkx, float *gky, float *gvm, float *gcdrag, float gtsum,
    struct ico sico, struct flow sflow ) { // function for dispersion

    int go, gp;
    float *gdisp, ggtlength[3], gh, galpha[3], guu[3]; // gu[3], gv[3], gw[3], guvw[3], ggamma[3], ggamma2[3], gcvm[3], gdrag[3], galphassfs, galphassff, galphafsff;
    gdisp = (float*) calloc( 6, sizeof(float));

    if ( sico.MODEL == 7 && ( sico.NOOSC == 3 || sico.NOOSC == 4 || sico.NOOSC == 6 || sico.NOOSC == 7 )) {

        gh = ghflow + ghflow2 + ghflow3;
        //gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
        //gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;
        
        if ( gh > sico.HFLOWMIN ) {

            if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
            if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
            if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

        /*if ( galpha[1] != 0 ) galphassfs = galpha[0] / galpha[1]; else galphassfs = 0;
        if ( galpha[2] != 0 ) galphassff = galpha[0] / galpha[2]; else galphassff = 0;
        if ( galpha[2] != 0 ) galphafsff = galpha[1] / galpha[2]; else galphafsff = 0;

        if ( galpha[0] > 0 && galpha[2] > 0 ) ggamma[0] = sflow.RHO3 / sflow.RHO1; else ggamma[0] = 0;
        if ( galpha[0] > 0 && galpha[1] > 0 ) ggamma[1] = sflow.RHO2 / sflow.RHO1; else ggamma[1] = 0;
        if ( galpha[1] > 0 && galpha[2] > 0 ) ggamma[2] = sflow.RHO3 / sflow.RHO2; else ggamma[2] = 0;

        if ( galpha[0] > 0 && galpha[2] > 0 ) ggamma2[0] = sflow.RHO3 / sflow.RHO1; else ggamma2[0] = 1;
        if ( galpha[0] > 0 && galpha[1] > 0 ) ggamma2[1] = sflow.RHO2 / sflow.RHO1; else ggamma2[1] = 1;
        if ( galpha[1] > 0 && galpha[2] > 0 ) ggamma2[2] = sflow.RHO3 / sflow.RHO2; else ggamma2[2] = 1;

        gcvm[0] = gvm[9]; gcvm[1] = gvm[10]; gcvm[2] = gvm[11];*/

        for ( go = 0; go < 3; go++ ) {

            //gw[go] = gu[go] * tan(gbetax) + gv[go] * tan(gbetay);
            //guvw[go] = pow( pow( gu[go], 2 ) + pow( gv[go], 2 ) + pow( gw[go], 2 ), 0.5 );
            guu[go] = ggux[go] + ggvy[go];
            ggtlength[go] = 2 * pow( guu[go], 2 ) - 2 * ggux[go] * ggvy[go] + 2 * ggvx[go] * gguy[go];
        }

        /*gdrag[0] = gcdrag[0] * pow( fabs( guvw[2] - guvw[0] ), sflow.DRAG_je - 1 );
        gdrag[1] = gcdrag[1] * pow( fabs( guvw[1] - guvw[0] ), sflow.DRAG_je - 1 );
        gdrag[2] = gcdrag[2] * pow( fabs( guvw[2] - guvw[1] ), sflow.DRAG_je - 1 );*/

        if ( galpha[0] > 0 )
            gdisp[0] = -gkx[0] * gh * gh / 3 * (( ggtlength[0] ));//- ggamma[0] * gcvm[0] * ( pow( guu[2], 2 ) - pow(guu[0], 2 )) - ggamma[1] * gcvm[1]
                //* ( pow( guu[1], 2 ) - pow(guu[0], 2 ))));
                //+ 1 / galpha[0] * gkx[0] * gh * gh / 6 * ( gdrag[0] * ( guu[2] - guu[0] ) + gdrag[1] * ( guu[1] - guu[0] ));
        else gdisp[0] = 0;

        if ( galpha[1] > 0 )
            gdisp[1] = -gh * gh / 3 * (( ggtlength[1] ));//- ggamma[2] * gcvm[2] * ( pow( guu[2], 2 ) - pow(guu[1], 2 )) + galphassfs * gcvm[1]
                //* ( pow( guu[1], 2 ) - pow(guu[0], 2 ))))
                //+ 1 / galpha[1] * gh * gh / 6 * ( - 1 / ggamma2[1] * gdrag[1] * ( guu[1] - guu[0] ) + gdrag[2] * ( guu[2] - guu[1] ));
        else gdisp[1] = 0;

        if ( galpha[2] > 0 )
            gdisp[2] = -gh * gh / 3 * (( ggtlength[2]));// + galphassff * gcvm[0] * ( pow( guu[2], 2 ) - pow(guu[0], 2 )) + galphafsff * gcvm[2]
                //* ( pow( guu[2], 2 ) - pow(guu[1], 2 ))))
                //- 1 / galpha[2] * gh * gh / 6 * ( 1 / ggamma2[0] * gdrag[0] * ( guu[2] - guu[0] ) + 1 / ggamma2[2] * gdrag[2] * ( guu[2] - guu[1] ));
        else gdisp[2] = 0;

        if ( galpha[0] > 0 )
            gdisp[3] = -gky[0] * gh * gh / 3 * (( ggtlength[0]));// - ggamma[0] * gcvm[0] * ( pow( guu[2], 2 ) - pow(guu[0], 2 )) - ggamma[1] * gcvm[1]
                //* ( pow( guu[1], 2 ) - pow(guu[0], 2 ))))
                //+ 1 / galpha[0] * gky[0] * gh * gh / 6 * ( gdrag[0] * ( guu[2] - guu[0] ) + gdrag[1] * ( guu[1] - guu[0] ));
        else gdisp[3] = 0;

        if ( galpha[1] > 0 )
            gdisp[4] = -gh * gh / 3 * (( ggtlength[1]));// - ggamma[2] * gcvm[2] * ( pow( guu[2], 2 ) - pow(guu[1], 2 )) + galphassfs * gcvm[1]
                //* ( pow( guu[1], 2 ) - pow(guu[0], 2 ))))
                //+ 1 / galpha[1] * gh * gh / 6 * ( - 1 / ggamma2[1] * gdrag[1] * ( guu[1] - guu[0] ) + gdrag[2] * ( guu[2] - guu[1] ));
        else gdisp[4] = 0;

        if ( galpha[2] > 0 )
            gdisp[5] = -gh * gh / 3 * (( ggtlength[2]));// + galphassff * gcvm[0] * ( pow( guu[2], 2 ) - pow(guu[0], 2 )) + galphafsff * gcvm[2]
                //* ( pow( guu[2], 2 ) - pow(guu[1], 2 ))))
                //- 1 / galpha[2] * gh * gh / 6 * ( 1 / ggamma2[0] * gdrag[0] * ( guu[2] - guu[0] ) + 1 / ggamma2[2] * gdrag[2] * ( guu[2] - guu[1] ));
        else gdisp[5] = 0;

        for ( gp=0; gp<6; gp++ ) gdisp[gp] *= sflow.disp_MULT;
        
    } else { gdisp[0] = 0; gdisp[1] = 0; gdisp[2] = 0; gdisp[3] = 0; gdisp[4] = 0; gdisp[5] = 0; }

    return gdisp;
}

float *ff ( float *ggh, float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy,
    float gvflowy2, float gvflowy3, float *gkx, float *ggz, float *gkappau, float *gvm, float gghflowx, float ggalphax, float *gcdrag, float *ggrav, float *gdisp, 
    float ggbetax, struct flow sflow, struct ico sico ) { // function for fluxes in x direction

    int gp;
    float *ggf, gh, gu[3], gv[3], galpha[3], gbetax[3], ggamma, glambdas, glambdaf, gseprate_x, gsepflux_xs, gsepflux_xf, ggtest;
    ggf = (float*) calloc( sico.NVECTMIN, sizeof(float));

    gh = ghflow + ghflow2 + ghflow3;
    gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
    gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

    if ( gh > sico.HFLOWMIN ) {

        if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
        if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
        if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

    } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }


    // Hydraulic pressure coefficients

    for ( gp = 0; gp < sico.PMAX; gp++ ) {

        if ( sico.PHASES[gp] == 1 && gp == 0 ) gbetax[gp] = gkx[gp] * ggz[0];
        else if ( sico.PHASES[gp] == 1 && gp == 1 ) gbetax[gp] = gkx[gp] * ggz[8];
        else gbetax[gp] = 0;
    }


    // Separation fluxes

    if ( sico.MODEL == 7 ) {

        if ( sico.SEPFLUX == 1 && ghflow > sico.HFLOWMIN && ghflow3 > sico.HFLOWMIN ) {

            if ( ggh[2] > sico.HFLOWMIN && ggh[5] > sico.HFLOWMIN ) ggalphax = ggalphax; else ggalphax = 0;

            ggamma = sflow.RHO3 / sflow.RHO1;

            glambdas = ggamma / ( galpha[0] + ggamma * ( 1 - galpha[0] ));
            glambdaf = 1 / ( galpha[0] + ggamma * ( 1 - galpha[0] ));
   
            if ( gcdrag[0] > 0 ) {

                gseprate_x = -1 / sflow.RHO1 / gcdrag[0] * ( 
                    -fsign( gu[0] ) * tan( sflow.DELTA[0] ) * ( 1 - ggamma ) * galpha[0] * sflow.RHO1 * ggrav[0]
                    + 0.5 * gkx[0] * ( 1 - ggamma ) * sflow.RHO1 * ggrav[0] * gh * ggalphax + gkx[0] * ( 1 - ggamma ) * sflow.RHO1 * ggrav[0] * galpha[0] * gghflowx
                    + ggamma * sflow.RHO1 * ggrav[0] * galpha[0] * gghflowx - sflow.RHO1 * ggrav[0] * galpha[0] * tan( ggbetax )
                    + sflow.sep_SHEAR * galpha[0] * sflow.NY[0] * sflow.RHO1 * gvflowx / ( ghflow * ghflow)); // separation rate in x direction
                
            } else gseprate_x = 0;

            gsepflux_xs = glambdas * gseprate_x * galpha[0] * ( 1 - galpha[0] ) * ( ghflow + ghflow3 ); // solid separation flux in x direction
            gsepflux_xf = glambdaf * gseprate_x * galpha[0] * ( 1 - galpha[0] ) * ( ghflow + ghflow3 ); // fluid separation flux in x direction

        } else { gsepflux_xs = 0; gsepflux_xf = 0; }

    } else { gsepflux_xs = 0; gsepflux_xf = 0; }


    // Flux terms

    if ( sico.MODEL == 0 ) {

        ggf[0] = gh * gu[0];
        ggf[1] = gh * gu[0] * gu[0] + gbetax[0] * gh * gh * 0.5;
        ggf[2] = gh * gu[0] * gv[0];

    } else {

        if ( galpha[0] > 0 ) {

            ggf[0] = galpha[0] * gh * gu[0] + gsepflux_xs;
            ggf[1] = galpha[0] * gh * ( gu[0] * gu[0] - gvm[1] + gbetax[0] * gh * 0.5 + gdisp[0] );

            if ( gdisp[0] != 0 ) {
            
                ggtest = galpha[0] * gh * ( gu[0] * gu[0] - gvm[1] + gbetax[0] * gh * 0.5 );
                if ( gdisp[0] != 0 && fsign( ggtest) != fsign( ggf[1] )) ggf[1] = 0;
                else if ( fabs( ggf[1] ) > sflow.CCONST * fabs( ggtest )) ggf[1] = sflow.CCONST * ggtest;
            }
            
            ggf[2] = galpha[0] * gh * ( gu[0] * gv[0] - gvm[2] );

        } else { ggf[0] = 0; ggf[1] = 0; ggf[2] = 0; }
    }

    if ( sico.MODEL == 7 ) {

        if ( galpha[1] > 0 ) {

            ggf[3] = galpha[1] * gh * gu[1];
            ggf[4] = galpha[1] * gh * ( gu[1] * gu[1] - gvm[4] + gbetax[1] * gh * 0.5 + gdisp[1] );
            
            if ( gdisp[1] != 0 ) {
            
                ggtest = galpha[1] * gh * ( gu[1] * gu[1] - gvm[4] + gbetax[1] * gh * 0.5 );
                if ( gdisp[1] != 0 && fsign( ggtest) != fsign( ggf[4] )) ggf[4] = 0;
                else if ( fabs( ggf[4] ) > sflow.CCONST * fabs( ggtest )) ggf[4] = sflow.CCONST * ggtest;
            }
            
            ggf[5] = galpha[1] * gh * ( gu[1] * gv[1] - gvm[5] );

        } else { ggf[3] = 0; ggf[4] = 0; ggf[5] = 0; }

        if ( galpha[2] > 0 ) {

            ggf[6] = galpha[2] * gh * gu[2] - gsepflux_xf;
            ggf[7] = galpha[2] * gh * ( gu[2] * gu[2] + gvm[7] + gbetax[2] * gh * 0.5 + gdisp[2] );
            
            if ( gdisp[2] != 0 ) {
            
                ggtest = galpha[2] * gh * ( gu[2] * gu[2] + gvm[7] + gbetax[2] * gh * 0.5 );
                if ( gdisp[2] != 0 && fsign( ggtest) != fsign( ggf[7] )) ggf[7] = 0;
                else if ( fabs( ggf[7] ) > sflow.CCONST * fabs( ggtest )) ggf[7] = sflow.CCONST * ggtest;
            }
          
            ggf[8] = galpha[2] * gh * ( gu[2] * gv[2] + gvm[8] ); 

        } else { ggf[6] = 0; ggf[7] = 0; ggf[8] = 0; }
    }

    return ggf;
}

float *fg ( float *ggh, float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, 
    float gvflowy2, float gvflowy3, float *gky, float *ggz, float *gkappau, float *gvm, float gghflowy, float ggalphay, float *gcdrag, float *ggrav, float *gdisp, 
    float ggbetay, struct flow sflow, struct ico sico ) { // function for fluxes in y direction

    int gp;
    float *ggg, gh, gu[3], gv[3], galpha[3], gbetay[3], ggamma, glambdas, glambdaf, gseprate_y, gsepflux_ys, gsepflux_yf, ggtest;
    ggg = (float*) calloc( sico.NVECTMIN, sizeof(float));

    gh = ghflow + ghflow2 + ghflow3;

    gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
    gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

    if ( gh > sico.HFLOWMIN ) {

        if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
        if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
        if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

    } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }


    // Hydraulic pressure coefficients

    for ( gp = 0; gp < sico.PMAX; gp++ ) {

        if ( sico.PHASES[gp] == 1 && gp == 0 ) gbetay[gp] = gky[gp] * ggz[4];
        else if ( sico.PHASES[gp] == 1 && gp == 1 ) gbetay[gp] = gky[gp] * ggz[10];
        else gbetay[gp] = 0;
    }
   

    // Separation fluxes

    if ( sico.MODEL == 7 ) {

        if ( sico.SEPFLUX == 1 && ghflow > sico.HFLOWMIN && ghflow3 > sico.HFLOWMIN ) {

            if ( ggh[1] > sico.HFLOWMIN && ggh[4] > sico.HFLOWMIN ) ggalphay = ggalphay; else ggalphay = 0;

            ggamma = sflow.RHO3 / sflow.RHO1;

            glambdas = ggamma / ( galpha[0] + ggamma * ( 1 - galpha[0] ));
            glambdaf = 1 / ( galpha[0] + ggamma * ( 1 - galpha[0] ));
      
            if ( gcdrag[0] > 0 ) {
         
                gseprate_y = -1 / sflow.RHO1 / gcdrag[0] * ( 
                    -fsign( gv[0] ) * tan( sflow.DELTA[0] ) * ( 1 - ggamma ) * galpha[0] * sflow.RHO1 * ggrav[1]
                    + 0.5 * gky[0] * ( 1 - ggamma ) * sflow.RHO1 * ggrav[1] * gh * ggalphay + gky[0] * ( 1 - ggamma ) * sflow.RHO1 * ggrav[1] * galpha[0] * gghflowy
                    + ggamma * sflow.RHO1 * ggrav[1] * galpha[0] * gghflowy - sflow.RHO1 * ggrav[1] * galpha[0] * tan( ggbetay )
                    + sflow.sep_SHEAR * galpha[0] * sflow.NY[0] * sflow.RHO1 * gvflowy / ( ghflow * ghflow )); // separation rate in y direction
                
            } else  gseprate_y = 0;

            gsepflux_ys = glambdas * gseprate_y * galpha[0] * ( 1 - galpha[0] ) * ( ghflow + ghflow3 ); // solid separation flux in y direction
            gsepflux_yf = glambdaf * gseprate_y * galpha[0] * ( 1 - galpha[0] ) * ( ghflow + ghflow3 ); // fluid separation flux in y direction

        } else { gsepflux_ys = 0; gsepflux_yf = 0; }

    } else { gsepflux_ys = 0; gsepflux_yf = 0; }


    // Flux terms

    if ( sico.MODEL == 0 ) {

        ggg[0] = gh * gv[0];
        ggg[1] = gh * gu[0] * gv[0];
        ggg[2] = gh * gv[0] * gv[0] + gbetay[0] * gh * gh * 0.5;

    } else {

        if ( galpha[0] > 0 ) {

            ggg[0] = galpha[0] * gh * gv[0] + gsepflux_ys;
            ggg[1] = galpha[0] * gh * ( gu[0] * gv[0] - gvm[14] );
            ggg[2] = galpha[0] * gh * ( gv[0] * gv[0] - gvm[13] + gbetay[0] * gh * 0.5 + gdisp[3] );
            
            if ( gdisp[3] != 0 ) {
            
                ggtest = galpha[0] * gh * ( gv[0] * gv[0] - gvm[13] + gbetay[0] * gh * 0.5 );
                if ( gdisp[3] != 0 && fsign( ggtest) != fsign( ggg[2] )) ggg[2] = 0;
                else if ( fabs( ggg[2] ) > sflow.CCONST * fabs( ggtest )) ggg[2] = sflow.CCONST * ggtest;
            }

        } else { ggg[0] = 0; ggg[1] = 0; ggg[2] = 0; }
    }

    if ( sico.MODEL == 7 ) {

        if ( galpha[1] > 0 ) {

            ggg[3] = galpha[1] * gh * gv[1];
            ggg[4] = galpha[1] * gh * ( gu[1] * gv[1] - gvm[17] ); 
            ggg[5] = galpha[1] * gh * ( gv[1] * gv[1] - gvm[16] + gbetay[1] * gh * 0.5 + gdisp[4] );
            
            if ( gdisp[4] != 0 ) {
            
                ggtest = galpha[1] * gh * ( gv[1] * gv[1] - gvm[16] + gbetay[1] * gh * 0.5 );
                if ( gdisp[4] != 0 && fsign( ggtest) != fsign( ggg[5] )) ggg[5] = 0;
                else if ( fabs( ggg[5] ) > sflow.CCONST * fabs( ggtest )) ggg[5] = sflow.CCONST * ggtest;
            }

        } else { ggg[3] = 0; ggg[4] = 0; ggg[5] = 0; }

        if ( galpha[2] > 0 ) {

            ggg[6] = galpha[2] * gh * gv[2] - gsepflux_yf;
            ggg[7] = galpha[2] * gh * ( gu[2] * gv[2] + gvm[20] );
            ggg[8] = galpha[2] * gh * ( gv[2] * gv[2] + gvm[19] + gbetay[2] * gh * 0.5 + gdisp[5] );

            if ( gdisp[5] != 0 ) {
            
                ggtest = galpha[2] * gh * ( gv[2] * gv[2] + gvm[19] + gbetay[2] * gh * 0.5 );
                if ( gdisp[5] != 0 && fsign( ggtest) != fsign( ggg[8] )) ggg[8] = 0;
                else if ( fabs( ggg[8] ) > sflow.CCONST * fabs( ggtest )) ggg[8] = sflow.CCONST * ggtest;
            }

        } else { ggg[6] = 0; ggg[7] = 0; ggg[8] = 0; }
    }

    return ggg;
}

float *fs ( float *gh, float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, float gvflowy2, float gvflowy3, 
    float gghx, float gghy, float ggalphassx, float ggalphassy, float ggalphafsx, float ggalphafsy, float ggalphaffx, float ggalphaffy, float *ggrav, float *ggz,
    float gdx, float gdy, float *gkappau, float *gcdrag, struct ico sico, struct flow sflow ) { // function for source terms (accelerating components)

    int gl, go;
    float *ggs, gu[3], gv[3], gw[3], galpha[3], ggalphax[3], ggalphay[3], ggamma0, ggamma1, ggamma2,
        gpbx[3], gpby[3], guvw[3], gdragx[3], gdragy[3], gcompx[3], gcompy[3];
    ggs = (float*) calloc( sico.NVECTMIN, sizeof(float));

    gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
    gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

    if ( sico.NOOSC == 4 || sico.NOOSC == 7 ) { // non-hydrostatic effects

        gw[0] = gu[0] * tan(gdx) + gv[0] * tan(gdy);
        gw[1] = gu[1] * tan(gdx) + gv[1] * tan(gdy);
        gw[2] = gu[2] * tan(gdx) + gv[2] * tan(gdy);

    } else { gw[0] = 0; gw[1] = 0; gw[2] = 0; }

    if ( gh[0] > sico.HFLOWMIN ) {

        if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh[0]; else galpha[0] = 0;
        if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh[0]; else galpha[1] = 0;
        if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh[0]; else galpha[2] = 0;

    } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

    if ( galpha[0] > 0  && galpha[2] > 0 ) ggamma0 = sflow.RHO3 / sflow.RHO1; else ggamma0 = 1;
    if ( galpha[0] > 0  && galpha[1] > 0 ) ggamma1 = sflow.RHO2 / sflow.RHO1; else ggamma1 = 1;
    if ( galpha[1] > 0  && galpha[2] > 0 ) ggamma2 = sflow.RHO3 / sflow.RHO2; else ggamma2 = 1;

    for ( go = 0; go < 3; go++ ) guvw[go] = pow( pow( gu[go], 2 ) + pow( gv[go], 2 ) + pow( gw[go], 2 ), 0.5 );

    if ( gh[2] > sico.HFLOWMIN && gh[5] > sico.HFLOWMIN ) ggalphax[0] = ggalphassx; else ggalphax[0] = 0;
    if ( gh[2] > sico.HFLOWMIN && gh[5] > sico.HFLOWMIN ) ggalphax[1] = ggalphafsx; else ggalphax[1] = 0;
    if ( gh[2] > sico.HFLOWMIN && gh[5] > sico.HFLOWMIN ) ggalphax[2] = ggalphaffx; else ggalphax[2] = 0;
    if ( gh[1] > sico.HFLOWMIN && gh[4] > sico.HFLOWMIN ) ggalphay[0] = ggalphassy; else ggalphay[0] = 0;
    if ( gh[1] > sico.HFLOWMIN && gh[4] > sico.HFLOWMIN ) ggalphay[1] = ggalphafsy; else ggalphay[1] = 0;
    if ( gh[1] > sico.HFLOWMIN && gh[4] > sico.HFLOWMIN ) ggalphay[2] = ggalphaffy; else ggalphay[2] = 0;


    // Effective basal pressures and phase-dependent components

    if ( sico.MODEL > 0 ) {

        for ( gl = 0; gl < sico.PMAX; gl++ ) {

            if ( gl == 0 && sico.PHASES[gl] == 1 ) {

                gcompx[gl] = galpha[gl] * ( ggrav[2] - ggz[3] * gghx );
                gcompy[gl] = galpha[gl] * ( ggrav[3] - ggz[7] * gghy );

            } else if ( gl == 1 && sico.PHASES[gl] == 1 ) {

                gcompx[gl] = galpha[gl] * ( ggrav[2] - ggz[9] * gghx );
                gcompy[gl] = galpha[gl] * ( ggrav[3] - ggz[11] * gghy );

            } else if ( sico.PHASES[gl] == 2 ) {

                gpbx[gl] = ggz[gl];
                gpby[gl] = ggz[gl+4];

                gcompx[gl] = galpha[gl] * ( ggrav[2] - ggz[gl] * gghx - ( -0.5 * gpbx[gl] * gh[0] / galpha[gl] * ggalphax[gl] ));
                gcompy[gl] = galpha[gl] * ( ggrav[3] - ggz[gl+4] * gghy - ( -0.5 * gpby[gl] * gh[0] / galpha[gl] * ggalphay[gl] ));

            } else if ( sico.PHASES[gl] == 3 ) {

                gpbx[gl] = ggz[gl];
                gpby[gl] = ggz[gl+4];

                gcompx[gl] = galpha[gl] * ( ggrav[2] - ggz[gl] * gghx - ( -0.5 * gpbx[gl] * gh[0] / galpha[gl] * ggalphax[gl] ));
                gcompy[gl] = galpha[gl] * ( ggrav[3] - ggz[gl+4] * gghy - ( -0.5 * gpby[gl] * gh[0] / galpha[gl] * ggalphay[gl] ));
            }
        }
    }


    // Drag terms

    if ( sico.MODEL == 7 ) {

        gdragx[0] = gcdrag[0] * ( gu[2] - gu[0] ) * pow( fabs( guvw[2] - guvw[0] ), sflow.DRAG_je - 1 );
        gdragx[1] = gcdrag[1] * ( gu[1] - gu[0] ) * pow( fabs( guvw[1] - guvw[0] ), sflow.DRAG_je - 1 );
        gdragx[2] = gcdrag[2] * ( gu[2] - gu[1] ) * pow( fabs( guvw[2] - guvw[1] ), sflow.DRAG_je - 1 );

        gdragy[0] = gcdrag[0] * ( gv[2] - gv[0] ) * pow( fabs( guvw[2] - guvw[0] ), sflow.DRAG_je - 1 );
        gdragy[1] = gcdrag[1] * ( gv[1] - gv[0] ) * pow( fabs( guvw[1] - guvw[0] ), sflow.DRAG_je - 1 );
        gdragy[2] = gcdrag[2] * ( gv[2] - gv[1] ) * pow( fabs( guvw[2] - guvw[1] ), sflow.DRAG_je - 1 );

    } else { gdragx[0] = 0; gdragx[1] = 0; gdragx[2] = 0; gdragy[0] = 0; gdragy[1] = 0; gdragy[2] = 0; }


    // Source terms

    if ( sico.MODEL == 0 ) {

        ggs[0] = 0; 
        ggs[1] = ggrav[2] * gh[0];
        ggs[2] = ggrav[3] * gh[0];

    } else {

        if ( galpha[0] > 0 ) {
 
            ggs[0] = 0;
            ggs[1] = gh[0] * ( gcompx[0] + gdragx[0] + gdragx[1] );
            ggs[2] = gh[0] * ( gcompy[0] + gdragy[0] + gdragy[1] );

        } else { ggs[0] = 0; ggs[1] = 0; ggs[2] = 0; }
    }

    if ( sico.MODEL == 7 ) {

        if ( galpha[1] > 0 ) {

            ggs[3] = 0; 
            ggs[4] = gh[0] * ( gcompx[1] - 1 / ggamma1 * gdragx[1] + gdragx[2] );
            ggs[5] = gh[0] * ( gcompy[1] - 1 / ggamma1 * gdragy[1] + gdragy[2] );

        } else { ggs[3] = 0; ggs[4] = 0; ggs[5] = 0; }

        if ( galpha[2] > 0 ) {

            ggs[6] = 0; 
            ggs[7] = gh[0] * ( gcompx[2] - 1 / ggamma0 * gdragx[0] - 1 / ggamma2 * gdragx[2] );
            ggs[8] = gh[0] * ( gcompy[2] - 1 / ggamma0 * gdragy[0] - 1 / ggamma2 * gdragy[2] );

        } else { ggs[6] = 0; ggs[7] = 0; ggs[8] = 0; }
    }

    return ggs;
}

float *fd ( float *gh, float gghx, float gghy, float *gnuss, float *gnvss, float *gnufs, float *gnvfs,
    float *gnuff, float *gnvff, float *ggalpha10, float *ggalpha20, float *ggalpha30, float *gnbetax, float *gnbetay, float *ggz, float gdx, float gdy,
    float *gd0, float *gkappau, float gekin, float *ggrav, struct ico sico, struct flow sflow ) { // function for source terms (decelerating components)

    int gg, gj, gl, go, gq;
    float *ggd, gd[3][6], gw[3][6], guvw[3][6], ggux[3][6], gguy[3][6], ggvx[3][6], ggvy[3][6], ggalpha[3][9], ggalphax[3][6], 
        ggalphay[3][6], gfricx[3], gfricy[3], gvisx[3], gvisy[3], gcambdragx, gcambdragy, gp, gtauy[3], gnye[3][6], gcuf[3][6], gguzb[3][6], gcvf[3][6], 
        ggvzb[3][6], gtaunnx[3], gtaunny[3], gflufri[3], gdelta, gphi, gny[6], difuxx[6], difvxx[6], difuyx[6], difvyy[6], difuyy[6], difvxy[6], difax, difay, gpbx[3], gpby[3], 
        guratio[3], gvratio[3], gnu[3][9], gnv[3][9], gxterm1, gyterm1, gxterm2, gyterm2, gxterm3, gyterm3;

    ggd = (float*) calloc( sico.NVECTMIN+18, sizeof(float));

    for( gj=0; gj<9; gj++ ) { 

        gnu[0][gj] = gnuss[gj]; gnu[1][gj] = gnufs[gj]; gnu[2][gj] = gnuff[gj]; gnv[0][gj] = gnvss[gj]; gnv[1][gj] = gnvfs[gj]; gnv[2][gj] = gnvff[gj];
        ggalpha[0][gj] = ggalpha10[gj]; ggalpha[1][gj] = ggalpha20[gj]; ggalpha[2][gj] = ggalpha30[gj];
    }

    for ( gq=0; gq<6; gq++ ) {

        gd[0][gq] = gd0[3*gq]; gd[1][gq] = gd0[3*gq+1]; gd[2][gq] = gd0[3*gq+2];

        for ( go = 0; go < 3; go++ ) {

            if ( sico.NOOSC == 4 || sico.NOOSC == 7 ) gw[0][gq] = gnu[go][gq] * tan( gnbetax[gq] ) + gnv[go][gq] * tan( gnbetay[gq] ); // non-hydrostatic effects
            else gw[go][gq] = 0;

            guvw[go][gq] = pow( pow( gnu[go][gq], 2 ) + pow( gnv[go][gq], 2 ) + pow( gw[go][gq], 2 ), 0.5 );
            if ( guvw[go][0] > 0 ) guratio[go] = gnu[go][0] / guvw[go][0]; else guratio[go] = 0;
            if ( guvw[go][0] > 0 ) gvratio[go] = gnv[go][0] / guvw[go][0]; else gvratio[go] = 0;
        }
    }

    if ( sico.MODEL > 0 ) {


        // Gradients of velocities and fractions

        for ( go=0; go<3; go++ ) {

            ggux[go][0] = ( gnu[go][2] - gnu[go][5] ) / ( 2 * gdx );
            gguy[go][0] = ( gnu[go][1] - gnu[go][4] ) / ( 2 * gdy );
            ggvx[go][0] = ( gnv[go][2] - gnv[go][5] ) / ( 2 * gdx );
            ggvy[go][0] = ( gnv[go][1] - gnv[go][4] ) / ( 2 * gdy );

            ggux[go][1] = ( gnu[go][3] - gnu[go][8] ) / ( 2 * gdx );
            gguy[go][1] = ( gnu[go][1] - gnu[go][0] ) / ( 1 * gdy );
            ggvx[go][1] = ( gnv[go][3] - gnv[go][8] ) / ( 2 * gdx );
            ggvy[go][1] = ( gnv[go][1] - gnv[go][0] ) / ( 1 * gdy );

            ggux[go][2] = ( gnu[go][2] - gnu[go][0] ) / ( 1 * gdx );
            gguy[go][2] = ( gnu[go][3] - gnu[go][7] ) / ( 2 * gdy );
            ggvx[go][2] = ( gnv[go][2] - gnv[go][0] ) / ( 1 * gdx );
            ggvy[go][2] = ( gnv[go][3] - gnv[go][7] ) / ( 2 * gdy );

            ggux[go][4] = ( gnu[go][7] - gnu[go][6] ) / ( 2 * gdx );
            gguy[go][4] = ( gnu[go][0] - gnu[go][4] ) / ( 1 * gdy );
            ggvx[go][4] = ( gnv[go][7] - gnv[go][6] ) / ( 2 * gdx );
            ggvy[go][4] = ( gnv[go][0] - gnv[go][4] ) / ( 1 * gdy );

            ggux[go][5] = ( gnu[go][0] - gnu[go][5] ) / ( 1 * gdx );
            gguy[go][5] = ( gnu[go][8] - gnu[go][6] ) / ( 2 * gdy );
            ggvx[go][5] = ( gnv[go][0] - gnv[go][5] ) / ( 1 * gdx );
            ggvy[go][5] = ( gnv[go][8] - gnv[go][6] ) / ( 2 * gdy );

            ggalphax[go][0] = ( ggalpha[go][2] - ggalpha[go][5] ) / ( 2 * gdx );
            ggalphay[go][0] = ( ggalpha[go][1] - ggalpha[go][4] ) / ( 2 * gdy );

            ggalphax[go][1] = ( ggalpha[go][3] - ggalpha[go][8] ) / ( 2 * gdx );
            ggalphay[go][1] = ( ggalpha[go][1] - ggalpha[go][0] ) / ( 1 * gdy );

            ggalphax[go][2] = ( ggalpha[go][2] - ggalpha[go][0] ) / ( 1 * gdx );
            ggalphay[go][2] = ( ggalpha[go][3] - ggalpha[go][7] ) / ( 2 * gdy );

            ggalphax[go][4] = ( ggalpha[go][7] - ggalpha[go][6] ) / ( 2 * gdx );
            ggalphay[go][4] = ( ggalpha[go][0] - ggalpha[go][4] ) / ( 1 * gdy );

            ggalphax[go][5] = ( ggalpha[go][0] - ggalpha[go][5] ) / ( 1 * gdx );
            ggalphay[go][5] = ( ggalpha[go][8] - ggalpha[go][6] ) / ( 2 * gdy );
        }


        // Effective ambient drag coefficients

        if ( gh[0] > sico.HFLOWMIN ) {

            if ( fsign( gnu[0][0] + gnu[1][0] + gnu[2][0] ) == fsign( gghx )) gcambdragx = 0;
            else gcambdragx = gdx * fabs( gghx ) / gh[0] * sflow.AMBDRAG;
            gcambdragx = ffmin( sflow.AMBDRAG, gcambdragx );

            if ( fsign( gnv[0][0] + gnv[1][0] + gnv[2][0] ) == fsign( gghy )) gcambdragy = 0;
            else gcambdragy = gdy * fabs( gghy ) / gh[0] * sflow.AMBDRAG;
            gcambdragy = ffmin( sflow.AMBDRAG, gcambdragy );

        } else { gcambdragx = 0; gcambdragy = 0; }

        for ( gl = 0; gl < sico.PMAX; gl++ ) {


            // Solid phase

            if ( sico.PHASES[gl] == 1 ) {

                if ( sico.DYNFRIC == 1 ) gdelta = fdynfric( sflow.DELTA[gl], sico.FRICMIN, gekin, ggalpha[gl][0], sico ); else gdelta = sflow.DELTA[gl];

                if ( gl == 0 ) {

                    gpbx[gl] = ggz[0];
                    gpby[gl] = ggz[4];

                } else if ( gl == 1 ) {

                    gpbx[gl] = ggz[8];
                    gpby[gl] = ggz[10]; // effective basal pressures
                }

                if ( sico.CURVCTRL == 1 || sico.CURVCTRL == 4 ) { gpbx[gl] += gkappau[gl]; gpby[gl] += gkappau[gl]; } // curvature effects

                gfricx[gl] = guratio[gl] * tan( gdelta ) * gpbx[gl];
                gfricy[gl] = gvratio[gl] * tan( gdelta ) * gpby[gl];
                gflufri[gl] = 0; // friction terms
                
                gvisx[gl] = 0;
                gvisy[gl] = 0; // viscosity terms


            // Fine-solid phase

            } else if ( sico.PHASES[gl] == 2 ) {

                gfricx[gl] = 0;
                gfricy[gl] = 0;
                gflufri[gl] = 0; // friction terms
                
                if ( sico.DYNFRIC == 1 ) gdelta = fdynfric( sflow.DELTA[gl], sico.FRICMIN, gekin, ggalpha[gl][0], sico ); else gdelta = sflow.DELTA[gl];
                if ( sico.DYNFRIC == 1 ) gphi = fdynfric( sflow.PHI[gl], sico.FRICMIN, gekin, ggalpha[gl][0], sico ); else gphi = sflow.PHI[gl];
                if ( gphi == sico.UNDEF ) gphi = gdelta;
            
                for ( gq=0; gq<6; gq++ ) {

                    if ( gl == 1 || sico.PMAX < 3 ) gp = ggrav[4];
                    else gp = ggrav[4] * sflow.RHO3 / sflow.RHO2;
                    
                    if ( sico.CURVCTRL == 1 || sico.CURVCTRL == 4 ) gp += gkappau[gl]; // curvature effects
                    
                    if ( sflow.TAUY[gl] == -9999 ) gtauy[gl] = sin( gphi ) * gp;
                    else gtauy[gl] = sflow.TAUY[gl]; // yield strength

                    gny[gq] = 1 / sico.SLOMO * sflow.NY[gl] * pow( ggalpha[gl][gq], sflow.NY_exp );
                    if ( gd[gl][gq] != 0 ) gnye[gl][gq] = gny[gq] + gtauy[gl] / gd[gl][gq] * ( 1 - exp( -sflow.VIS_ry * gd[gl][gq] )); else gnye[gl][gq] = gny[gq]; // effective kinematic viscosity

                    if ( gphi > 0 && gdelta > 0 ) {

                        if ( guvw[gl][gq] != 0 ) gcuf[gl][gq] = -gnu[gl][gq] / fabs( guvw[gl][gq] ) * tan( gdelta ); else gcuf[gl][gq] = 0;
                        if ( gnye[gl][gq] != 0 ) gguzb[gl][gq] = gcuf[gl][gq] / gnye[gl][gq] * gp + 2 * gcuf[gl][gq] * ggux[gl][gq]; else gguzb[gl][gq] = 0;

                        if ( guvw[gl][gq] != 0 ) gcvf[gl][gq] = -gnv[gl][gq] / fabs( guvw[gl][gq] ) * tan( gdelta ); else gcvf[gl][gq] = 0;
                        if ( gnye[gl][gq] != 0 ) ggvzb[gl][gq] = gcvf[gl][gq] / gnye[gl][gq] * gp + 2 * gcvf[gl][gq] * ggvy[gl][gq]; else ggvzb[gl][gq] = 0;
                            // fluid-type basal shear stresses (slip condition)

                    } else { 

                        if ( gh[gq] > sico.HFLOWMIN ) gguzb[gl][gq] = sflow.VIS_chi[gl] * gnu[gl][gq] / gh[gq]; else gguzb[gl][gq] = 0;
                        if ( gh[gq] > sico.HFLOWMIN ) ggvzb[gl][gq] = sflow.VIS_chi[gl] * gnv[gl][gq] / gh[gq]; else ggvzb[gl][gq] = 0; // fluid-type basal shear stresses (no-slip condition)
                    }

                    if ( gh[gq] > sico.HFLOWMIN )

                        gd[gl][gq] = 1 / sico.SLOMO * 
                            pow( fabs( 4 * ggux[gl][gq] * ggvy[gl][gq] - pow( gguy[gl][gq] + ggvx[gl][gq], 2 ) - pow( gguzb[gl][gq], 2 ) - pow( ggvzb[gl][gq], 2 )), 0.5 );
                            // second invariant of the strain-rate tensor

                    else gd[gl][gq] = 0;
                }

                if ( ggalpha[gl][0] > 0 ) {

                    gtaunnx[gl] = 0; gtaunny[gl] = 0;

                    for ( gg = 0; gg < sico.PMAX; gg++ ) {

                        if ( gg != gl && ( sico.PHASES[gg] == 1 )) {

                            difuxx[5] = ggalphax[gg][5] * ( gnu[gl][5] - gnu[gg][5] );
                            difuxx[2] = ggalphax[gg][2] * ( gnu[gl][2] - gnu[gg][2] );

                            difvxx[4] = ggalphax[gg][4] * ( gnv[gl][4] - gnv[gg][4] );
                            difvxx[1] = ggalphax[gg][1] * ( gnv[gl][1] - gnv[gg][1] );

                            difuyx[4] = ggalphay[gg][4] * ( gnu[gl][4] - gnu[gg][4] );
                            difuyx[1] = ggalphay[gg][1] * ( gnu[gl][1] - gnu[gg][1] );

                            difvyy[4] = ggalphay[gg][4] * ( gnv[gl][4] - gnv[gg][4] );
                            difvyy[1] = ggalphay[gg][1] * ( gnv[gl][1] - gnv[gg][1] );

                            difuyy[5] = ggalphay[gg][5] * ( gnu[gl][5] - gnu[gg][5] );
                            difuyy[2] = ggalphay[gg][2] * ( gnu[gl][2] - gnu[gg][2] );

                            difvxy[5] = ggalphax[gg][5] * ( gnv[gl][5] - gnv[gg][5] );
                            difvxy[2] = ggalphax[gg][2] * ( gnv[gl][2] - gnv[gg][2] );

                            difax = gnu[gl][0] - gnu[gg][0];
                            difay = gnv[gl][0] - gnv[gg][0];

                            gxterm1 = gnye[gl][2] * difuxx[2] - gnye[gl][5] * difuxx[5];
                            gxterm2 = gnye[gl][1] * ( difvxx[1] + difuyx[1] ) - gnye[gl][4] * ( difvxx[4] - difuyx[4] );
                            gyterm1 = gnye[gl][1] * difvyy[1] - gnye[gl][4] * difvyy[4];
                            gyterm2 = gnye[gl][2] * ( difuyy[2] + difvxy[2] ) - gnye[gl][5] * ( difuyy[5] - difvxy[5] );

                            gtaunnx[gl] += pow( 1 - ggalpha[gg][0], sflow.VIS_a[gg] ) / ggalpha[gl][0] * ( 2 / gdx * gxterm1 + 1 / gdy * gxterm2 ) 
                                - pow( 1 - ggalpha[gg][0], sflow.VIS_a[gg] ) / ggalpha[gl][0] * sflow.VIS_xi[gg] * ggalpha[gg][0] * gnye[gl][0] * difax / pow( gh[0], 2);

                            gtaunny[gl] += pow( 1 - ggalpha[gg][0], sflow.VIS_a[gg] ) / ggalpha[gl][0] * ( 2 / gdy * gyterm1 + 1 / gdx * gyterm2 ) 
                                - pow( 1 - ggalpha[gg][0], sflow.VIS_a[gg] ) / ggalpha[gl][0] * sflow.VIS_xi[gg] * ggalpha[gg][0] * gnye[gl][0] * difay / pow( gh[0], 2);
                                // enhanced non-Newtonan viscous stresses
                        }
                    }

                } else { gtaunnx[gl] = 0; gtaunny[gl] = 0; }

                if ( ggalpha[gl][0] > 0 && guvw[gl][0] != 0 ) {

                    gxterm1 = 2 / gdx * ( gnye[gl][2] * ggux[gl][2] - gnye[gl][5] * ggux[gl][5] );
                    gxterm2 = 1 / gdy * ( gnye[gl][1] * ggvx[gl][1] - gnye[gl][4] * ggvx[gl][4] );
                    gxterm3 = 1 / gdy * ( gnye[gl][1] * gguy[gl][1] - gnye[gl][4] * gguy[gl][4] );

                    gyterm1 = 2 / gdy * ( gnye[gl][1] * ggvy[gl][1] - gnye[gl][4] * ggvy[gl][4] );
                    gyterm2 = 1 / gdx * ( gnye[gl][2] * gguy[gl][2] - gnye[gl][5] * gguy[gl][5] );
                    gyterm3 = 1 / gdx * ( gnye[gl][2] * ggvx[gl][2] - gnye[gl][5] * ggvx[gl][5] );

                    gvisx[gl] = -( gxterm1 + gxterm2 + gxterm3 - gnye[gl][0] * gguzb[gl][0] / gh[0] ) + gtaunnx[gl];
                    gvisy[gl] = -( gyterm1 + gyterm2 + gyterm3 - gnye[gl][0] * ggvzb[gl][0] / gh[0] ) + gtaunny[gl]; // viscosity terms

                } else { gvisx[gl] = 0; gvisy[gl] = 0; }


            // Fluid phase

            } else if ( sico.PHASES[gl] == 3 ) {

                gfricx[gl] = 0;
                gfricy[gl] = 0;
                if ( sico.DYNFRIC == 1 ) gflufri[gl] = fdynfric( sflow.FLUFRI, 0.0, gekin, ggalpha[gl][0], sico ); else gflufri[gl] = sflow.FLUFRI; // friction terms

                if ( sflow.TAUY[gl] == -9999 ) gtauy[gl] = 0;
                else gtauy[gl] = sflow.TAUY[gl]; // yield strength

                for ( gq=0; gq<6; gq++ ) {

                    gny[gq] = 1 / sico.SLOMO * sflow.NY[gl] * pow( ggalpha[gl][gq], sflow.NY_exp );
                    if ( gd[gl][gq] != 0 ) gnye[gl][gq] = gny[gq] + gtauy[gl] / gd[gl][gq] * ( 1 - exp( -sflow.VIS_ry * gd[gl][gq] )); else gnye[gl][gq] = gny[gq]; // effective kinematic viscosity

                    if ( gh[gq] > sico.HFLOWMIN ) gguzb[gl][gq] = sflow.VIS_chi[gl] * gnu[gl][gq] / gh[gq]; else gguzb[gl][gq] = 0;
                    if ( gh[gq] > sico.HFLOWMIN ) ggvzb[gl][gq] = sflow.VIS_chi[gl] * gnv[gl][gq] / gh[gq]; else ggvzb[gl][gq] = 0; // fluid-type basal shear stresses (no-slip conditions)
                
                    if ( gh[gq] > sico.HFLOWMIN )

                        gd[gl][gq] = 1 / sico.SLOMO *  
                            pow( fabs( 4 * ggux[gl][gq] * ggvy[gl][gq] - pow( gguy[gl][gq] + ggvx[gl][gq], 2 ) - pow( gguzb[gl][gq], 2 ) - pow( ggvzb[gl][gq], 2 )), 0.5 );
                            // second invariant of the strain-rate tensor

                    else gd[gl][gq] = 0;
                }

                if ( ggalpha[gl][0] > 0 ) {

                    gtaunnx[gl] = 0; gtaunny[gl] = 0;

                    for ( gg = 0; gg < sico.PMAX; gg++ ) {

                        if ( gg != gl && ( sico.PHASES[gg] == 1 || sico.PHASES[gg] == 2 )) {

                            difuxx[5] = ggalphax[gg][5] * ( gnu[gl][5] - gnu[gg][5] );
                            difuxx[2] = ggalphax[gg][2] * ( gnu[gl][2] - gnu[gg][2] );

                            difvxx[4] = ggalphax[gg][4] * ( gnv[gl][4] - gnv[gg][4] );
                            difvxx[1] = ggalphax[gg][1] * ( gnv[gl][1] - gnv[gg][1] );

                            difuyx[4] = ggalphay[gg][4] * ( gnu[gl][4] - gnu[gg][4] );
                            difuyx[1] = ggalphay[gg][1] * ( gnu[gl][1] - gnu[gg][1] );

                            difvyy[4] = ggalphay[gg][4] * ( gnv[gl][4] - gnv[gg][4] );
                            difvyy[1] = ggalphay[gg][1] * ( gnv[gl][1] - gnv[gg][1] );

                            difuyy[5] = ggalphay[gg][5] * ( gnu[gl][5] - gnu[gg][5] );
                            difuyy[2] = ggalphay[gg][2] * ( gnu[gl][2] - gnu[gg][2] );

                            difvxy[5] = ggalphax[gg][5] * ( gnv[gl][5] - gnv[gg][5] );
                            difvxy[2] = ggalphax[gg][2] * ( gnv[gl][2] - gnv[gg][2] );

                            difax = gnu[gl][0] - gnu[gg][0];
                            difay = gnv[gl][0] - gnv[gg][0];

                            gxterm1 = gnye[gl][2] * difuxx[2] - gnye[gl][5] * difuxx[5];
                            gxterm2 = gnye[gl][1] * ( difvxx[1] + difuyx[1] ) - gnye[gl][4] * ( difvxx[4] - difuyx[4] );
                            gyterm1 = gnye[gl][1] * difvyy[1] - gnye[gl][4] * difvyy[4];
                            gyterm2 = gnye[gl][2] * ( difuyy[2] + difvxy[2] ) - gnye[gl][5] * ( difuyy[5] - difvxy[5] );

                            gtaunnx[gl] += pow( 1 - ggalpha[gg][0], sflow.VIS_a[sico.PHASES[gg]] ) / ggalpha[gl][0] * ( 2 / gdx * gxterm1 + 1 / gdy * gxterm2 ) 
                                - pow( 1 - ggalpha[gg][0], sflow.VIS_a[sico.PHASES[gg]] ) / ggalpha[gl][0] * sflow.VIS_xi[gg] * ggalpha[gg][0] * gnye[gl][0] * difax / pow( gh[0], 2);

                            gtaunny[gl] += pow( 1 - ggalpha[gg][0], sflow.VIS_a[sico.PHASES[gg]] ) / ggalpha[gl][0] * ( 2 / gdy * gyterm1 + 1 / gdx * gyterm2 ) 
                                - pow( 1 - ggalpha[gg][0], sflow.VIS_a[sico.PHASES[gg]] ) / ggalpha[gl][0] * sflow.VIS_xi[gg] * ggalpha[gg][0] * gnye[gl][0] * difay / pow( gh[0], 2);
                                // enhanced non-Newtonan viscous stresses
                        }
                    }

                } else { gtaunnx[gl] = 0; gtaunny[gl] = 0; }

                if ( ggalpha[gl][0] > 0 ) {

                    gxterm1 = 2 / gdx * ( gnye[gl][2] * ggux[gl][2] - gnye[gl][5] * ggux[gl][5] );
                    gxterm2 = 1 / gdy * ( gnye[gl][1] * ggvx[gl][1] - gnye[gl][4] * ggvx[gl][4] );
                    gxterm3 = 1 / gdy * ( gnye[gl][1] * gguy[gl][1] - gnye[gl][4] * gguy[gl][4] );

                    gyterm1 = 2 / gdy * ( gnye[gl][1] * ggvy[gl][1] - gnye[gl][4] * ggvy[gl][4] );
                    gyterm2 = 1 / gdx * ( gnye[gl][2] * gguy[gl][2] - gnye[gl][5] * gguy[gl][5] );
                    gyterm3 = 1 / gdx * ( gnye[gl][2] * ggvx[gl][2] - gnye[gl][5] * ggvx[gl][5] );

                    gvisx[gl] = -( gxterm1 + gxterm2 + gxterm3 - gnye[gl][0] * gguzb[gl][0] / gh[0] ) + gtaunnx[gl];
                    gvisy[gl] = -( gyterm1 + gyterm2 + gyterm3 - gnye[gl][0] * ggvzb[gl][0] / gh[0] ) + gtaunny[gl]; // viscosity terms

                } else { gvisx[gl] = 0; gvisy[gl] = 0; }
            }
        }


        // Deceleration terms

        ggd[0] = 0;

        if ( gnu[0][0] != 0 )
            ggd[1] = gh[0] * ggalpha[0][0] * (( gcambdragx + pow( gflufri[0], 2 ) / pow( ffmax(sico.HFLOWMIN, gh[0] ), 1.3333 )) 
                * gnu[0][0] * guvw[0][0] + gfricx[0] + gvisx[0] );
        else ggd[1] = 0;

        if ( gnv[0][0] != 0 )
            ggd[2] = gh[0] * ggalpha[0][0] * (( gcambdragy + pow( gflufri[0], 2 ) / pow( ffmax(sico.HFLOWMIN, gh[0] ), 1.3333 )) 
                * gnv[0][0] * guvw[0][0] + gfricy[0] + gvisy[0] );
        else ggd[2] = 0;

        if ( sico.MODEL == 7 ) {

            ggd[3] = 0; 

            if ( gnu[1][0] != 0 )
                ggd[4] = gh[0] * ggalpha[1][0] * (( gcambdragx + pow( gflufri[1], 2 ) / pow( ffmax(sico.HFLOWMIN, gh[0] ), 1.3333 ))
                    * gnu[1][0] * guvw[1][0] + gfricx[1] + gvisx[1] );
            else ggd[4] = 0;

            if ( gnv[1][0] != 0 )
                ggd[5] = gh[0] * ggalpha[1][0] * (( gcambdragy + pow( gflufri[1], 2 ) / pow( ffmax(sico.HFLOWMIN, gh[0] ), 1.3333 ))
                    * gnv[1][0] * guvw[1][0] + gfricy[1] + gvisy[1] );
            else ggd[5] = 0;

            ggd[6] = 0;

            if ( gnu[2][0] != 0 )
                ggd[7] = gh[0] * ggalpha[2][0] * (( gcambdragx + pow( gflufri[2], 2 ) / pow( ffmax(sico.HFLOWMIN, gh[0] ), 1.3333 ))
                    * gnu[2][0] * guvw[2][0] + gfricx[2] + gvisx[2] );
            else ggd[7] = 0;

            if ( gnv[2][0] != 0 )
                ggd[8] = gh[0] * ggalpha[2][0] * (( gcambdragy + pow( gflufri[2], 2 ) / pow( ffmax(sico.HFLOWMIN, gh[0] ), 1.3333 ))
                    * gnv[2][0] * guvw[2][0] + gfricy[2] + gvisy[2] );
            else ggd[8] = 0;
        }

        for ( gq=0; gq<6; gq++ ) {

            ggd[sico.NVECTMIN+3*gq] = gd[0][gq];
            ggd[sico.NVECTMIN+3*gq+1] = gd[1][gq];
            ggd[sico.NVECTMIN+3*gq+2] = gd[2][gq];
        }

    } else {

        gdelta = sflow.DELTA[0];

        gpbx[0] = ggz[0];
        gpby[0] = ggz[4];
        
        if ( sico.CURVCTRL == 1 || sico.CURVCTRL == 4 ) { gpbx[0] += gkappau[0]; gpby[0] += gkappau[0]; } // curvature effects

        ggd[0] = 0; 
        ggd[1] = guratio[0] * ( tan( gdelta ) * gpbx[0] * gh[0] + sico.GRAVITY * pow( guvw[0][0], 2 ) / sflow.TUFRI );
        ggd[2] = gvratio[0] * ( tan( gdelta ) * gpby[0] * gh[0] + sico.GRAVITY * pow( guvw[0][0], 2 ) / sflow.TUFRI );

        for ( go=3; go<sico.NVECTMIN+18; go++ ) ggd[go] = 0;
    }

    return ggd;
}

float fsigma ( float gv, float gvm, float gvp, int gdir, float *gdx, float *gdy, int gi, struct ico sico ) { // function for slope of vector

    float gsigma, gd = 0;

    if ( gdir == 1 ) gd = gdx[gi]; // cell spacing in appropriate direction
    else if ( gdir == 2 ) gd = gdy[gi];

    if ( gvp == gv ) gsigma = 0;
    else {

        gsigma = ( gv - gvm ) / ( gvp - gv ); // input for limiter

        if ( sico.LIMITER == 1 ) gsigma = ffmax( 0, ffmin( 1, gsigma )); // Minmod limiter

        else if ( sico.LIMITER == 2 ) {

            gsigma = ffmax( 0, ffmin( 1, 2 * gsigma ));
            gsigma = ffmax( gsigma, ffmin( gsigma, 2 )); // Superbee limiter
            
        } else if ( sico.LIMITER == 3 ) {

            gsigma = ffmin( 2 * gsigma, 0.5 * ( 1 + gsigma ));
            gsigma = ffmax( 0, ffmin ( 2, gsigma )); // Woodward limiter
            
        } else if ( sico.LIMITER == 4 ) gsigma = ( gsigma + fabs( gsigma )) / ( 1 + fabs( gsigma )); // van Leer limiter

        gsigma *= ( gvp - gv ) / gd; // slope of vector
    }
    return gsigma;
}

float fcvelr2 ( float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowy, float gvflowx2, float gvflowy2,
    float gvflowx3, float gvflowy3, struct flow sflow, float *gkx, float *gky, float *gbetaxy, int gi, struct ico sico ) { // function for determining time step length

    int gl;
    float gcelx = 0, gcely = 0, gcel, gcelr, gh, galpha[3], ggamma[3], ggamma2[3], gu[3], gv[3];

    gh = ghflow + ghflow2 + ghflow3;
    gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
    gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

    if ( gh > sico.HFLOWMIN ) {

        if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
        if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
        if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        if ( galpha[0] > 0  && galpha[2] > 0 ) ggamma[0] = sflow.RHO3 / sflow.RHO1; else ggamma[0] = 0;
        if ( galpha[0] > 0  && galpha[1] > 0 ) ggamma[1] = sflow.RHO2 / sflow.RHO1; else ggamma[1] = 0;
        if ( galpha[1] > 0  && galpha[2] > 0 ) ggamma[2] = sflow.RHO3 / sflow.RHO2; else ggamma[2] = 0;

        if ( galpha[0] > 0  && galpha[2] > 0 ) ggamma2[0] = sflow.RHO3 / sflow.RHO1; else ggamma2[0] = 1;
        if ( galpha[0] > 0  && galpha[1] > 0 ) ggamma2[1] = sflow.RHO2 / sflow.RHO1; else ggamma2[1] = 1;
        if ( galpha[1] > 0  && galpha[2] > 0 ) ggamma2[2] = sflow.RHO3 / sflow.RHO2; else ggamma2[2] = 1;

        gcel = 0;

        for ( gl = 0; gl < sico.PMAX; gl++ ) {

            if ( gl == 0 && sico.MODEL == 0 ) {

                gcelx = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * gkx[gl] * ghflow, 0.5 );
                gcely = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * gky[gl] * ghflow, 0.5 );

            } else if ( sico.PHASES[gl] == 1 && gl == 0 ) {

                gcelx  = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * ( 1 - ggamma[0] ) * gkx[gl] * ( gh * ( galpha[gl] + 0.5 * ( 1 - galpha[gl] ))), 0.5 );
                gcely  = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * ( 1 - ggamma[0] ) * gky[gl] * ( gh * ( galpha[gl] + 0.5 * ( 1 - galpha[gl] ))), 0.5 );

            } else if ( sico.PHASES[gl] == 1 && gl == 1 ) {

                gcelx  = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * ( 1 - ggamma[2] ) * gkx[gl] * ( gh * ( galpha[gl] + 0.5 * ( 1 - galpha[gl] ))), 0.5 );
                gcely  = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * ( 1 - ggamma[2] ) * gky[gl] * ( gh * ( galpha[gl] + 0.5 * ( 1 - galpha[gl] ))), 0.5 );

            } else if ( sico.PHASES[gl] == 2 ) {

                gcelx = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * ggamma2[2] * ( gh * ( galpha[gl] + 0.5 * ( 1 - galpha[gl] ))), 0.5 );
                gcely = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * ggamma2[2] * ( gh * ( galpha[gl] + 0.5 * ( 1 - galpha[gl] ))), 0.5 );

            } else if ( sico.PHASES[gl] == 3 ) {

                gcelx = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * ( gh * ( galpha[gl] + 0.5 * ( 1 - galpha[gl] ))), 0.5 );
                gcely = pow( sico.GRAVITY * fabs( cos( gbetaxy[gi] )) * ( gh * ( galpha[gl] + 0.5 * ( 1 - galpha[gl] ))), 0.5 );
            }

            gcelx = fabs( gu[gl] ) + ffmin( gcelx, 0.0 + 20 * gcelx * fabs( -0.5 * tanh( 1.0 * gcelx  )));
            gcely = fabs( gv[gl] ) + ffmin( gcely, 0.0 + 20 * gcely * fabs( -0.5 * tanh( 1.0 * gcely  )));
            if ( gcelx > gcely ) gcelr = gcelx;  else gcelr =  gcely;
            if ( gcelr > gcel ) gcel = gcelr;
        }

    } else gcel = 0;

    return gcel;
}

float fconvout ( int ggi, float **gdflow, int gmaterial, int gtype, float gbetaxyi, struct ico sico ) { // function for converting depths into heights

    int gm1 = 0, gm2 = 0, gm3 = 0;
    float ghflowi = 0, gdflowi = 0, grflowi = 0;

    if ( sico.MODEL <= 3 ) {

        if ( gtype == 1 ) gm1 = 0;
        else if ( gtype == 2 ) gm1 = 3;
        else if ( gtype == 3 ) gm1 = 7;
        gdflowi = gdflow[ggi][gm1];

        if ( gdflowi != 0 ) ghflowi = gdflowi / cos ( gbetaxyi );
        else ghflowi = 0;

        grflowi = 1;

    } else if ( sico.MODEL == 7 ) {

        if ( gtype == 1 ) { gm1 = 0; gm2 = 3; gm3 = 6; }
        else if ( gtype == 2 ) { gm1 = 9; gm2 = 10; gm3 = 11; }
        else if ( gtype == 3 ) { gm1 = 25; gm2 = 27; gm3 = 29; }
        gdflowi = gdflow[ggi][gm1] + gdflow[ggi][gm2] + gdflow[ggi][gm3];

        if ( gdflowi != 0 ) {
            if ( gmaterial == 1 ) grflowi = gdflow[ggi][gm1] / gdflowi;
            else if ( gmaterial == 2 ) grflowi = gdflow[ggi][gm2] / gdflowi;
            else if ( gmaterial == 3 ) grflowi = gdflow[ggi][gm3] / gdflowi;
            else if ( gmaterial == 4 ) grflowi = 1;
            ghflowi = gdflowi * grflowi / cos ( gbetaxyi );
        }
        else ghflowi = 0;
    }

    if ( sico.CORRHEIGHT == 0 ) ghflowi = gdflowi * grflowi;

    return ghflowi;
}

void foutasc ( float **gparam, int *gpx, int *gpy, char *goutmaps, char *gname, float *gbetaxy, int gk, struct ico sico ) { // function for output of ascii raster maps

    FILE *gfascii;
    char gpath[200];
    int gi, gx, gy;
    float gpout;

    sprintf(gpath, "%s%s.asc", goutmaps, gname ); // writing name of ascii file to string
    gfascii=fopen(gpath, "w"); // opening ascii file

    fprintf( gfascii, "ncols %i\n", sico.N ); // writing header to ascii file
    fprintf( gfascii, "nrows %i\n", sico.M );
    fprintf( gfascii, "xllcenter %.2f\n", sico.BDWEST );
    fprintf( gfascii, "yllcenter %.2f\n", sico.BDSOUTH );
    fprintf( gfascii, "cellsize %.2f\n", sico.CSZ );
    fprintf( gfascii, "NODATA_value %.3f\n", sico.UNDEF );

    gi = 0;
    for ( gx=0; gx<sico.M; gx++ ) { for ( gy=0; gy<sico.N; gy++ ) {

        if ( gpx[gi] == gx && gpy[gi] == gy ) {

            if ( sico.MODEL <= 3 && gk == 0 ) gpout = fconvout( gi, gparam, 1, 1, gbetaxy[gi], sico ); // computing output cell value (height-to-depth conversion where necessary)
            else if ( sico.MODEL <= 3 && gk == 3 ) gpout = fconvout( gi, gparam, 1, 2, gbetaxy[gi], sico );
            else if ( sico.MODEL <= 3 && gk == 7 ) gpout = fconvout( gi, gparam, 1, 3, gbetaxy[gi], sico );

            else if ( sico.MODEL <= 3 && gk == 1 ) gpout = fdiv( gparam[gi][gk+1], gparam[gi][0], sico.HFLOWMIN );
            else if ( sico.MODEL <= 3 && gk == 2 ) gpout = fdiv( -gparam[gi][gk-1], gparam[gi][0], sico.HFLOWMIN );

            else if ( sico.MODEL == 7 && gk == 0 ) gpout = fconvout( gi, gparam, 1, 1, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 3 ) gpout = fconvout( gi, gparam, 2, 1, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 6 ) gpout = fconvout( gi, gparam, 3, 1, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 9 ) gpout = fconvout( gi, gparam, 1, 2, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 10 ) gpout = fconvout( gi, gparam, 2, 2, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 11 ) gpout = fconvout( gi, gparam, 3, 2, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 15 ) gpout = fconvout( gi, gparam, 4, 1, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 24 ) gpout = fconvout( gi, gparam, 4, 2, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 25 ) gpout = fconvout( gi, gparam, 1, 3, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 27 ) gpout = fconvout( gi, gparam, 2, 3, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 29 ) gpout = fconvout( gi, gparam, 3, 3, gbetaxy[gi], sico );
            else if ( sico.MODEL == 7 && gk == 31 ) gpout = fconvout( gi, gparam, 4, 3, gbetaxy[gi], sico );

            else if ( sico.MODEL == 7 && gk == 1 ) gpout = fdiv( gparam[gi][gk+1], gparam[gi][0], sico.HFLOWMIN );
            else if ( sico.MODEL == 7 && gk == 2 ) gpout = fdiv( -gparam[gi][gk-1], gparam[gi][0], sico.HFLOWMIN );
            else if ( sico.MODEL == 7 && gk == 4 ) gpout = fdiv( gparam[gi][gk+1], gparam[gi][3], sico.HFLOWMIN );
            else if ( sico.MODEL == 7 && gk == 5 ) gpout = fdiv( -gparam[gi][gk-1], gparam[gi][3], sico.HFLOWMIN );
            else if ( sico.MODEL == 7 && gk == 7 ) gpout = fdiv( gparam[gi][gk+1], gparam[gi][6], sico.HFLOWMIN );
            else if ( sico.MODEL == 7 && gk == 8 ) gpout = fdiv( -gparam[gi][gk-1], gparam[gi][6], sico.HFLOWMIN );

            else gpout = gparam[gi][gk];

            fprintf( gfascii, "%.3f ", gpout ); // writing output to ascii file
            gi += 1;
            
        } else fprintf( gfascii, "-9999.000 " ); // writing no data value to file if cell is not in area of interest
        
    } fprintf(gfascii, "\n"); }

    fclose(gfascii);

    return;
}

void foutdircoord ( FILE *gf_directions, int gintmode, int *gpx, int *gpy, struct ico sico ) { // function for coordinates for output flow direction file

    int gi, git, gy0 = 0, gx, gy, intv, intv2 = 0;
    float gxmetric, gymetric;

    intv = (int)( sico.N / 75 + 1 ); // cell interval for writing output

    if ( gintmode == 1 ) { intv2 = 1; gy0 = 0; } // defining shift depending on type of parameter
    else if ( gintmode == 2 ) { intv2 = 2; gy0 = 0; }
    else if ( gintmode == 3 ) { intv2 = 2; gy0 = intv; }
    else if ( gintmode == 4 ) { intv2 = 3; gy0 = 0; }
    else if ( gintmode == 5 ) { intv2 = 3; gy0 = intv; }
    else if ( gintmode == 6 ) { intv2 = 3; gy0 = 2*intv; }

    gi = 0;
    for ( gx=0; gx<sico.M; gx+=intv ) { 
      for ( gy=gy0; gy<sico.N; gy+=(intv2*intv)) {

        git = gi;
        while (( gpx[git] != gx || gpy[git] != gy ) && git < sico.IMAX ) git+= 1;
        if ( git < sico.IMAX ) {

            gi = git;
            gxmetric = ( float ) gpy[gi] * sico.CSZ + sico.BDWEST;
            fprintf( gf_directions, "%.0f\t", gxmetric ); // writing x coordinate to file
        }
      }

      if ( sico.MODEL <= 3 ) { if ( gy0 == intv ) gy0 = 0; else gy0 = intv; }
      else { if ( gy0 == 2*intv ) gy0 = intv; else if ( gy0 == intv ) gy0 = 0; else gy0 = 2*intv; }
    }

    fprintf( gf_directions, "\n" );

    gi = 0;
    for ( gx=0; gx<sico.M; gx+=intv ) { 
      for ( gy=gy0; gy<sico.N; gy+=(intv2*intv)) {

        git = gi;
        while (( gpx[git] != gx || gpy[git] != gy ) && git < sico.IMAX ) git+= 1;
        if ( git < sico.IMAX ) {

            gi = git;
            gymetric = sico.BDNORTH - ( float ) gpx[gi] * sico.CSZ;
            fprintf( gf_directions, "%.0f\t", gymetric ); // writing y coordinate to file
        }
      }

      if ( sico.MODEL <= 3 ) { if ( gy0 == intv ) gy0 = 0; else gy0 = intv; }
      else { if ( gy0 == 2*intv ) gy0 = intv; else if ( gy0 == intv ) gy0 = 0; else gy0 = 2*intv; }
    }

    fprintf( gf_directions, "\n" );

    return;
}

void foutdir ( FILE *gf_directions, float **gparam, int gl, int gintmode, int *gpx, int *gpy, struct ico sico ) { // function for output flow direction file

    int gi, git, gy0 = 0, gx, gy, intv, intv2 = 0;
    float ghflowi = 0, gparami = 0;

    intv = (int)( sico.N / 75 + 1 ); // cell interval for writing output

    if ( gintmode == 1 ) { intv2 = 1; gy0 = 0; } // defining shift depending on type of parameter
    else if ( gintmode == 2 ) { intv2 = 2; gy0 = 0; }
    else if ( gintmode == 3 ) { intv2 = 2; gy0 = intv; }
    else if ( gintmode == 4 ) { intv2 = 3; gy0 = 0; }
    else if ( gintmode == 5 ) { intv2 = 3; gy0 = intv; }
    else if ( gintmode == 6 ) { intv2 = 3; gy0 = 2*intv; }

    gi = 0;
    for ( gx=0; gx<sico.M; gx+=intv ) {
      for ( gy=gy0; gy<sico.N; gy+=(intv2*intv)) {

        git = gi;
        while (( gpx[git] != gx || gpy[git] != gy ) && git < sico.IMAX ) git+= 1;
        if ( git < sico.IMAX ) {

            gi = git;

            if ( gintmode == 1 ) ghflowi = gparam[gi][0];
            else if ( gintmode == 2 ) ghflowi = gparam[gi][0];
            else if ( gintmode == 3 ) ghflowi = gparam[gi][3];
            else if ( gintmode == 4 ) ghflowi = gparam[gi][0];
            else if ( gintmode == 5 ) ghflowi = gparam[gi][3];
            else if ( gintmode == 6 ) ghflowi = gparam[gi][6];

            if ( ghflowi >= sico.IMPTHR[0] ) {

                if ( gintmode == 1 || gintmode == 2 ) gparami = gparam[gi][gl] / gparam[gi][0];
                else if ( gintmode == 3 ) gparami = gparam[gi][gl] / gparam[gi][3];
                else if ( gintmode == 4 ) gparami = gparam[gi][gl] / gparam[gi][0];
                else if ( gintmode == 5 ) gparami = gparam[gi][gl] / gparam[gi][3];
                else if ( gintmode == 6 ) gparami = gparam[gi][gl] / gparam[gi][6];
                fprintf( gf_directions, "%.2f\t", gparami );
            }
            else fprintf( gf_directions, "0\t" ); // writing parameter to file
        }
      }

      if ( sico.MODEL <= 3 ) { if ( gy0 == intv ) gy0 = 0; else gy0 = intv; }
      else { if ( gy0 == 2*intv ) gy0 = intv; else if ( gy0 == intv ) gy0 = 0; else gy0 = 2*intv; }
    }

    fprintf( gf_directions, "\n" );

    return;
}


#ifdef WITHGRASS


    SEGMENT finrasti ( char *gm, struct ico sico ) { // input of integer GRASS raster maps

        int gf, gx, gy, gin;
        CELL *gc;
        SEGMENT gseg_in;

        Segment_open ( &gseg_in, G_tempfile(), sico.M, sico.N, sico.NSEGRC, sico.NSEGRC, sizeof(int), sico.NSEGS );

        gc = Rast_allocate_c_buf();
        gf = Rast_open_old ( gm, sico.MAINMAPSET );

        for ( gx = 0; gx < sico.M; gx++ ) {

            Rast_get_c_row ( gf, gc, gx );

            for ( gy = 0; gy < sico.N; gy++ ) {

                if ( !Rast_is_c_null_value ( gc + gy ) ) gin = ( int ) gc[gy];
                else gin = (int)sico.UNDEF;
                Segment_put(&gseg_in, &gin, gx, gy);
	    }
        }
        Segment_flush(&gseg_in);
        G_free ( gc );
        Rast_close ( gf );

        return gseg_in;
    }

    SEGMENT finrastd ( char *gm, struct ico sico ) { // input of float GRASS raster maps

        int gf, gx, gy;
        float gin;
        DCELL *gc;
        SEGMENT gseg_in;

        Segment_open ( &gseg_in, G_tempfile(), sico.M, sico.N, sico.NSEGRC, sico.NSEGRC, sizeof(float), sico.NSEGS );

        gc = Rast_allocate_d_buf();
        gf = Rast_open_old ( gm, sico.MAINMAPSET );

        for ( gx = 0; gx < sico.M; gx++ ) {

            Rast_get_d_row ( gf, gc, gx );

            for ( gy = 0; gy < sico.N; gy++ ) {

                if ( !Rast_is_d_null_value ( gc + gy ) ) gin = ( float ) gc[gy];
                else gin = sico.UNDEF;
                Segment_put(&gseg_in, &gin, gx, gy);
	    }
        }

        Segment_flush(&gseg_in);
        G_free ( gc );
        Rast_close ( gf );

        return gseg_in;
    }

    void foutrast ( char *gm, float ***goutv, struct ico sico, int gk, float gtsum ) { // output of GRASS raster maps

        int gf, gx, gy;
        DCELL *gc;

        if ( gm != NULL ) {
            gc = Rast_allocate_d_buf();
            gf = Rast_open_new( gm, DCELL_TYPE );

            for ( gx = 0; gx < sico.M; gx++ ) {

                if ( gm != NULL ) {

	            for ( gy = 0; gy < sico.N; gy++ ) {

                        gc[gy] = (DCELL) goutv[gx][gy][gk];
                    }
	            Rast_put_d_row( gf, gc );
                }
            }
            Rast_close( gf );
            G_free( gc );

        }

        return;
    }


#else


    float *finaschdr ( char *gname ) { // function for input of ascii raster map header

        FILE *gfascii;
        char gpath[200], *tok;
        char delim[] = " \t\n";
        int gi, gx, gy;
        float gtest3, *gparam, gtestthrsh;
        char *gtest, *gtest2;
        char **strtodhlp=NULL;

        gparam = (float*) calloc( 7, sizeof(float));

        sprintf(gpath, "%s", gname ); // writing name of ascii file to string
        gfascii=fopen(gpath, "r"); // opening ascii file

        for ( gi=0; gi<6; gi++ ) {
            gtest = fcparam ( gfascii );
            tok=strtok(gtest, delim);
            while ( tok != NULL ) {
                gtest2=tok;
                tok=strtok(NULL, delim);
            }
            gparam[gi] = strtod (gtest2, strtodhlp); // reading header lines and writing to array
            free(gtest);
        }

        gtestthrsh = gparam[5]; // threshold for validity

        gi = 0;
        for ( gx=0; gx<(int)gparam[1]; gx++ ) { // counting number of cells

            gtest = flparam ( gfascii );

            tok = strtok(gtest, delim);
            for ( gy=0; gy<(int)gparam[0]; gy++ ) {
                gtest2 = tok;
                tok = strtok(NULL, delim);
                gtest3 = strtod (gtest2, strtodhlp);
                if ( gtest3 >= gtestthrsh ) gi += 1;
            }
            free(gtest);
        }
        gparam[6] = (float)gi;
        fclose( gfascii );
        return gparam;
    }

    int *finascxy ( char *gname, int gcontrol, struct ico sico ) { // function for input of ascii raster x and y coordinates

        FILE *gfascii;
        char *gpath, *tok;
        char delim[] = " \t\n";
        int gi, gx, gy, *gpxy;
        float gtest3, gtestthrsh;
        char **strtodhlp=NULL;

        char *gtest;
        char *gtest2;

        gpath = (char*) calloc ( 1000, sizeof(char));
        gpxy = (int*) calloc ( sico.IMAX, sizeof(int));

        gtestthrsh = sico.UNDEF; // threshold for validity

        sprintf(gpath, "%s", gname ); // writing name of ascii file to string
        gfascii=fopen(gpath, "r"); // opening ascii file
        for ( gi=0; gi<6; gi++ ) {
            gtest = fcparam ( gfascii ); // reading header lines
            free(gtest);
        }

        gi = 0;
        for ( gx=0; gx<sico.M; gx+=sico.MESH ) {

            gtest = flparam ( gfascii );
            tok = strtok(gtest, delim);
            for ( gy=0; gy<sico.N; gy+=sico.MESH ) {

                gtest2 = tok;
                tok = strtok(NULL, delim);
                gtest3 = strtod (gtest2, strtodhlp);

                if ( gtest3 >= gtestthrsh ) {

                    if ( gcontrol == 1 ) gpxy[gi] = gx; // writing internal coordinates of all rows and columns to array
                    else gpxy[gi] = gy;

                    gi += 1;
                }
            }
            free(gtest);
        }
        fclose( gfascii );
        free(gpath);
        return gpxy;
    }

    float *finascval ( char *gname, int *gpx, int *gpy, struct ico sico ) { // function for input of ascii raster map values

        FILE *gfascii;
        char *gpath, *tok;
        char delim[] = " \t\n";
        int gi, gx, gy;
        float gtest3, *gparam;
        char **strtodhlp=NULL;

        char *gtest;
        char *gtest2;

        gpath = (char*) calloc ( 1000, sizeof(char));
        gparam = (float*) calloc( sico.IMAX, sizeof(float));

        sprintf(gpath, "%s", gname ); // writing name of ascii file to string
        gfascii=fopen(gpath, "r"); // opening ascii file
        for ( gi=0; gi<6; gi++ ) {
            gtest = fcparam ( gfascii ); // reading header lines
            free(gtest);
        }

        gi = 0;
        for ( gx=0; gx<sico.M; gx+=sico.MESH ) {

            gtest = flparam ( gfascii );
            tok = strtok(gtest, delim);
            for ( gy=0; gy<sico.N; gy+=sico.MESH ) {

                gtest2 = tok;
                tok = strtok(NULL, delim);
                gtest3 = strtod(gtest2, strtodhlp);

                if ( gx == gpx[gi] && gy == gpy[gi] ) {
                    gparam[gi] = gtest3; // writing raster values of all rows and columns to array
                    gi += 1;
                }
            }
            free(gtest);
        }
        fclose( gfascii );
        free(gpath);
        return gparam;
    }


#endif


// -- STOP --- Functions ----------------------------------------------------------------------------------------


int main ( int argc, char *argv[] ) { // calling main function


// -- START -- Declaration of variables -------------------------------------------------------------------------


    FILE *fparam;

    char *prefix, *outmaps, *outfiles, *elevname,
        *hreleasename = 0, *hreleasename2 = 0, *hreleasename3 = 0, *hentrmaxname = 0, *hentrmaxname2 = 0, *hentrmaxname3 = 0, *vinxname = 0, *vinyname = 0, *vinxname2 = 0, *vinyname2 = 0,
        *vinxname3 = 0, *vinyname3 = 0, *zonesname = 0, *centrname = 0, *cvshearname = 0, *phiname = 0, *phi2name = 0, *phi3name = 0, *deltabname = 0, *tufriname = 0, *deltaname = 0, 
        *delta2name = 0, *delta3name = 0, *nyssname = 0, *nyfsname = 0, *nyffname = 0, *ambdragname = 0, *flufriname = 0, *transssfsname = 0, *transssffname = 0, *transfsffname = 0, 
        *hydreleasename0 = 0, *frictioname = 0, *transformoname = 0, *treleasename = 0, *trelstopname = 0, *stoptimename = 0, *tslidename = 0, 
        *mv, **mv0; char *madd = (char*) calloc(1000, sizeof(char));

    char *path = (char*) calloc(1000, sizeof(char));
    char *wkdir = (char*) calloc(800, sizeof(char));

    int *px, *py, *iin, *inn, *innn, hydrograph, frictiograph, transformograph, hydnin, hydnout, *hydi = 0, hydj = 0, hydk = 0, *hydx = 0, *hydy = 0, hydt = 0, *hydtmax = 0,
        hydtmaxx = 0, *hydp0 = 0, **hydp = 0, fritmax = 0, frik, frit = 0, tratmax = 0, trak, tratx = 0, nvect_all, nvect_red, lmax, xint, ccontinue, csuccess, nsum, nout, time_start,
        time_stop, ctrlr, ctrlv, ctrlvv, ctrlvvv, ccfl, fj[4][2], i, j, jj, jjj, k, l, ll, p, z, prec_hflow, prec_vol, prec_ekin, nzones,
        imax = 0, iloop, ctrl_hydout, ix, i2, hydcols, ctrl_lake, ctrl_noosc, cflowpre, anctr, anid, andist, cslide, cflow, iy, iz, anwin0;

    float *pelev, gkx[3], gky[3], *kappau, *vm, *cdrag, *disp, *gze, *gf, *gg, *gs, *flowpar, tout, tmax, vflowx = 0, vflowy = 0, vflow, vflowx1, vflowy1, vflowx2, vflowy2, 
        vflowx3, vflowy3, time_elapsed, tlength, tlength0, tlengthx, tlengthpre, tint, tsum, cfl, cflmax, grav[5], dw_dt[9], hflow_maxmax, hflow_max0, hflow_max, hflow_max02, 
        hflow_max2, hflow_max03, hflow_max3, vflow_max, vflow_max2, vflow_max3, **hydhyd0 = 0, ***hydhyd = 0, **frifri = 0, **tratra = 0, vol_flow0, vol_flow, vol_flow02, vol_flow2, vol_flow03,
        vol_flow3, vol_entr, vol_entr2, vol_entr3, vol_edge, vol_edge2, vol_edge3, whx, why, vcelr0, vcelr, *hydelev = 0, *hydalpha = 0, *hydx_metric = 0, *hydy_metric = 0, *hydl = 0, 
        hhyd0 = 0, hhyd = 0, hyde = 0, hyde2 = 0, hyde3 = 0, hydfalpha = 0, hydout_xmin_metric = 0, hydout_xmax_metric = 0, hydout_ymin_metric = 0, hydout_ymax_metric = 0, 
        hydh = 0, hydh2 = 0, hydh3 = 0, hydv = 0, hydv2 = 0, hydv3 = 0, hydq = 0, hydq2 = 0, hydq3 = 0, hydlmax = 0, ekin, ekin_flow, vol_noflux, vol_noflux2, vol_noflux3, carea = 0, 
        hydmalpha, hydnalpha, hydbeta, hydm0, hydmx, hydmy, hydm, hydfcorr, wu[9], wv[9], wu2[9], wv2[9], wu3[9], wv3[9], nbetax[9], nbetay[9], walpha[9], walpha2[9], walpha3[9], welev[9], 
        vflowxj[9], vflowyj[9], corrfact, qentr, qentr1, qentr2, qentr3, qentrtest, 
        alpha, alphav, mom, walphax, walphay, walphax2, walphay2, walphax3, walphay3, wbetax, wbetay, wbetaxy,
        wdx, wdy, wgrav[5], vol_hydbef, vol_hydaft, pelevhydtest, vol_hydbef2, vol_hydaft2, vol_hydbef3, vol_hydaft3, vol_hyd, vol_hyd2, vol_hyd3,
        *gdecel, hflow = 0, hflown = 0, trat, ttrelease, ttrelstop, tfact, trans, wwd[18], wdu[3], wdv[3], xwdu[3], xwdv[3], wh[9], whflow = 0, whflow2 = 0, whflow3 = 0,
        dw_dttest[9], wintbtest, wintctest, awttest[9], hekin, hekin_sum, hekin_max, hydbef[9], minlevel, wmx = 0, wmy = 0, dux, duy, duxy, dumain,
        rhrem, qvol_test, qh_test, qtrelstart = 0, qtrelstop = 0, qtrelspan, ctrl_pressthr, vmax, elevtot[9], betav, lambdam, lambdab,
        alpha1 = 0, alpha2 = 0, alpha3 = 0, alphab1 = 0, alphab2 = 0, alphab3 = 0, rhom, rhob, gammam, gammab, gz, 
        momaddsx, momaddfsx, momaddfx, momaddsy, momaddfsy, momaddfy, mym, myb, alphasfs, alphabsfs, qentrup, qentrdown, momfact = 0, momfacttest = 0, hflowj[9], hentrmaxx,
        phreleaseall = 0, qhreleaseall, vol_flow0all, andelta, anslopex, anslopey, anslope, angx, angy, angz, 
        anpx, anpy, anpx1, anpy1, anpx2, anpy2, anwht[4], ansumh = 0, ansumslopex = 0, ansumslopey = 0, anavgslopex, anavgslopey, 
        andhdx = 0, andhdy = 0, vflowxratio, vflowyratio, anrad, anwhtd, anupx, anupy; // andh, anaspectdh, aneta, anzeta, hflows, anetax, anetay, anzetax, anzetay, analpha, anbeta, ansy0, anaspect, momx = 0, momy = 0, ancdv, anavgh, angamma, ank, andhdxy, analphas, angxy;

    struct ico sico;
    struct flow sflow;


    #ifdef WITHGRASS


        SEGMENT seg_elev, seg_hrelease, seg_hrelease2, seg_hrelease3, seg_vinx, seg_viny, seg_vinx2, seg_viny2, seg_vinx3, seg_viny3,
        seg_hentrmax, seg_hentrmax2, seg_hentrmax3, seg_zones, seg_centr, seg_cvshear, seg_phi, seg_phi2, seg_phi3, seg_deltab, seg_tufri, seg_delta, seg_delta2, seg_delta3, seg_nyss, seg_nyfs, 
        seg_nyff, seg_ambdrag, seg_flufri, seg_transssfs, seg_transssff, seg_transfsff, seg_trelease, seg_trelstop, seg_stoptime, seg_tslide;
        char *yint = (char*) calloc(1000, sizeof(char));
        int x, y, aoi, zones;
        float hflowi = 0, hflowi2 = 0, hflowi3 = 0, hentri = 0, hentri2 = 0, hentri3 = 0, pout, elev, hrelease, hrelease2, hrelease3, vinx, viny, vinx2, viny2, vinx3, viny3,
        hentrmax, hentrmax2, hentrmax3, centr, cvshear, phi, phi2, phi3, deltab, tufri, delta, delta2, delta3, nyss, nyfs, nyff, ambdrag, flufri,
        transssfs, transssff, transfsff, trelease, trelstop, stoptime, tslide;

        struct GModule *module;
        struct Cell_head cellhd;
        struct Option *input_opt1;


    #else


        float *aschdr;


    #endif


// -- STOP --- Declaration of variables -------------------------------------------------------------------------


    time_start = clock(); // start time


    #ifdef WITHGRASS


        printf("Starting model execution.\n\n");


    #else


        printf("Starting model execution.\n");


    #endif


    fflush(stdout); // forcing immediate display


// Preparing environment


    #ifdef WITHGRASS


        xint=1;
        yint=getenv("XINT");
        if ( yint != NULL ) xint=atoi(yint);
        
        G_gisinit( argv[0] );

        module = G_define_module();
        G_add_keyword(_("raster"));
        G_add_keyword(_("landslide"));
        G_add_keyword(_("numerical simulation"));
        module->description =
	    _("The mass flow simulation tool");

        sprintf( wkdir, "%s/%s/.tmp/rtemp", G_location_path(),G_mapset());
        sprintf( path, "%s/param%d.txt", wkdir, xint ); // writing name of parameter file to string

        input_opt1 = G_define_standard_option(G_OPT_F_BIN_INPUT);
        input_opt1->key = "input1";
        input_opt1->description = _("Parameter file 1");
        input_opt1->required = NO;
        input_opt1->answer = path;

        if (G_parser(argc, argv))
            exit(EXIT_FAILURE);

        G_get_set_window( &cellhd );
        
        G_init_tempfile();
        sico.NSEGRC = 16;
        sico.NSEGS = 16;

        fparam=fopen( input_opt1->answer, "r" ); // opening parameter file


    #else


        if ( argc > 1 ) sprintf( wkdir, "%s", argv[1] );
        else {
            printf( "ERROR: Please provide a parameter file.\n" );
            fflush(stdout);
            exit( EXIT_SUCCESS );
        }
        xint=1;

        sprintf( path, "%s/param%d.txt", wkdir, xint ); // writing name of parameter file to string
        fparam=fopen( path, "r" ); // opening parameter file


    #endif


// -- START -- Reading and preprocessing parameter file input ---------------------------------------------------


    sico.PI = 3.1415926536; // pi

    if( !fparam ) { // managing lack of parameter file
        printf( "ERROR: Unable to open parameter file: '%s'\n", path );
        fflush(stdout);
        exit( EXIT_SUCCESS );
    }

    sico.MAINMAPSET = fcparam ( fparam ); // name of main mapset (relevant with GRASS)
    sico.MULT = fiparam ( fparam ); // identifier for single model run (0) or multiple model runs (1)
    sico.MESH = fiparam ( fparam ); // mesh spacing
    sico.MODEL = fiparam ( fparam ); // model (0=mixture, 1=one-phase, 7=multi-phase)

    sico.PMAX = fiparam ( fparam ); // number of phases
    for ( l=0; l<sico.PMAX; l++ ) sico.PHASES[l] = fiparam ( fparam ); // phases

    sico.AFLAG = fiparam ( fparam ); // control for additional output raster maps
    sico.GRAVITY = fdparam ( fparam ); // gravity
    sico.LIMITER = fiparam ( fparam ); // numerical limiter (1=Minmod, 2=Superbee, 3=Woodward, 4=van Leer)
    
    sico.CORRHEIGHT = fiparam ( fparam ); // conversion control (0=no, 1=yes)
    sico.CURVCTRL = fiparam ( fparam ); // curvature control (0=no, 1=decelerating, 2=all, 3-5=combined with diffusion control)
    sico.NOOSC = fiparam ( fparam ); // surface control (0=no, 1-6=different combinations)
    sico.ENTRAINMENT = fiparam ( fparam ); // entrainment and deposition control (0=no, 1-4=type of approach)
    sico.STOPPING = fiparam ( fparam ); // stopping control (0=no, 1-3=type of approach)
    sico.DYNFRIC = fiparam ( fparam ); // friction control (0=no, 1=yes)
    sico.BETACTRL = 1; // control for adapting slope to shifted coordinate system (0=no, 1=yes)

    if ( sico.NOOSC < 0 || sico.NOOSC == 8 ) sico.SEPFLUX = 1; else sico.SEPFLUX = 0; // identifier whether to use phase separation (1) or not (0)
    if ( sico.NOOSC < 0 ) sico.NOOSC *= -1;
    if ( sico.NOOSC == 8 ) sico.NOOSC = 0;

    if ( sico.MODEL <= 3 ) { // size of arrays of state variables

        nvect_all = 19;
        nvect_red = 7;
        sico.NVECTMIN = 3;
        
    } else {

        nvect_all = 48;
        nvect_red = 25;
        sico.NVECTMIN = 9;
    }

    sico.RELM = fiparam ( fparam ); // identifier for use of raster maps (1) or not (0)
    if ( sico.MODEL == 7 ) sico.RELM2 = fiparam ( fparam );
    if ( sico.MODEL == 7 ) sico.RELM3 = fiparam ( fparam );
    sico.RELV = fiparam ( fparam );
    if ( sico.MODEL == 7 ) sico.RELV2 = fiparam ( fparam );
    if ( sico.MODEL == 7 ) sico.RELV3 = fiparam ( fparam );
    sico.ENTR = fiparam ( fparam );
    if ( sico.MODEL == 7 ) sico.ENTR2 = fiparam ( fparam );
    if ( sico.MODEL == 7 ) sico.ENTR3 = fiparam ( fparam );
    sico.ZONES = fiparam ( fparam );
    sico.CENTR = fiparam ( fparam );
    sico.CVSHEAR = fiparam ( fparam );
    sico.PHI = fiparam ( fparam );
    sico.PHI2 = fiparam ( fparam );
    sico.PHI3 = fiparam ( fparam );
    sico.DELTAB = fiparam ( fparam );
    sico.TUFRI = fiparam ( fparam );
    sico.DELTA = fiparam ( fparam );
    sico.DELTA2 = fiparam ( fparam );
    sico.DELTA3 = fiparam ( fparam );
    sico.NYSS = fiparam ( fparam );
    sico.NYFS = fiparam ( fparam );
    sico.NYFF = fiparam ( fparam );
    sico.AMBDRAG = fiparam ( fparam );
    sico.FLUFRI = fiparam ( fparam );
    sico.TRANSSSFS = fiparam ( fparam );
    sico.TRANSSSFF = fiparam ( fparam );
    sico.TRANSFSFF = fiparam ( fparam );
    sico.TRELEASE = fiparam ( fparam );
    sico.TRELSTOP = fiparam ( fparam );
    sico.STOPTIME = fiparam ( fparam );
    sico.TSLIDE = fiparam ( fparam );

    hydrograph = fiparam ( fparam ); // identifier for use of hydrograph(s), frictiograph, and transformograph (1) or not (0)
    frictiograph = fiparam ( fparam );
    transformograph = fiparam ( fparam );

    elevname = fcparam ( fparam ); // name of elevation map


// Setting spatial extent and cell size


    #ifdef WITHGRASS // if GRASS rasters are used


        sico.CSZ = cellhd.ew_res; // cell size
        sico.BDWEST = cellhd.west; // western boundary
        sico.BDNORTH = cellhd.north; // northern boundary
        sico.BDSOUTH = cellhd.south; // southern boundary
        sico.N = cellhd.cols; // number of cells in y direction (W-E)
        sico.M = cellhd.rows; // number of cells in x direction (S-N)
        sico.UNDEF = -9999; // no data value


    #else // if ascii rasters are used


        aschdr = finaschdr( elevname ); // reading header of ascii raster
        sico.N = (int)aschdr[0] / sico.MESH; // number of cells in y direction (W-E)
        sico.M = (int)aschdr[1] / sico.MESH; // number of cells in x direction (S-N)
        sico.CSZ = aschdr[4] * sico.MESH; // cell size
        sico.BDWEST = aschdr[2]; // western boundary
        sico.BDSOUTH = aschdr[3]; // southern boundary
        sico.BDNORTH = aschdr[3] + sico.M * sico.CSZ; // northern boundary
        sico.UNDEF = aschdr[5]; // no data value
        sico.IMAX = (int)aschdr[6] / pow( sico.MESH, 2 ); // number of cells

        free(aschdr);


    #endif


    // Raster maps

    if ( sico.RELM == 1 ) hreleasename = fcparam ( fparam ); // name of mixture or PHASE 1 release height map
    if ( sico.MODEL == 7 && sico.RELM2 == 1 ) hreleasename2 = fcparam ( fparam ); // name of PHASE 2 release height map
    if ( sico.MODEL == 7 && sico.RELM3 == 1 ) hreleasename3 = fcparam ( fparam ); // name of PHASE 3 release height map
    if ( sico.RELV == 1 ) vinxname = fcparam ( fparam ); // name of mixture or PHASE 1 release velocity in x direction map
    if ( sico.MODEL == 7 && sico.RELV2 == 1 ) vinxname2 = fcparam ( fparam ); // name of PHASE 2 release velocity in x direction map
    if ( sico.MODEL == 7 && sico.RELV3 == 1 ) vinxname3 = fcparam ( fparam ); // name of PHASE 3 release velocity in x direction map
    if ( sico.RELV == 1 ) vinyname = fcparam ( fparam ); // name of mixture or PHASE 1 release velocity in y direction map
    if ( sico.MODEL == 7 && sico.RELV2 == 1 ) vinyname2 = fcparam ( fparam ); // name of PHASE 2 release velocity in y direction map
    if ( sico.MODEL == 7 && sico.RELV3 == 1 ) vinyname3 = fcparam ( fparam ); // name of PHASE 3 release velocity in y direction map
    if ( sico.ENTR == 1 ) hentrmaxname = fcparam ( fparam ); // name of maximum height of mixture or PHASE 1 entrainment map
    if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) hentrmaxname2 = fcparam ( fparam ); // name of maximum height of PHASE 2 entrainment map
    if ( sico.MODEL == 7 && sico.ENTR3 == 1 ) hentrmaxname3 = fcparam ( fparam ); // name of maximum height of PHASE 3 entrainment map
    if ( sico.ZONES == 1 ) zonesname = fcparam ( fparam ); // name of zones map
    if ( sico.CENTR == 1 ) centrname = fcparam ( fparam ); // name of entrainment coefficient map
    if ( sico.CVSHEAR == 1 ) cvshearname = fcparam ( fparam ); // name of shear velocity coefficient map
    if ( sico.PHI == 1 ) phiname = fcparam ( fparam ); // name of internal friction angle of mixture or PHASE 1 map
    if ( sico.PHI2 == 1 ) phi2name = fcparam ( fparam ); // name of internal friction angle of PHASE 2 map
    if ( sico.PHI3 == 1 ) phi3name = fcparam ( fparam ); // name of internal friction angle of PHASE 3 map
    if ( sico.DELTAB == 1 ) deltabname = fcparam ( fparam ); // name of basal friction difference map
    if ( sico.TUFRI == 1 ) tufriname = fcparam ( fparam ); // name of turbulent friction coefficient map
    if ( sico.DELTA == 1 ) deltaname = fcparam ( fparam ); // name of basal friction angle of mixture or PHASE 1 map
    if ( sico.DELTA2 == 1 ) delta2name = fcparam ( fparam ); // name of basal friction angle of PHASE 2 map
    if ( sico.DELTA3 == 1 ) delta3name = fcparam ( fparam ); // name of basal friction angle of PHASE 3 map
    if ( sico.NYSS == 1 ) nyssname = fcparam ( fparam ); // name of kinematic viscosity of PHASE 1 map
    if ( sico.NYFS == 1 ) nyfsname = fcparam ( fparam ); // name of kinematic viscosity of PHASE 2 map
    if ( sico.NYFF == 1 ) nyffname = fcparam ( fparam ); // name of kinematic viscosity of PHASE 3 map
    if ( sico.AMBDRAG == 1 ) ambdragname = fcparam ( fparam ); // name of ambient drag map
    if ( sico.FLUFRI == 1 ) flufriname = fcparam ( fparam ); // name of fluid friction number map
    if ( sico.TRANSSSFS == 1 ) transssfsname = fcparam ( fparam ); // name of PHASE 1 - PHASE 2 transformation coefficient map
    if ( sico.TRANSSSFF == 1 ) transssffname = fcparam ( fparam ); // name of PHASE 1 - PHASE 3 transformation coefficient map
    if ( sico.TRANSFSFF == 1 ) transfsffname = fcparam ( fparam ); // name of PHASE 2 - PHASE 3 transformation coefficient map
    if ( sico.TRELEASE == 1 ) treleasename = fcparam ( fparam ); // name of release time map
    if ( sico.TRELSTOP == 1 ) trelstopname = fcparam ( fparam ); // name of release stop time map
    if ( sico.STOPTIME == 1 ) stoptimename = fcparam ( fparam ); // name of stopping time map
    if ( sico.TSLIDE == 1 ) tslidename = fcparam ( fparam ); // name of time of initial sliding map


    // Hydrograph, frictiograph, and transformograph

    if ( hydrograph == 1 ) {

        hydnin = fiparam ( fparam ); // number of input hydrographs
        hydnout = fiparam ( fparam ); // number of output hydrographs

        hydtmax = (int*) calloc( hydnin, sizeof(int)); // allocating memory to variable for number(s) of hydrograph time steps
        hydi = (int*) calloc(( hydnin + hydnout ), sizeof(int));
            // allocating memory to variable for internal coordinate of hydrograph point
        hydelev = (float*) calloc(( hydnin + hydnout ), sizeof(float));
            // allocating memory to variable for elevation of hydrograph point
        hydalpha = (float*) calloc(( hydnin + hydnout ), sizeof(float));
            // allocating memory to variable for aspect of hydrograph point
        hydx_metric = (float*) calloc(( hydnin + hydnout ), sizeof(float));
            // allocating memory to variable for x coordinate of hydrograph
        hydy_metric = (float*) calloc(( hydnin + hydnout ), sizeof(float));
            // allocating memory to variable for y coordinate of hydrograph
        hydl = (float*) calloc(( hydnin + hydnout ), sizeof(float)); // allocating memory to variable for length of hydrograph profile
        hydx = (int*) calloc(( hydnin + hydnout ), sizeof(int)); // allocating memory to variable for internal x coordinate of hydrograph
        hydy = (int*) calloc(( hydnin + hydnout ), sizeof(int)); // allocating memory to variable for internal y coordinate of hydrograph

        hydtmaxx = 0; // initializing maximum number of hydrograph time steps
    }
    else hydnin = 1;

    char *hydreleasename[hydnin]; // declaring pointers to input hydrograph files

    if ( hydrograph == 1 ) {

        for ( hydj = 0; hydj < hydnin; hydj++ ) { // loop over all input hydrographs

            hydreleasename[hydj] = (char*) calloc( 1000 * hydnin, sizeof(char)); // allocating memory to pointers to input hydrograph files

            hydreleasename0 = fcparam ( fparam ); // name of hydrograph file
            strcpy( hydreleasename[hydj], hydreleasename0 );
            free( hydreleasename0 );
            hydtmax[hydj] = fiparam ( fparam ); // number of hydrograph time steps
            if ( hydtmax[hydj] > hydtmaxx ) hydtmaxx = hydtmax[hydj]; // updating maximum number of hydrograph time steps
        }

        for ( hydj = 0; hydj < hydnin + hydnout; hydj++ ) { // loop over all hydrographs

            hydx_metric[hydj] = fdparam ( fparam ); // x coordinate of hydrograph
            hydy_metric[hydj] = fdparam ( fparam ); // y coordinate of hydrograph
            hydl[hydj] = fdparam ( fparam ); // length of hydrograph profile
            hydalpha[hydj] = fdparam ( fparam ) * sico.PI / 180; // aspect of hydrograph profile
        }

        hydhyd = alloc_dmatrix3( hydtmaxx+1, 7, hydnin ); // allocating memory to array of input hydrograph data
    }
    else { hydnin = 1; hydnout = 1; }

    if ( frictiograph == 1 ) {

        frictioname = fcparam ( fparam ); // name of frictiograph file
        fritmax = fiparam ( fparam ); // number of frictiograph time steps
        frifri = alloc_dmatrix( fritmax+1, 7 ); // allocating memory to array of input frictiograph data
    }

    if ( transformograph == 1 ) {

        transformoname = fcparam ( fparam ); // name of transformograph file
        tratmax = fiparam ( fparam ); // number of transformograph time steps
        tratra = alloc_dmatrix( tratmax+1, 7 ); // allocating memory to array of input transformograph data
    }


    // Flow parameters

    lmax = fiparam ( fparam ); // number of flow parameters
    flowpar = (float*) calloc( lmax, sizeof(float)); // allocating memory for storage of flow parameter values
    for ( l=0; l<lmax; l++ ) flowpar[l] = fdparam ( fparam ); // flow parameters

    if ( sico.MODEL == 0 ) {

        sflow.RHO1 = flowpar[0]; // density
        
        sflow.PHI0[0] = flowpar[1] * sico.PI / 180; // internal friction angle
        sflow.DELTA0[0] = flowpar[2] * sico.PI / 180; // basal friction angle
        sflow.TUFRI0 = pow(10, flowpar[3]); // turbulent friction
        
        sflow.CENTR0 = flowpar[4]; // entrainment coefficient
        sflow.CVSHEAR0 = flowpar[5]; // shear velocity coefficient
        sflow.DELTAB0 = flowpar[6]; // basal friction difference
        sflow.CSTOP = flowpar[7]; // stopping criterion

        sflow.KPMAX[0] = 9999; // maximum of earth pressure coefficient (hard-coded)
        sico.HYDADD = flowpar[8]; //2; identifier whether to impose hydrograph upon flow (1), reset flow before putting hydrograph (0), or impose entire discharge on centre of hydrograph (2)
        sico.ORIGINAL = flowpar[9]; //0; identifier whether to suppress change of direction (0) or not (1)

        sico.FRICMIN = flowpar[10] * sico.PI / 180; //0.0;
        sico.EKINCOEF = pow( 10, flowpar[11] ); //-6.0;
        sico.FRI_exp = 0; //0.0;

    } else if ( sico.MODEL == 1 ) {

        sflow.RHO1 = flowpar[0]; // PHASE 1 density

        sflow.PHI0[0] = flowpar[1] * sico.PI / 180; // internal friction angle
        sflow.DELTA0[0] = flowpar[2] * sico.PI / 180; // basal friction angle
        sflow.FLUFRI0 = flowpar[3]; // fluid friction number

        sflow.NY0[0] = pow(10, flowpar[4]); // viscosity
        sflow.TAUY[0] = flowpar[5]; // yield strength

        sflow.CENTR0 = flowpar[6]; // entrainment coefficient
        sflow.CVSHEAR0 = flowpar[7]; // shear velocity coefficient
        sflow.DELTAB0 = flowpar[8]; // basal friction difference
        sflow.CSTOP = flowpar[9]; // stopping criterion

        sflow.AMBDRAG0 = flowpar[10]; //0.0;

        sflow.VIS_chi[0] = flowpar[11]; //1;
        sflow.VIS_ry = flowpar[12]; //10;
        sflow.NY_exp = 0; //0;
        
        sflow.KPMAX[0] = flowpar[13]; //1;

        sico.HYDADD = flowpar[14]; //2; identifier whether to impose hydrograph upon flow (1), reset flow before putting hydrograph (0), or impose entire discharge on centre of hydrograph (2)
        sico.ORIGINAL = flowpar[15]; //0; identifier whether to suppress change of direction (0) or not (1)

        sico.FRICMIN = flowpar[16] * sico.PI / 180; //0.0;
        sico.EKINCOEF = pow( 10, flowpar[17] ); //-6.0;
        sico.FRI_exp = 0; //0.0;

    } else if ( sico.MODEL == 7 ) {

        sflow.RHO1 = flowpar[0]; // PHASE 1 density
        sflow.RHO2 = flowpar[1]; // PHASE 2 density
        sflow.RHO3 = flowpar[2]; // PHASE 3 density

        sflow.PHI0[0] = flowpar[3] * sico.PI / 180; // internal friction angle of PHASE 1
        sflow.DELTA0[0] = flowpar[4] * sico.PI / 180; // basal friction angle of PHASE 1
        sflow.PHI0[1] = flowpar[5] * sico.PI / 180; // internal friction angle of PHASE 2
        sflow.DELTA0[1] = flowpar[6] * sico.PI / 180; // basal friction angle of PHASE 2
        sflow.PHI0[2] = flowpar[7] * sico.PI / 180; // internal friction angle of PHASE 3
        sflow.DELTA0[2] = flowpar[8] * sico.PI / 180; // basal friction angle of PHASE 3
        sflow.FLUFRI0 = flowpar[9]; // fluid friction number

        sflow.NY0[0] = pow(10, flowpar[10]); // viscosity of PHASE 1
        sflow.TAUY[0] = flowpar[11]; // yield strength of PHASE 1
        sflow.NY0[1] = pow(10, flowpar[12]); // viscosity of PHASE 2
        sflow.TAUY[1] = flowpar[13]; // yield strength of PHASE 2
        sflow.NY0[2] = pow(10, flowpar[14]); // viscosity of PHASE 3
        sflow.TAUY[2] = flowpar[15]; // yield strength of PHASE 3

        sflow.CENTR0 = flowpar[16]; // entrainment coefficient
        sflow.CVSHEAR0 = flowpar[17]; // shear velocity coefficient
        sflow.DELTAB0 = flowpar[18]; // basal friction difference
        sflow.THETAS = 1 / ( 1 - flowpar[19] ) - 1; // maximum water content of deposit
        sflow.CSTOP = flowpar[20]; // stopping criterion

        sflow.TRANSSSFS0 = flowpar[21]; // transformation coefficient PHASE 1 - PHASE 2
        sflow.TRANSSSFF0 = flowpar[22]; // transformation coefficient PHASE 1 - PHASE 3
        sflow.TRANSFSFF0 = flowpar[23]; // transformation coefficient PHASE 2 - PHASE 3

        sflow.AMBDRAG0 = flowpar[24]; //0.0;

        sflow.VM_n0 = flowpar[25]; //10;
        sflow.VM_l = flowpar[26]; //0.12;
        sflow.VM_n = flowpar[27]; //1;
        sflow.VM_fact = flowpar[28]; //1;

        sflow.DRAG_k = flowpar[29]; //1;
        sflow.DRAG_m = flowpar[30]; //3;
        sflow.DRAG_n = flowpar[31]; //1;
        sflow.DRAG_vterm = flowpar[32]; //0.1;
        sflow.DRAG_rep = flowpar[33]; //1;
        sflow.DRAG_je = flowpar[34]; //1;

        sflow.VIS_chi[0] = flowpar[35]; //1;
        sflow.VIS_chi[1] = flowpar[36]; //1;
        sflow.VIS_chi[2] = flowpar[37]; //1;
        sflow.VIS_xi[0] = flowpar[38]; //0;
        sflow.VIS_xi[1] = flowpar[39]; //0;
        sflow.VIS_xi[2] = flowpar[40]; //0;
        sflow.VIS_a[0] = flowpar[41]; //1;
        sflow.VIS_a[1] = flowpar[42]; //1;
        sflow.VIS_a[2] = flowpar[43]; //1;
        sflow.VIS_ry = flowpar[44]; //10;
        sflow.NY_exp = flowpar[45]; //0;
        
        sflow.KPMAX[0] = flowpar[46]; //1;
        sflow.KPMAX[1] = flowpar[47]; //1;
        sflow.KPMAX[2] = flowpar[48]; //1;
        
        sico.HYDADD = flowpar[49]; //2; identifier whether to impose hydrograph upon flow (1), reset flow before putting hydrograph (0),
        //or impose entire discharge on centre of hydrograph (2)
        sico.ORIGINAL = flowpar[50]; //0; identifier whether to suppress change of direction (0) or not (1)

        sico.FRICMIN = flowpar[51] * sico.PI / 180; //0.0;
        sico.EKINCOEF = pow( 10, flowpar[52] ); //-6.0;
        sico.FRI_exp = flowpar[53]; //0.0;
    }

    sflow.RHOB1 = sflow.RHO1; // densities of basal layer (set equal to densities of corresponding flow phases)
    sflow.RHOB2 = sflow.RHO2;
    sflow.RHOB3 = sflow.RHO3;

    sflow.sep_SHEAR = 0.0; // shear factor for phase separation
    sflow.disp_MULT = 1.0; // multiplication factor for dispersion   
    sflow.CCONST = 1.0; // coefficient for constraining dispersion and phase separation


    // Thresholds and further parameters

    sico.HFLOWMIN = fdparam ( fparam ); // threshold flow depth for simulation
    sico.IMPTHR[0] = fdparam ( fparam ); // threshold for display of flow height
    sico.IMPTHR[1] = fdparam ( fparam ); // threshold for display of flow kinetic energy
    sico.IMPTHR[2] = fdparam ( fparam ); // threshold for display of flow pressure

    tout = fdparam ( fparam ); // time interval for writing output
    tmax = fdparam ( fparam ); // time at which to stop simulation
    sico.SLOMO = fdparam ( fparam ); // factor for slow motion
    sico.GRAVITY = sico.GRAVITY / pow( sico.SLOMO, 2 ); // adapting gravity to slow motion

    sico.SLIDERAD = fdparam ( fparam ); // search radius for initial sliding
    sico.SLIDEEXP = fdparam ( fparam ); // weighting exponent for initial sliding
    sico.SLIDEDEF = fdparam ( fparam ); // coefficient for deformation during initial sliding
    anwin0 = (int)( sico.SLIDERAD / sico.CSZ + 0.5 ); // search radius for initial sliding in cells

    trat = round( tmax / tout );
    tmax = trat * tout; // ensuring that time to stop is multiple of output time step length

    sico.CFL[0] = fdparam ( fparam ); // CFL criterion
    sico.CFL[1] = fdparam ( fparam ); // length of first time step

    prefix = fcparam ( fparam ); // prefix for output
    outmaps = fcparam ( fparam ); // path and prefix for storing output maps
    outfiles = fcparam ( fparam ); // path and prefix for storing output files

    fclose(fparam); // closing parameter file


// -- STOP --- Reading and preprocessing parameter file input ---------------------------------------------------


// -- START -- Reading hydrograph, frictiograph, and transformograph data ---------------------------------------


    if ( hydrograph == 1 ) {

        if ( sico.MODEL <= 3 ) hydcols = 3; else hydcols = 7;

        hydlmax = 0;
        for ( hydj = 0; hydj < hydnin + hydnout; hydj++ ) if ( fabs(hydl[hydj]) > hydlmax ) hydlmax = fabs(hydl[hydj]);

        for ( hydj = 0; hydj < hydnin; hydj++ ) {
            hydhyd0 = finhyd( hydreleasename[hydj], hydtmax[hydj], hydtmaxx ); // reading data from all input hydrographs
            for ( i=0; i<hydtmaxx+1; i++ ) { for ( j=0; j<hydcols; j++ ) hydhyd[i][j][hydj] = hydhyd0[i][j]; }
            free( hydhyd0[0] ); free( hydhyd0 );
        }

        for ( hydj = 0; hydj < hydnin + hydnout; hydj++ ) { // internal x and y coordinates of all hydrographs

            hydy[hydj] = (int)(( hydx_metric[hydj] - sico.BDWEST ) / sico.CSZ + 0.5 );
            hydx[hydj] = (int)(( sico.BDNORTH - hydy_metric[hydj] ) / sico.CSZ + 0.5 );
        }
    }

    if ( frictiograph == 1 ) {

        frifri = finhyd( frictioname, fritmax, fritmax ); // reading data from input frictiograph
    }

    if ( transformograph == 1 ) {

        tratra = finhyd( transformoname, tratmax, tratmax ); // reading data from input transformograph
    }


// -- STOP --- Reading hydrograph, frictiograph, and transformograph data ---------------------------------------


// -- START -- Reading input raster maps (GRASS) ----------------------------------------------------------------


    #ifdef WITHGRASS


        seg_elev = finrastd( elevname, sico );
        
        if ( sico.RELM == 1 ) seg_hrelease = finrastd( hreleasename, sico );
        if ( sico.MODEL == 7 && sico.RELM2 == 1 ) seg_hrelease2 = finrastd( hreleasename2, sico );
        if ( sico.MODEL == 7 && sico.RELM3 == 1 ) seg_hrelease3 = finrastd( hreleasename3, sico );
        
        if ( sico.RELV == 1 ) seg_vinx = finrastd( vinxname, sico );
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) seg_vinx2 = finrastd( vinxname2, sico );
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) seg_vinx3 = finrastd( vinxname3, sico );
        
        if ( sico.RELV == 1 ) seg_viny = finrastd( vinyname, sico );
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) seg_viny2 = finrastd( vinyname2, sico );
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) seg_viny3 = finrastd( vinyname3, sico );
        
        if ( sico.ENTR == 1 ) seg_hentrmax = finrastd( hentrmaxname, sico );
        if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) seg_hentrmax2 = finrastd( hentrmaxname2, sico );
        if ( sico.MODEL == 7 && sico.ENTR3 == 1 ) seg_hentrmax3 = finrastd( hentrmaxname3, sico );

        if ( sico.ZONES == 1 ) seg_zones = finrasti( zonesname, sico );
        if ( sico.CENTR == 1 ) seg_centr = finrastd( centrname, sico );
        if ( sico.CVSHEAR == 1 ) seg_cvshear = finrastd( cvshearname, sico );
        if ( sico.PHI == 1 ) seg_phi = finrastd( phiname, sico );
        if ( sico.PHI2 == 1 ) seg_phi2 = finrastd( phi2name, sico ); 
        if ( sico.PHI3 == 1 ) seg_phi3 = finrastd( phi3name, sico );
        if ( sico.DELTAB == 1 ) seg_deltab = finrastd( deltabname, sico );
        if ( sico.TUFRI == 1 ) seg_tufri = finrastd( tufriname, sico );
        if ( sico.DELTA == 1 ) seg_delta = finrastd( deltaname, sico );
        if ( sico.DELTA2 == 1 ) seg_delta2 = finrastd( delta2name, sico );
        if ( sico.DELTA3 == 1 ) seg_delta3 = finrastd( delta3name, sico );
        if ( sico.NYSS == 1 ) seg_nyss = finrastd( nyssname, sico );
        if ( sico.NYFS == 1 ) seg_nyfs = finrastd( nyfsname, sico );
        if ( sico.NYFF == 1 ) seg_nyff = finrastd( nyffname, sico );
        if ( sico.AMBDRAG == 1 ) seg_ambdrag = finrastd( ambdragname, sico );
        if ( sico.FLUFRI == 1 ) seg_flufri = finrastd( flufriname, sico );
        if ( sico.TRANSSSFS == 1 ) seg_transssfs = finrastd( transssfsname, sico );
        if ( sico.TRANSSSFF == 1 ) seg_transssff = finrastd( transssffname, sico );
        if ( sico.TRANSFSFF == 1 ) seg_transfsff = finrastd( transfsffname, sico );
        if ( sico.TRELEASE == 1 ) seg_trelease = finrastd( treleasename, sico );
        if ( sico.TRELSTOP == 1 ) seg_trelstop = finrastd( trelstopname, sico );
        if ( sico.STOPTIME == 1 ) seg_stoptime = finrastd( stoptimename, sico );
        if ( sico.TSLIDE == 1 ) seg_tslide = finrastd( tslidename, sico );

        sico.IMAX = 0;
        for ( x = 0; x < sico.M; x++ ) {
            for ( y = 0; y < sico.N; y++ ) sico.IMAX += 1; // number of raster cells
        }


    #endif


// -- STOP --- Reading input raster maps (GRASS) ----------------------------------------------------------------


// -- START -- Preparing arrays for input -----------------------------------------------------------------------


    int **cedge = 0, **cedge2 = 0, **cedge3 = 0, *cedge0 = 0, *cedge02 = 0, *cedge03 = 0, *cneighbours = 0, *cneighbours2 = 0, *cneighbours3 = 0, *pzones;
    float *phrelease = 0, *phrelease2 = 0, *phrelease3 = 0, *qhrelease = 0, *qhrelease2 = 0, *qhrelease3 = 0, *pvinx = 0, *pvinx2 = 0, *pvinx3 = 0, *pviny = 0, *pviny2 = 0, *pviny3 = 0,
    *phentrmax = 0, *phentrmax2 = 0, *phentrmax3 = 0, *pcentr = 0, *pcvshear = 0, *pphi = 0, *pphi2 = 0, *pphi3 = 0, *pdeltab = 0, *ptufri = 0, *pdelta = 0, *pdelta2 = 0, *pdelta3 = 0, 
    *pnyss = 0, *pnyfs = 0, *pnyff = 0, *pambdrag = 0, *pflufri = 0, *ptransssfs = 0, *ptransssff = 0, *ptransfsff = 0, *ptrelease = 0, *ptrelstop = 0, *pstoptime = 0, *ptslide = 0, 
    **cready = 0, **cready2 = 0, **cready3 = 0, qtinit[hydnin+hydnout];

    int *cdomain = (int*) calloc( sico.IMAX, sizeof(int));
    int *cdomain2 = (int*) calloc( sico.IMAX, sizeof(int));
    int *cflux = (int*) calloc( sico.IMAX, sizeof(int));
    int *cstopped = (int*) calloc( sico.IMAX, sizeof(int));
    int *pxslide = (int*) calloc( sico.IMAX, sizeof(int));

    float *betax = (float*) calloc( sico.IMAX, sizeof(float));
    float *betay = (float*) calloc( sico.IMAX, sizeof(float));
    float *betaxy = (float*) calloc( sico.IMAX, sizeof(float));
    float *pelev0 = (float*) calloc( sico.IMAX, sizeof(float));
    float *anx = (float*) calloc( sico.IMAX, sizeof(float));
    float *anu = (float*) calloc( sico.IMAX, sizeof(float));

    if ( sico.ZONES != 1 ) pzones = (int*) calloc( sico.IMAX, sizeof(int));
    

    #ifdef WITHGRASS


        float *v;
        float ***outv = alloc_dmatrix3( sico.M, sico.N, nvect_all );

        v = (float*) calloc( nvect_all, sizeof(float));
        px = (int*) calloc( sico.IMAX, sizeof(int));
        py = (int*) calloc( sico.IMAX, sizeof(int));

        pelev = (float*) calloc( sico.IMAX, sizeof(float));
        phrelease = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRELEASE == 1 ) qhrelease = (float*) calloc( sico.IMAX, sizeof(float));
        pvinx = (float*) calloc( sico.IMAX, sizeof(float));
        pviny = (float*) calloc( sico.IMAX, sizeof(float));
        phentrmax = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.MODEL == 7 ) {

            phrelease2 = (float*) calloc( sico.IMAX, sizeof(float));
            if ( sico.TRELEASE == 1 ) qhrelease2 = (float*) calloc( sico.IMAX, sizeof(float));
            pvinx2 = (float*) calloc( sico.IMAX, sizeof(float));
            pviny2 = (float*) calloc( sico.IMAX, sizeof(float));
            phentrmax2 = (float*) calloc( sico.IMAX, sizeof(float));

            phrelease3 = (float*) calloc( sico.IMAX, sizeof(float));
            if ( sico.TRELEASE == 1 ) qhrelease3 = (float*) calloc( sico.IMAX, sizeof(float));
            pvinx3 = (float*) calloc( sico.IMAX, sizeof(float));
            pviny3 = (float*) calloc( sico.IMAX, sizeof(float));
            phentrmax3 = (float*) calloc( sico.IMAX, sizeof(float));
        }

        if ( sico.ZONES == 1 ) pzones = (int*) calloc( sico.IMAX, sizeof(int));
        if ( sico.CENTR == 1 ) pcentr = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.CVSHEAR == 1 ) pcvshear = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.PHI == 1 ) pphi = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.PHI2 == 1 ) pphi2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.PHI3 == 1 ) pphi3 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.DELTAB == 1 ) pdeltab = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TUFRI == 1 ) ptufri = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.DELTA == 1 ) pdelta = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.DELTA2 == 1 ) pdelta2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.DELTA3 == 1 ) pdelta3 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.NYSS == 1 ) pnyss = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.NYFS == 1 ) pnyfs = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.NYFF == 1 ) pnyff = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.AMBDRAG == 1 ) pambdrag = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.FLUFRI == 1 ) pflufri = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRANSSSFS == 1 ) ptransssfs = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRANSSSFF == 1 ) ptransssff = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRANSFSFF == 1 ) ptransfsff = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRELEASE == 1 ) ptrelease = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRELSTOP == 1 ) ptrelstop = (float*) calloc( sico.IMAX, sizeof(float));
        pstoptime = (float*) calloc( sico.IMAX, sizeof(float));
        ptslide = (float*) calloc( sico.IMAX, sizeof(float));


    #endif
   

    int **in = alloc_imatrix( sico.IMAX, 9 );
    int **ibasket = alloc_imatrix( 2, sico.IMAX );
    int **icheck = alloc_imatrix( sico.IMAX, 2 );
    int *ib = (int*) calloc( 2, sizeof(int));

    float *asigma_xelev = calloc( sico.IMAX, sizeof(float));
    float *asigma_yelev = calloc( sico.IMAX, sizeof(float));

    float **aw = alloc_dmatrix( sico.IMAX, nvect_all );
    float **awt = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **af = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **ag = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **as = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **ad = alloc_dmatrix( sico.IMAX, sico.NVECTMIN+18 );

    float **asigma_x = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **asigma_y = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **asigma_f = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **asigma_g = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );

    float **wintelev = alloc_dmatrix( sico.IMAX, 4 );
    float ***winta = alloc_dmatrix3( sico.IMAX, 4, sico.NVECTMIN );
    float ***wintb = alloc_dmatrix3( sico.IMAX, 6, sico.NVECTMIN );
    float ***wintc = alloc_dmatrix3( sico.IMAX, 4, sico.NVECTMIN );
    float **wintd = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **wintdtest = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );

    float **f = alloc_dmatrix( 4, sico.NVECTMIN );
    float **g = alloc_dmatrix( 4, sico.NVECTMIN );
    float **s = alloc_dmatrix( 4, sico.NVECTMIN );
    float ***d = alloc_dmatrix3( sico.IMAX, 6, sico.NVECTMIN+18 );

    float *dx = (float*) calloc( sico.IMAX, sizeof(float));
    float *dy = (float*) calloc( sico.IMAX, sizeof(float));
    
    if ( sico.CURVCTRL > 2 ) {

        cedge = alloc_imatrix( sico.IMAX, 9 );
        cready = alloc_dmatrix( sico.IMAX, 9 );
        cedge2 = alloc_imatrix( sico.IMAX, 9 );
        cready2 = alloc_dmatrix( sico.IMAX, 9 );
        cedge3 = alloc_imatrix( sico.IMAX, 9 );
        cready3 = alloc_dmatrix( sico.IMAX, 9 );
        cedge0 = (int*) calloc( sico.IMAX, sizeof(int));
        cedge02 = (int*) calloc( sico.IMAX, sizeof(int));
        cedge03 = (int*) calloc( sico.IMAX, sizeof(int));
        cneighbours = (int*) calloc( sico.IMAX, sizeof(int));
        cneighbours2 = (int*) calloc( sico.IMAX, sizeof(int));
        cneighbours3 = (int*) calloc( sico.IMAX, sizeof(int));
    }


// -- STOP --- Preparing arrays for input -----------------------------------------------------------------------


    mv0 = (char**) calloc( 1000, sizeof(char*));
    mv = (char*) calloc( 1000, sizeof(char));

    if ( sico.MODEL <= 3 ) {

        mv0[0] = "hflow"; mv0[1] = "vflowx"; mv0[2] = "vflowy"; mv0[3] = "basechange"; mv0[4] = "vflow";
        mv0[5] = "tflow"; mv0[6] = "pflow"; mv0[7] = "hflow_max"; mv0[8] = "vflow_max";
        mv0[9] = "tflow_max"; mv0[10] = "pflow_max"; mv0[11] = "vfront"; mv0[12] = "r1front"; mv0[13] = "r3front"; mv0[14] = "r1max"; mv0[15] = "r3max"; mv0[16] = "vhmax"; mv0[17] = "treach";
        
    } else if ( sico.MODEL == 7 ) {

        mv0[0]  = "hflow1";     mv0[1]  = "vflowx1";    mv0[2]  = "vflowy1";    mv0[3]  = "hflow2";      mv0[4]  = "vflowx2";     mv0[5] = "vflowy2";
        mv0[6]  = "hflow3";     mv0[7]  = "vflowx3";    mv0[8]  = "vflowy3";    mv0[9]  = "basechange1"; mv0[10] = "basechange2"; mv0[11] = "basechange3";
        mv0[12] = "vflow1";     mv0[13] = "vflow2";     mv0[14] = "vflow3";     mv0[15] = "hflow";       mv0[16] = "tflow1";      mv0[17] = "tflow2";
        mv0[18] = "tflow3";     mv0[19] = "tflow";      mv0[20] = "pflow1";     mv0[21] = "pflow2";      mv0[22] = "pflow3";      mv0[23] = "pflow";
        mv0[24] = "basechange"; mv0[25] = "hflow1_max"; mv0[26] = "vflow1_max"; mv0[27] = "hflow2_max";  mv0[28] = "vflow2_max";  mv0[29] = "hflow3_max";
        mv0[30] = "vflow3_max"; mv0[31] = "hflow_max";  mv0[32] = "tflow1_max"; mv0[33] = "pflow1_max";  mv0[34] = "tflow2_max";  mv0[35] = "pflow2_max";
        mv0[36] = "tflow3_max"; mv0[37] = "pflow3_max"; mv0[38] = "tflow_max";  mv0[39] = "pflow_max";   mv0[40] = "vfront";      mv0[41] = "r1front"; 
        mv0[42] = "r3front"; mv0[43] = "r1max"; mv0[44] = "r3max"; mv0[45] = "vhmax"; mv0[46] = "treach"; // variables for names of output raster maps
    }


// -- START -- Defining, opening, and pre-processing output files -----------------------------------------------


    FILE *f_summary, *f_volumes, *f_directions = 0, *f_directions2 = 0, *f_directions3 = 0, *f_nout, *f_hydout, *f_hydinfo[hydnin+hydnout], *f_hydtrans[hydnin+hydnout];

    if ( sico.MULT == 0 ) sprintf(path, "%s%ssummary.txt", outfiles, prefix); // summary file
    else sprintf(path, "%s%ssummary%d.txt", outfiles, prefix, xint);
    f_summary=fopen(path, "w");
    
    if ( sico.MULT == 0 ) sprintf(path, "%s%svolumes.txt", outfiles, prefix); // volumes file
    else sprintf(path, "%s%svolumes%d.txt", outfiles, prefix, xint);      
    f_volumes=fopen(path, "w");

    if ( hydrograph == 1 ) {

        sprintf(path, "%s%shydprofiles.txt", outfiles, prefix); // hydrograph profiles file
        f_hydout=fopen(path, "w");

        fprintf(f_hydout, "ID\tx1\txC\tx2\ty1\tyC\ty2\n"); // printing header of hydrograph profiles file

        if ( sico.MULT == 0 ) {

            for ( i = 0; i < hydnin + hydnout; i++ ) { // loop over all hydrographs

                sprintf(path, "%s%shydinfo%i.txt", outfiles, prefix, i+1); // hydrograph info file
                f_hydinfo[i]=fopen(path, "w");

                if ( i >= hydnin ) fprintf(f_hydinfo[i],
                    "T\tH1\tV1\tE1\tQ1\tH2\tV2\tE2\tQ2\tH3\tV3\tE3\tQ3\n0.0\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\n");
                else fprintf(f_hydinfo[i], "T\tH1\tV1\tE1\tQ1\tH2\tV2\tE2\tQ2\tH3\tV3\tE3\tQ3\n"); // printing header of hydrograph info file
                
                sprintf(path, "%s%shydtrans%i.txt", outfiles, prefix, i+1); // hydrograph transfer file
                f_hydtrans[i]=fopen(path, "w");

                fprintf(f_hydtrans[i], "T\tQ1\tV1\tQ2\tV2\tQ3\tV3\n"); // printing header of hydrograph transfer file
            }

        }
    }

    if ( sico.MULT == 0 ) {
    
        sprintf(path, "%s%sdirections1.txt", outfiles, prefix); // mixture or PHASE 1 flow direction file (for display as arrows)
        f_directions=fopen(path, "w");

        if ( sico.MODEL == 7 ) {
        
            sprintf(path, "%s%sdirections2.txt", outfiles, prefix); // PHASE 2 flow direction file (for display as arrows)
            f_directions2=fopen(path, "w");

            sprintf(path, "%s%sdirections3.txt", outfiles, prefix); // PHASE 3 flow direction file (for display as arrows)
            f_directions3=fopen(path, "w");
        }
    }


// -- STOP --- Defining, opening, and pre-processing output files -----------------------------------------------


// -- START -- Preparing arrays (if GRASS is used) --------------------------------------------------------------


    #ifdef WITHGRASS


        i = 0;

        for ( x = 0; x < sico.M; x++ ) { 
          for ( y = 0; y < sico.N; y++ ) {

            aoi = 1;
            Segment_get(&seg_elev, &elev, x, y); // reading data from segmentation files
            
            if ( sico.RELM == 1 ) Segment_get(&seg_hrelease, &hrelease, x, y);
            if ( sico.MODEL == 7 && sico.RELM2 == 1 ) Segment_get(&seg_hrelease2, &hrelease2, x, y);
            if ( sico.MODEL == 7 && sico.RELM3 == 1 ) Segment_get(&seg_hrelease3, &hrelease3, x, y);
            
            if ( sico.RELV == 1 ) Segment_get(&seg_vinx, &vinx, x, y);
            if ( sico.MODEL == 7 && sico.RELV2 == 1 ) Segment_get(&seg_vinx2, &vinx2, x, y);
            if ( sico.MODEL == 7 && sico.RELV3 == 1 ) Segment_get(&seg_vinx3, &vinx3, x, y);
            
            if ( sico.RELV == 1 ) Segment_get(&seg_viny, &viny, x, y);
            if ( sico.MODEL == 7 && sico.RELV2 == 1 ) Segment_get(&seg_viny2, &viny2, x, y);
            if ( sico.MODEL == 7 && sico.RELV3 == 1 ) Segment_get(&seg_viny3, &viny3, x, y);
            
            if ( sico.ENTR == 1 ) Segment_get(&seg_hentrmax, &hentrmax, x, y);
            if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) Segment_get(&seg_hentrmax2, &hentrmax2, x, y);
            if ( sico.MODEL == 7 && sico.ENTR3 == 1 ) Segment_get(&seg_hentrmax3, &hentrmax3, x, y);
       
            if ( sico.ZONES == 1 ) Segment_get(&seg_zones, &zones, x, y);     
            if ( sico.CENTR == 1 ) Segment_get(&seg_centr, &centr, x, y);
            if ( sico.CVSHEAR == 1 ) Segment_get(&seg_cvshear, &cvshear, x, y);
            if ( sico.PHI == 1 ) Segment_get(&seg_phi, &phi, x, y);
            if ( sico.PHI2 == 1 ) Segment_get(&seg_phi2, &phi2, x, y);
            if ( sico.PHI3 == 1 ) Segment_get(&seg_phi3, &phi3, x, y);
            if ( sico.DELTAB == 1 ) Segment_get(&seg_deltab, &deltab, x, y);
            if ( sico.TUFRI == 1 ) Segment_get(&seg_tufri, &tufri, x, y);
            if ( sico.DELTA == 1 ) Segment_get(&seg_delta, &delta, x, y);
            if ( sico.DELTA2 == 1 ) Segment_get(&seg_delta2, &delta2, x, y);
            if ( sico.DELTA3 == 1 ) Segment_get(&seg_delta3, &delta3, x, y);
            if ( sico.NYSS == 1 ) Segment_get(&seg_nyss, &nyss, x, y);
            if ( sico.NYFS == 1 ) Segment_get(&seg_nyfs, &nyfs, x, y);
            if ( sico.NYFF == 1 ) Segment_get(&seg_nyff, &nyff, x, y);
            if ( sico.AMBDRAG == 1 ) Segment_get(&seg_ambdrag, &ambdrag, x, y);
            if ( sico.FLUFRI == 1 ) Segment_get(&seg_flufri, &flufri, x, y);
            if ( sico.TRANSSSFS == 1 ) Segment_get(&seg_transssfs, &transssfs, x, y);
            if ( sico.TRANSSSFF == 1 ) Segment_get(&seg_transssff, &transssff, x, y);
            if ( sico.TRANSFSFF == 1 ) Segment_get(&seg_transfsff, &transfsff, x, y);
            if ( sico.TRELEASE == 1 ) Segment_get(&seg_trelease, &trelease, x, y);
            if ( sico.TRELSTOP == 1 ) Segment_get(&seg_trelstop, &trelstop, x, y);
            if ( sico.STOPTIME == 1 ) Segment_get(&seg_stoptime, &stoptime, x, y);
            if ( sico.TSLIDE == 1 ) Segment_get(&seg_tslide, &tslide, x, y);

            if ( aoi >= 1 ) {

                px[i] = x; // internal x and y coordinates of cell
                py[i] = y;

                pelev[i] = elev; pelev0[i] = elev;

                if ( sico.RELM == 1 && hrelease != sico.UNDEF ) phrelease[i] = hrelease; // writing data to arrays
                else phrelease[i] = 0;
                if ( sico.MODEL == 7 && sico.RELM2 == 1  && hrelease2 != sico.UNDEF ) phrelease2[i] = hrelease2;
                else if ( sico.MODEL == 7 ) phrelease2[i] = 0;
                if ( sico.MODEL == 7 && sico.RELM3 == 1  && hrelease3 != sico.UNDEF ) phrelease3[i] = hrelease3;
                else if ( sico.MODEL == 7 ) phrelease3[i] = 0;

                if ( sico.RELV == 1 && vinx != sico.UNDEF ) pvinx[i] = vinx;
                else pvinx[i] = 0;
                if ( sico.MODEL == 7 && sico.RELV2 == 1  && vinx2 != sico.UNDEF ) pvinx2[i] = vinx2;
                else if ( sico.MODEL == 7 ) pvinx2[i] = 0;
                if ( sico.MODEL == 7 && sico.RELV3 == 1  && vinx3 != sico.UNDEF ) pvinx3[i] = vinx3;
                else if ( sico.MODEL == 7 ) pvinx3[i] = 0;
                
                if ( sico.RELV == 1 && viny != sico.UNDEF ) pviny[i] = viny;
                else pviny[i] = 0;
                if ( sico.MODEL == 7 && sico.RELV2 == 1  && viny2 != sico.UNDEF ) pviny2[i] = viny2;
                else if ( sico.MODEL == 7 ) pviny2[i] = 0;
                if ( sico.MODEL == 7 && sico.RELV3 == 1  && viny3 != sico.UNDEF ) pviny3[i] = viny3;
                else if ( sico.MODEL == 7 ) pviny3[i] = 0;
                
                if ( sico.ENTR == 1 ) phentrmax[i] = hentrmax;
                else phentrmax[i] = 9999;
                if ( sico.MODEL == 7 && sico.ENTR2 == 1 && hentrmax2 != sico.UNDEF ) phentrmax2[i] = hentrmax2;
                else if ( sico.MODEL == 7 ) phentrmax2[i] = 9999;
                if ( sico.MODEL == 7 && sico.ENTR3 == 1 && hentrmax3 != sico.UNDEF ) phentrmax3[i] = hentrmax3;
                else if ( sico.MODEL == 7 ) phentrmax3[i] = 9999;

                if ( sico.ZONES == 1 && zones != sico.UNDEF ) pzones[i] = 0;
                else if ( sico.ZONES== 1 ) pzones[i] = 0;        
                if ( sico.CENTR == 1 && centr != sico.UNDEF ) pcentr[i] = centr;
                else if ( sico.CENTR == 1 ) pcentr[i] = sico.UNDEF;
                if ( sico.CVSHEAR == 1 && cvshear != sico.UNDEF ) pcvshear[i] = cvshear;
                else if ( sico.CVSHEAR == 1 ) pcvshear[i] = sico.UNDEF;
                if ( sico.PHI == 1 && phi != sico.UNDEF ) pphi[i] = phi * sico.PI / 180;
                else if ( sico.PHI == 1 ) pphi[i] = sico.UNDEF;
                if ( sico.PHI2 == 1 && phi2 != sico.UNDEF ) pphi2[i] = phi2 * sico.PI / 180;
                else if ( sico.PHI2 == 1 ) pphi2[i] = sico.UNDEF;
                if ( sico.PHI3 == 1 && phi3 != sico.UNDEF ) pphi3[i] = phi3 * sico.PI / 180;
                else if ( sico.PHI3 == 1 ) pphi3[i] = sico.UNDEF;
                if ( sico.DELTAB == 1 && deltab != sico.UNDEF ) pdeltab[i] = deltab * sico.PI / 180;
                else if ( sico.DELTAB == 1 ) pdeltab[i] = sico.UNDEF;
                if ( sico.TUFRI == 1 && tufri != sico.UNDEF ) ptufri[i] = pow( 10, tufri );
                else if ( sico.TUFRI == 1 ) ptufri[i] = sico.UNDEF;
                if ( sico.DELTA == 1 && delta != sico.UNDEF ) pdelta[i] = delta * sico.PI / 180;
                else if ( sico.DELTA == 1 ) pdelta[i] = sico.UNDEF;
                if ( sico.DELTA2 == 1 && delta2 != sico.UNDEF ) pdelta2[i] = delta2 * sico.PI / 180;
                else if ( sico.DELTA2 == 1 ) pdelta2[i] = sico.UNDEF;
                if ( sico.DELTA3 == 1 && delta3 != sico.UNDEF ) pdelta3[i] = delta3 * sico.PI / 180;
                else if ( sico.DELTA3 == 1 ) pdelta3[i] = sico.UNDEF;
                if ( sico.NYSS == 1 && nyss != sico.UNDEF ) pnyss[i] = pow( 10, nyss );
                else if ( sico.NYSS == 1 ) pnyss[i] = sico.UNDEF;
                if ( sico.NYFS == 1 && nyfs != sico.UNDEF ) pnyfs[i] = pow( 10, nyfs );
                else if ( sico.NYFS == 1 ) pnyfs[i] = sico.UNDEF;
                if ( sico.NYFF == 1 && nyff != sico.UNDEF ) pnyff[i] = pow( 10, nyff );
                else if ( sico.NYFF == 1 ) pnyff[i] = sico.UNDEF;
                if ( sico.AMBDRAG == 1 && ambdrag != sico.UNDEF ) pambdrag[i] = ambdrag;
                else if ( sico.AMBDRAG == 1 ) pambdrag[i] = sico.UNDEF;
                if ( sico.FLUFRI == 1 && flufri != sico.UNDEF ) pflufri[i] = flufri;
                else if ( sico.FLUFRI == 1 ) pflufri[i] = sico.UNDEF;
                if ( sico.TRANSSSFS == 1 && transssfs != sico.UNDEF ) ptransssfs[i] = transssfs;
                else if ( sico.TRANSSSFS == 1 ) ptransssfs[i] = sico.UNDEF;
                if ( sico.TRANSSSFF == 1 && transssff != sico.UNDEF ) ptransssff[i] = transssff;
                else if ( sico.TRANSSSFF == 1 ) ptransssff[i] = sico.UNDEF;
                if ( sico.TRANSFSFF == 1 && transfsff != sico.UNDEF ) ptransfsff[i] = transfsff;
                else if ( sico.TRANSFSFF == 1 ) ptransfsff[i] = sico.UNDEF;
                if ( sico.TRELEASE == 1 && trelease != sico.UNDEF ) ptrelease[i] = trelease;
                else if ( sico.TRELEASE == 1 ) ptrelease[i] = sico.UNDEF;
                if ( sico.TRELSTOP == 1 && trelstop != sico.UNDEF ) ptrelstop[i] = trelstop;
                else if ( sico.TRELSTOP == 1 ) ptrelstop[i] = sico.UNDEF;
                if ( sico.STOPTIME == 1 && stoptime != sico.UNDEF ) pstoptime[i] = stoptime;
                else pstoptime[i] = 5;
                if ( sico.TSLIDE == 1 && tslide != sico.UNDEF ) ptslide[i] = tslide;
                else ptslide[i] = sico.UNDEF;

                i += 1;
            }
          }
        }


// -- STOP --- Preparing arrays (if GRASS is used) --------------------------------------------------------------


    #else


// -- START -- Reading ascii rasters and pre-processing arrays (if GRASS is not used) ---------------------------


        px = finascxy( elevname, 1, sico ); // internal x and y coordinates of cell
        py = finascxy( elevname, 2, sico );

        pelev = finascval( elevname, px, py, sico ); // reading data from ascii rasters
        
        if ( sico.RELM == 1 ) phrelease = finascval( hreleasename, px, py, sico );
        else phrelease = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELM2 == 1 ) phrelease2 = finascval( hreleasename2, px, py, sico );
        else if ( sico.MODEL == 7 ) phrelease2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELM3 == 1 ) phrelease3 = finascval( hreleasename3, px, py, sico );
        else if ( sico.MODEL == 7 ) phrelease3 = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.RELV == 1 ) pvinx = finascval( vinxname, px, py, sico );
        else pvinx = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) pvinx2 = finascval( vinxname2, px, py, sico );
        else if ( sico.MODEL == 7 ) pvinx2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) pvinx3 = finascval( vinxname3, px, py, sico );
        else if ( sico.MODEL == 7 ) pvinx3 = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.RELV == 1 ) pviny = finascval( vinyname, px, py, sico );
        else pviny = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) pviny2 = finascval( vinyname2, px, py, sico ); 
        else if ( sico.MODEL == 7 ) pviny2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) pviny3 = finascval( vinyname3, px, py, sico );
        else if ( sico.MODEL == 7 ) pviny3 = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.ENTR == 1 ) phentrmax = finascval( hentrmaxname, px, py, sico );
        else phentrmax = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) phentrmax2 = finascval( hentrmaxname2, px, py, sico );
        else if ( sico.MODEL == 7 ) phentrmax2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.ENTR3 == 1 ) phentrmax3 = finascval( hentrmaxname3, px, py, sico );
        else if ( sico.MODEL == 7 ) phentrmax3 = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.ZONES == 1 ) pzones = finascval( zonesname, px, py, sico );
        if ( sico.CENTR == 1 ) pcentr = finascval( centrname, px, py, sico );
        if ( sico.CVSHEAR == 1 ) pcvshear = finascval( cvshearname, px, py, sico );
        if ( sico.PHI == 1 ) pphi = finascval( phiname, px, py, sico );
        if ( sico.PHI2 == 1 ) pphi2 = finascval( phi2name, px, py, sico );
        if ( sico.PHI3 == 1 ) pphi3 = finascval( phi3name, px, py, sico );
        if ( sico.DELTAB == 1 ) pdeltab = finascval( deltabname, px, py, sico );
        if ( sico.TUFRI == 1 ) ptufri = finascval( tufriname, px, py, sico );
        if ( sico.DELTA == 1 ) pdelta = finascval( deltaname, px, py, sico );
        if ( sico.DELTA2 == 1 ) pdelta2 = finascval( delta2name, px, py, sico );
        if ( sico.DELTA3 == 1 ) pdelta3 = finascval( delta3name, px, py, sico );
        if ( sico.NYSS == 1 ) pnyss = finascval( nyssname, px, py, sico );
        if ( sico.NYFS == 1 ) pnyfs = finascval( nyfsname, px, py, sico );
        if ( sico.NYFF == 1 ) pnyff = finascval( nyffname, px, py, sico );
        if ( sico.AMBDRAG == 1 ) pambdrag = finascval( ambdragname, px, py, sico );
        if ( sico.FLUFRI == 1 ) pflufri = finascval( flufriname, px, py, sico );
        if ( sico.TRANSSSFS == 1 ) ptransssfs = finascval( transssfsname, px, py, sico );
        if ( sico.TRANSSSFF == 1 ) ptransssff = finascval( transssffname, px, py, sico );
        if ( sico.TRANSFSFF == 1 ) ptransfsff = finascval( transfsffname, px, py, sico );
        if ( sico.TRELEASE == 1 ) ptrelease = finascval( treleasename, px, py, sico );
        if ( sico.TRELSTOP == 1 ) ptrelstop = finascval( trelstopname, px, py, sico );
        if ( sico.STOPTIME == 1 ) pstoptime = finascval( stoptimename, px, py, sico );
        else pstoptime = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TSLIDE == 1 ) ptslide = finascval( tslidename, px, py, sico );
        else ptslide = (float*) calloc( sico.IMAX, sizeof(float));

        for ( i=0; i<sico.IMAX; i++ ) { // finalizing input arrays

            pelev0[i] = pelev[i];
            
            if ( sico.RELM != 1 || phrelease[i] == sico.UNDEF ) phrelease[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELM2 != 1 || phrelease2[i] == sico.UNDEF )) phrelease2[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELM3 != 1 || phrelease3[i] == sico.UNDEF )) phrelease3[i] = 0;

            if ( sico.RELV != 1 || pvinx[i] == sico.UNDEF ) pvinx[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELV2 != 1 || pvinx2[i] == sico.UNDEF )) pvinx2[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELV3 != 1 || pvinx3[i] == sico.UNDEF )) pvinx3[i] = 0;

            if ( sico.RELV != 1 || pviny[i] == sico.UNDEF ) pviny[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELV2 != 1 || pviny2[i] == sico.UNDEF )) pviny2[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELV3 != 1 || pviny3[i] == sico.UNDEF )) pviny3[i] = 0;

            if ( sico.ENTR != 1 || phentrmax[i] == sico.UNDEF ) phentrmax[i] = 9999;
            if ( sico.MODEL == 7 && ( sico.ENTR2 != 1 || phentrmax2[i] == sico.UNDEF )) phentrmax2[i] = 9999;
            if ( sico.MODEL == 7 && ( sico.ENTR3 != 1 || phentrmax3[i] == sico.UNDEF )) phentrmax3[i] = 9999;

            if ( sico.ZONES == 1 && pzones[i] != sico.UNDEF ) pzones[i] = pzones[i];
            else if ( sico.ZONES == 1 ) pzones[i] = 0;
            if ( sico.CENTR == 1 && pcentr[i] != sico.UNDEF ) pcentr[i] = pcentr[i];
            else if ( sico.CENTR == 1 ) pcentr[i] = sico.UNDEF;
            if ( sico.CVSHEAR == 1 && pcvshear[i] != sico.UNDEF ) pcvshear[i] = pcvshear[i];
            else if ( sico.CVSHEAR == 1 ) pcvshear[i] = sico.UNDEF;
            if ( sico.PHI == 1 && pphi[i] != sico.UNDEF ) pphi[i] = pphi[i] * sico.PI / 180;
            else if ( sico.PHI == 1 ) pphi[i] = sico.UNDEF;
            if ( sico.PHI2 == 1 && pphi2[i] != sico.UNDEF ) pphi2[i] = pphi2[i] * sico.PI / 180;
            else if ( sico.PHI2 == 1 ) pphi2[i] = sico.UNDEF;
            if ( sico.PHI3 == 1 && pphi3[i] != sico.UNDEF ) pphi3[i] = pphi3[i] * sico.PI / 180;
            else if ( sico.PHI3 == 1 ) pphi3[i] = sico.UNDEF;
            if ( sico.DELTAB == 1 && pdeltab[i] != sico.UNDEF ) pdeltab[i] = pdeltab[i] * sico.PI / 180;
            else if ( sico.DELTAB == 1 ) pdeltab[i] = sico.UNDEF;
            if ( sico.TUFRI == 1 && ptufri[i] != sico.UNDEF ) ptufri[i] = pow( 10, ptufri[i] );
            else if ( sico.TUFRI == 1 ) ptufri[i] = sico.UNDEF;
            if ( sico.DELTA == 1 && pdelta[i] != sico.UNDEF ) pdelta[i] = pdelta[i] * sico.PI / 180;
            else if ( sico.DELTA == 1 ) pdelta[i] = sico.UNDEF;
            if ( sico.DELTA2 == 1 && pdelta2[i] != sico.UNDEF ) pdelta2[i] = pdelta2[i] * sico.PI / 180;
            else if ( sico.DELTA2 == 1 ) pdelta2[i] = sico.UNDEF;
            if ( sico.DELTA3 == 1 && pdelta3[i] != sico.UNDEF ) pdelta3[i] = pdelta3[i] * sico.PI / 180;
            else if ( sico.DELTA3 == 1 ) pdelta3[i] = sico.UNDEF;
            if ( sico.NYSS == 1 && pnyss[i] != sico.UNDEF ) pnyss[i] = pow( 10, pnyss[i] );
            else if ( sico.NYSS == 1 ) pnyss[i] = sico.UNDEF;
            if ( sico.NYFS == 1 && pnyfs[i] != sico.UNDEF ) pnyfs[i] = pow( 10, pnyfs[i] );
            else if ( sico.NYFS == 1 ) pnyfs[i] = sico.UNDEF;
            if ( sico.NYFF == 1 && pnyff[i] != sico.UNDEF ) pnyff[i] = pow( 10, pnyff[i] );
            else if ( sico.NYFF == 1 ) pnyff[i] = sico.UNDEF;
            if ( sico.AMBDRAG == 1 && pambdrag[i] != sico.UNDEF ) pambdrag[i] = pambdrag[i];
            else if ( sico.AMBDRAG == 1 ) pambdrag[i] = sico.UNDEF;
            if ( sico.FLUFRI == 1 && pflufri[i] != sico.UNDEF ) pflufri[i] = pflufri[i];
            else if ( sico.FLUFRI == 1 ) pflufri[i] = sico.UNDEF;
            if ( sico.TRANSSSFS == 1 && ptransssfs[i] != sico.UNDEF ) ptransssfs[i] = ptransssfs[i];
            else if ( sico.TRANSSSFS == 1 ) ptransssfs[i] = sico.UNDEF;
            if ( sico.TRANSSSFF == 1 && ptransssff[i] != sico.UNDEF ) ptransssff[i] = ptransssff[i];
            else if ( sico.TRANSSSFF == 1 ) ptransssff[i] = sico.UNDEF;
            if ( sico.TRANSFSFF == 1 && ptransfsff[i] != sico.UNDEF ) ptransfsff[i] = ptransfsff[i];
            else if ( sico.TRANSFSFF == 1 ) ptransfsff[i] = sico.UNDEF;
            if ( sico.TRELEASE == 1 && ptrelease[i] != sico.UNDEF ) ptrelease[i] = ptrelease[i];
            else if ( sico.TRELEASE == 1 ) ptrelease[i] = sico.UNDEF;
            if ( sico.TRELSTOP == 1 && ptrelstop[i] != sico.UNDEF ) ptrelstop[i] = ptrelstop[i];
            else if ( sico.TRELSTOP == 1 ) ptrelstop[i] = sico.UNDEF;
            if ( sico.STOPTIME == 1 && pstoptime[i] != sico.UNDEF ) pstoptime[i] = pstoptime[i];
            else pstoptime[i] = 5;
            if ( sico.TSLIDE == 1 && ptslide[i] != sico.UNDEF ) ptslide[i] = ptslide[i];
            else ptslide[i] = sico.UNDEF;
        }


    #endif


// -- STOP --- Reading ascii rasters and pre-processing arrays (if GRASS is not used) ---------------------------


// -- START -- Computing maximum release heights and release volumes --------------------------------------------


    hflow_max0 = 0; // initializing maximum mixture or PHASE 1 release height
    vol_flow0 = 0; // initializing mixture or PHASE 1 release volume
    
    if ( sico.MODEL == 7 ) { 
    
        hflow_max02 = 0; hflow_max03 = 0; // initializing maximum PHASE 2 and PHASE 3 release heights
        vol_flow02 = 0; vol_flow03 = 0; // initializing PHASE 2 and PHASE 3 release volumes
    }

    for ( i=0; i<sico.IMAX; i++ ) { // updating release height maxima and volumes:

        if ( sico.RELM == 1 ) {
            if ( phrelease[i] > hflow_max0 ) hflow_max0 = phrelease[i]; // updating maximum mixture or PHASE 1 release height
            vol_flow0 += phrelease[i] * pow ( sico.CSZ, 2 ); // updating mixture or PHASE 1 release volume
        }
        if ( sico.MODEL == 7 && sico.RELM2 == 1 ) {
            if ( phrelease2[i] > hflow_max02 ) hflow_max02 = phrelease2[i]; // updating maximum PHASE 2 release height
            vol_flow02 += phrelease2[i] * pow ( sico.CSZ, 2 ); // updating PHASE 2 release volume
        }
        if ( sico.MODEL == 7 && sico.RELM3 == 1 ) {
            if ( phrelease3[i] > hflow_max03 ) hflow_max03 = phrelease3[i]; // updating maximum PHASE 3 release height
            vol_flow03 += phrelease3[i] * pow ( sico.CSZ, 2 ); // updating PHASE 3 release volume
        }
    }
    
    
// -- STOP --- Computing maximum release heights and release volumes --------------------------------------------


// -- START -- Definition of computational domains --------------------------------------------------------------


    // *** cdomain[i] = 0: 1st row edge cells
    // *** cdomain[i] = 1: all other cells
    // *** cdomain2[i] = 1: 1st, 2nd or 3rd row edge cells
    // *** cdomain2[i] = 0: all other cells
    // *** Further refinement and dynamic adaptation of cdomain[i] at the beginning of each time step

    ib[0] = 0;

    for ( i=0; i<sico.IMAX; i++ ) {

        iin = fin( i, pelev, sico ); // cell neighbourhood
        for ( j=0; j<9; j++ ) in[i][j] = iin[j];
        ctrlv = 1; ctrlvv = 1; ctrlvvv = 1; // resetting controls
        for ( j=1; j<9; j++ ) { // loop over all neighbour cells

            if ( in[i][j] < 0 || in[i][j] >= sico.IMAX ) ctrlv = 0; // 1st row edge cells

            inn = fin( in[i][j], pelev, sico ); // neighbourhood of neighbour cell
            for ( jj=1; jj<9; jj++ ) { // loop over neighbourhood
                if ( inn[jj] < 0 || inn[jj] >= sico.IMAX ) ctrlvv = 0; // 2nd row edge cells

                innn = fin( inn[jj], pelev, sico ); // neighbourhood of neighbour cell
                for ( jjj=1; jjj<9; jjj++ ) // loop over neighbourhood
                    if ( innn[jjj] < 0 || innn[jjj] >= sico.IMAX ) ctrlvvv = 0; // 3rd row edge cells

                free( innn );
            }
            free( inn );
        }
        free( iin );

        icheck[i][0] = 0;
        icheck[i][1] = 0;

        cflux[i] = 1;

        if ( ctrlv == 0 ) { cdomain[i] = 0; cdomain2[i] = 1; } // setting domain controls
        else if ( ctrlvv == 0 || ctrlvvv == 0 ) { cdomain[i] = 1; cdomain2[i] = 1; }
        else { cdomain[i] = 1; cdomain2[i] = 0; }

        if ( cdomain2[i] == 0 ) {

            ibasket[0][ib[0]] = i;
            ib[0] += 1;
        }
    }
    
    nzones = 0; // zones for material budgets
    
    for ( i=0; i<sico.IMAX; i++ ) {
    
        if ( sico.ZONES != 1 ) pzones[i] = 0;
        if ( pzones[i] > nzones ) nzones = pzones[i];
    }
    
    nzones += 1;
    
    float vol_zone1[nzones], vol_zone2[nzones], vol_zone3[nzones], vol_czone1[nzones], vol_czone2[nzones], vol_czone3[nzones];


// -- STOP --- Definition of computational domains --------------------------------------------------------------


// -- START -- Computing slopes and topography-following cell sizes ---------------------------------------------


    // *** Slopes are computed from the 3x3 environment of each raster cell
    // *** Cell sizes applied in the simulation depend on the definition of the height-to-depth conversion

    for ( i=0; i<sico.IMAX; i++ ) { // loop over all cells:

        if ( cdomain[i] != 0 ) { // if cell is no 1st row edge cell:

            betax[i] = fbeta( pelev[in[i][5]], pelev[in[i][2]], 2.0, sico ); // slopes
            betay[i] = fbeta( pelev[in[i][4]], pelev[in[i][1]], 2.0, sico );
            betaxy[i] = fbetaxy( betax[i], betay[i] );

        } else { betax[i] = 0; betay[i] = 0; betaxy[i] = 0; }

        if ( sico.CORRHEIGHT != 0 ) {

            dx[i] = sico.CSZ / cos( betax[i] ); // topography-following cell spacing
            dy[i] = sico.CSZ / cos( betay[i] );

        } else {

            dx[i] = sico.CSZ;
            dy[i] = sico.CSZ;
        }
    }


// -- STOP --- Computing slopes and topography-following cell sizes ---------------------------------------------


// -- START -- Converting heights into depths and deriving release momenta --------------------------------------


    for ( i=0; i<sico.IMAX; i++ ) {

        sico.COLLAPSE = 0;
        if ( sico.TRELSTOP == 1 ) { if ( ptrelstop[i] < 0 && ptrelstop[i] != sico.UNDEF ) { sico.COLLAPSE = 1; ptrelstop[i] *= -1; }} // checking whether to use collapse mode

        if ( sico.MODEL <= 3 ) {

            if ( sico.TRELEASE == 0 || ( sico.TRELSTOP == 1 && sico.COLLAPSE == 0 )) {

                aw[i][0] = fconvin( phrelease[i], betaxy[i], sico ); // release depth
                aw[i][1] = aw[i][0] * pviny[i] * -1; // release momenta in x and y directions
                aw[i][2] = aw[i][0] * pvinx[i];
                
            } else { aw[i][0] = 0; aw[i][1] = 0; aw[i][2] = 0; }
            
        } else if ( sico.MODEL == 7 ) {

            if ( sico.TRELEASE == 0 || ( sico.TRELSTOP == 1 && sico.COLLAPSE == 0 )) {

                aw[i][0] = fconvin( phrelease[i], betaxy[i], sico ); // PHASE 1 release depth
                aw[i][3] = fconvin( phrelease2[i], betaxy[i], sico ); // PHASE 2 release depth
                aw[i][6] = fconvin( phrelease3[i], betaxy[i], sico ); // PHASE 3 release depth

                aw[i][1] = aw[i][0] * pviny[i] * -1; // PHASE 1 release momenta in x and y directions
                aw[i][2] = aw[i][0] * pvinx[i];
                aw[i][4] = aw[i][3] * pviny2[i] * -1; // PHASE 2 release momenta in x and y directions
                aw[i][5] = aw[i][3] * pvinx2[i];
                aw[i][7] = aw[i][6] * pviny3[i] * -1; // PHASE 3 release momenta in x and y directions
                aw[i][8] = aw[i][6] * pvinx3[i];

            } else { aw[i][0] = 0; aw[i][1] = 0; aw[i][2] = 0; aw[i][3] = 0; aw[i][4] = 0; aw[i][5] = 0; aw[i][6] = 0; aw[i][7] = 0; aw[i][8] = 0; }
        }
    }


// -- STOP --- Converting heights into depths and deriving release momenta --------------------------------------


// -- START -- Initializing vector elements, volumes, and maximum depths ----------------------------------------


    hflow_max = 0; // initializing maximum mixture or PHASE 1 flow depth
    hflow_maxmax = 0; // initializing absolute maximum flow depth
    vol_flow = 0; // initializing mixture or PHASE 1 flow volume
    vol_entr = 0; // initializing mixture or PHASE 1 entrained or deposited volume
    
    for ( z=0; z<nzones; z++ ) vol_zone1[z] = 0; // zone-specific volumes

    if ( sico.MODEL == 7 ) {
    
        hflow_max2 = 0; hflow_max3 = 0; // initializing maximum PHASE 2 and PHASE 3 flow depths
        vol_flow2 = 0; vol_flow3 = 0;  // initializing PHASE 2 and PHASE 3 flow volumes
        vol_entr2 = 0; vol_entr3 = 0; // initializing PHASE 2 and PHASE 3 entrained or deposited volumes
        
        for ( z=0; z<nzones; z++ ) { vol_zone2[z] = 0; vol_zone3[z] = 0; } // zone-specific volumes
    }

    for ( i=0; i<sico.IMAX; i++ ) {

        if ( sico.CORRHEIGHT == 0 ) carea = pow ( sico.CSZ, 2 );
        else carea = pow ( sico.CSZ , 2 ) * pow( 1 - pow( sin( betax[i] ) , 2 ) * pow ( sin( betay[i] ) , 2 ) , 0.5 ) / ( cos( betax[i] ) * cos( betay[i] ) );
            // topography-following area of cell

        if ( aw[i][0] > hflow_max ) hflow_max = aw[i][0]; // updating maximum mixture or PHASE 1 release depth
        vol_flow += aw[i][0] * carea; // updating mixture or PHASE 1 release volume
        vol_zone1[pzones[i]] += aw[i][0] * carea;

        if ( sico.MODEL <= 3 ) { // one-phase models

            for ( k=3; k<7; k++ ) aw[i][k] = 0;
            if ( sico.RELM == 1 ) aw[i][7] = aw[i][0]; else aw[i][7] = 0;
            for ( k=8; k<11; k++ ) aw[i][k] = 0;
            
            if ( hflow_max > hflow_maxmax ) hflow_maxmax = hflow_max; // updating absolute maximum flow depth
            
        }  else if ( sico.MODEL == 7 ) { // multi-phase model

            for ( k=9; k<14; k++ ) aw[i][k] = 0;
            aw[i][15] = aw[i][0] + aw[i][3] + aw[i][6];
            for ( k=16; k<24; k++ ) aw[i][k] = 0;
            if ( sico.RELM == 1 ) aw[i][25] = aw[i][0]; else aw[i][25] = 0;
            aw[i][26] = 0;
            if ( sico.RELM2 == 1 ) aw[i][27] = aw[i][3]; else aw[i][27] = 0;
            aw[i][28] = 0;
            if ( sico.RELM3 == 1 ) aw[i][29] = aw[i][6]; else aw[i][29] = 0;
            aw[i][30] = 0;
            aw[i][31] = aw[i][0] + aw[i][3] + aw [i][6];
            for ( k=32; k<40; k++ ) aw[i][k] = 0;

            if ( aw[i][3] > hflow_max2 ) hflow_max2 = aw[i][3]; // updating maximum PHASE 2 release depth
            if ( aw[i][6] > hflow_max3 ) hflow_max3 = aw[i][6]; // updating maximum PHASE 3 release depth
            if ( hflow_max + hflow_max2 + hflow_max3 > hflow_maxmax ) hflow_maxmax = hflow_max + hflow_max2 + hflow_max3; // updating absolute maximum flow depth
            
            vol_flow2 += aw[i][3] * carea; // updating PHASE 2 release volume
            vol_flow3 += aw[i][6] * carea; // updating PHASE 3 release volume
            vol_zone2[pzones[i]] += aw[i][3] * carea;
            vol_zone3[pzones[i]] += aw[i][6] * carea;
        }

        for ( l=-8; l<=-1; l++ ) aw[i][nvect_all+l] = sico.UNDEF; // velocities and phase fractions

        if ( sico.CURVCTRL > 2 ) { for ( j=1; j<9; j++ ) { cready[i][j] = 0; if ( sico.MODEL == 7 ) { cready2[i][j] = 0; cready3[i][j] = 0;  }}}
            // degree of fill of cell with regard to flux in a given direction (ratio)
        cstopped[i] = 0; // control for stopping
    }

    for ( ix=0; ix<ib[0]; ix++ ) {

        i = ibasket[0][ix];
        if ( cdomain[i] != 0 ) {

            asigma_xelev[i] = fsigma( pelev[i], pelev[in[i][5]], pelev[in[i][2]], 1, dx, dy, i, sico ); // gradients of elevation
            asigma_yelev[i] = fsigma( pelev[i], pelev[in[i][4]], pelev[in[i][1]], 2, dx, dy, i, sico );
        }
    }

    for ( ix=0; ix<ib[0]; ix++ ) {

        i = ibasket[0][ix];
        if ( cdomain[i] != 0 ) {

            fj[0][0]=1; fj[0][1]=1; fj[1][0]=1; fj[1][1]=-1; fj[2][0]=-1; fj[2][1]=1; fj[3][0]=-1; fj[3][1]=-1; // factors for positive or negative gradients
            for ( j=0; j<4; j++ ) wintelev[i][j]=pelev[in[i][j]]+fj[j][0]*0.25*dx[i]*asigma_xelev[in[i][j]]+fj[j][1]*0.25*dy[i]*asigma_yelev[in[i][j]];
        }
    }


// -- STOP --- Initializing vector elements, volumes, and maximum depths ----------------------------------------


// -- START -- Definition of hydrograph profiles ----------------------------------------------------------------


    // *** Profile parameters needed for adding material through input hydrographs and for quantifying material at output hydrographs at each time step (if applicable)
    // *** Centre and end points of each hydrograph profile are written to text files

    if ( hydrograph == 1 ) { // if at least one hydrograph is provided
    
        for ( i=0; i<sico.IMAX; i++ ) { 

            for ( hydj = 0; hydj < ( hydnin + hydnout ); hydj++ ) {
                if ( px[i] == hydx[hydj] && py[i] == hydy[hydj] ) {

                    hydi[hydj] = i; // coordinates and elevation of centres of hydrographs
                    hydelev[hydj] = pelev[i];
                }
            }
            qtinit[hydj] = sico.UNDEF; // initializing start time of hydrograph discharge
        }

        hydp = alloc_imatrix( (int)(2*hydlmax/sico.CSZ+2), hydnin+hydnout ); // allocating memory to hydrograph profiles

        for ( hydj = 0; hydj < ( hydnin + hydnout ); hydj++ ) { // loop over all hydrographs

            if ( hydalpha[hydj] < 0 ) hydalpha[hydj] = falpha ( betax[hydi[hydj]], betay[hydi[hydj]], sico ) + sico.PI * 0.5;
                // cross-slope direction
            if ( hydalpha[hydj] > 2 * sico.PI ) hydalpha[hydj] -= 2 * sico.PI;
            hydp0 = fhydp ( hydx[hydj], hydy[hydj], fabs(hydl[hydj]), px, py, hydalpha[hydj], hydlmax, sico );
            for ( i=0; i<(int)(2*hydlmax/sico.CSZ+2); i++ ) hydp[i][hydj] = hydp0[i]; // hydrograph profile
            free( hydp0 );

            if ( hydl[hydj] > 0 ) { // Adjusting central point according to lowest elevation along hydrograph profile

                pelevhydtest = sico.UNDEF * -1;

                for ( hydk=1; hydk<= hydp[0][hydj]; hydk++ ) { // loop over all cells of hydrograph profile

                    if ( pelev[hydp[hydk][hydj]] < pelevhydtest ) {

                        pelevhydtest = pelev[hydp[hydk][hydj]];
                        hydi[hydj] = hydp[hydk][hydj];
                        hydx[hydj] = px[hydp[hydk][hydj]];
                        hydy[hydj] = py[hydp[hydk][hydj]];
                    }
                }

                hydelev[hydj] = pelevhydtest;
            }

            hydx_metric[hydj] = ( float ) hydy[hydj] * sico.CSZ + sico.BDWEST;
            hydy_metric[hydj] = sico.BDNORTH - ( float ) hydx[hydj] * sico.CSZ;

            hydout_xmin_metric = ( float ) py[hydp[1][hydj]] * sico.CSZ + sico.BDWEST; // external coordinates of profile terminal points
            hydout_xmax_metric = ( float ) py[hydp[hydp[0][hydj]][hydj]] * sico.CSZ + sico.BDWEST;
            hydout_ymin_metric = sico.BDNORTH - ( float ) px[hydp[1][hydj]] * sico.CSZ;
            hydout_ymax_metric = sico.BDNORTH - ( float ) px[hydp[hydp[0][hydj]][hydj]] * sico.CSZ;

            if ( hydj < hydnin ) fprintf(f_hydout, "I%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", hydj+1,
                hydout_xmin_metric, hydx_metric[hydj], hydout_xmax_metric, hydout_ymin_metric, hydy_metric[hydj], hydout_ymax_metric);
            else fprintf(f_hydout, "O%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", hydj-hydnin+1,
                hydout_xmin_metric, hydx_metric[hydj], hydout_xmax_metric, hydout_ymin_metric, hydy_metric[hydj], hydout_ymax_metric);
                // writing hydrograph profile terminal points and centre to text file

        }
    }


// -- STOP --- Definition of hydrograph profiles ----------------------------------------------------------------


// -- START -- Definition / initialization of control variables -------------------------------------------------


    tlength0 = sico.CFL[1]; // length of first time step
    tlengthpre = sico.CFL[1]; // length of previous time step
    tsum = 0; // total time
    tint = 0; // time since last output
    cflmax = 0; // maximum CFL value
    nsum = 1; // number of time steps
    nout = 1; // number of output time steps
    ccontinue = 1; // control whether simulation should be continued (1 or 0)
    csuccess = 1; // control for success of simulation (1 or 0)
    //coverflow = 0; // control for overflow i.e. part of the flow leaving the study area (1 or 0), deactivated
    if ( sico.MODEL <= 3 ) qh_test = hflow_max0; else qh_test = hflow_max0 + hflow_max02 + hflow_max03;  // control release height for collapse
    ctrl_noosc = 0; // control for effective application of surface control
    vmax = 0; // maximum flow velocity
    cflowpre = 0; // control for activation of flow model at previous time step
    vol_hyd = 0; vol_hyd2 = 0; vol_hyd3 = 0; ctrl_hydout = 1; // hydrograph volumes and control for hydrograph output
    hekin_max = 0; // initialization of maximum value of flow kinetic energy


// -- STOP --- Definition / initialization of control variables -------------------------------------------------


// -- START -- Summary of release and initial state of the flow -------------------------------------------------


    // *** Information is written to the display (immediately) and to the summary file

    if ( sico.MODEL <= 3 ) { // one-phase models

        if ( hflow_max0 >= 1000 ) prec_hflow = 0; // precision for summary of maximum flow heights and depths
        else if ( hflow_max0 >= 100 ) prec_hflow = 1;
        else if ( hflow_max0 >= 10 ) prec_hflow = 2;
        else if ( hflow_max0 >= 1 ) prec_hflow = 3;
        else if ( hflow_max0 >= 0.1 ) prec_hflow = 4;
        else if ( hflow_max0 == 0 ) prec_hflow = 2;
        else prec_hflow = 5;

        if ( vol_flow0 >= 100000000 ) prec_vol = 0; // precision for summary of flow volumes
        else if ( vol_flow0 >= 1000000 ) prec_vol = 1;
        else if ( vol_flow0 >= 100000 ) prec_vol = 2;
        else if ( vol_flow0 >= 10000 ) prec_vol = 3;
        else if ( vol_flow0 >= 1000 ) prec_vol = 4;
        else if ( vol_flow0 == 0 ) prec_vol = 3;
        else prec_vol = 5;

        prec_ekin = prec_vol-3; if ( prec_ekin < 0 ) prec_ekin = 0; // precision for summary of flow kinetic energies

        printf("   nout\tnsum\tcfl\ttlength\ttsum\tdmax\tvmax\tvolume\tekin\n");
        printf("   R\t-----\t-----\t-----\t-----\t*%.*f\t-----\t*%.*f\t-----\n",
            prec_hflow, hflow_max0, prec_vol, vol_flow0/1000 ); // displaying release parameters


        #ifdef WITHGRASS


            printf("   0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t-----\n",
                prec_hflow, hflow_max, prec_vol, vol_flow/1000 ); // displaying initial parameters


        #else


            printf("   0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t-----",
                prec_hflow, hflow_max, prec_vol, vol_flow/1000 ); // displaying initial parameters


        #endif


        fflush(stdout); // forcing immediate display

        fprintf(f_summary, "nout\tnsum\tcfl\ttlength\ttsum\tdmax\tvmax\tvolume\tekin\n");
        fprintf(f_summary, "R\t-----\t-----\t-----\t-----\t*%.*f\t-----\t*%.*f\t-----\n", prec_hflow, hflow_max0, prec_vol, vol_flow0 ); // writing release parameters to summary file
        fprintf(f_summary, "0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t\n", prec_hflow, hflow_max, prec_vol, vol_flow ); // writing initial parameters to summary file
        
    } else if ( sico.MODEL == 7 ) { // multi-phase model

        if ( hflow_max0 >= 1000 || hflow_max02 >= 1000 || hflow_max03 >= 1000 ) prec_hflow = 0;  // precision for summary of maximum flow heights and depths
        else if ( hflow_max0 >= 100 || hflow_max02 >= 100 || hflow_max03 >= 100 ) prec_hflow = 1;
        else if ( hflow_max0 >= 10 || hflow_max02 >= 10 || hflow_max03 >= 10 ) prec_hflow = 2;
        else if ( hflow_max0 >= 1 || hflow_max02 >= 1 || hflow_max03 >= 1 ) prec_hflow = 3;
        else if ( hflow_max0 >= 0.1 || hflow_max02 >= 0.1 || hflow_max03 >= 0.1 ) prec_hflow = 4;
        else if ( hflow_max0 == 0 && hflow_max02 == 0 && hflow_max03 == 0 ) prec_hflow = 2;
        else prec_hflow = 5;

        if ( vol_flow0 >= 100000000 || vol_flow02 >= 10000000 || vol_flow03 >= 10000000 ) prec_vol = 0; // precision for summary of flow volumes
        else if ( vol_flow0 >= 1000000 || vol_flow02 >= 1000000 || vol_flow03 >= 1000000 ) prec_vol = 1;
        else if ( vol_flow0 >= 100000 || vol_flow02 >= 100000 || vol_flow03 >= 100000 ) prec_vol = 2;
        else if ( vol_flow0 >= 10000 || vol_flow02 >= 10000 || vol_flow03 >= 10000 ) prec_vol = 3;
        else if ( vol_flow0 >= 1000 || vol_flow02 >= 1000  || vol_flow03 >= 1000 ) prec_vol = 4;
        else if ( vol_flow0 == 0 && vol_flow02 == 0 && vol_flow03 == 0 ) prec_vol = 3;
        else prec_vol = 5;

        prec_ekin = prec_vol-3; if ( prec_ekin < 0 ) prec_ekin = 0; // precision for summary of flow kinetic energies

        printf("   nout\tnsum\tcfl\ttlength\ttsum\tdmax1\tvmax1\tdmax2\tvmax2\tdmax3\tvmax3\tvolume1\tvolume2\tvolume3\tekin\n");
        printf("   R\t-----\t-----\t-----\t-----\t*%.*f\t-----\t*%.*f\t-----\t*%.*f\t-----\t*%.*f\t*%.*f\t*%.*f\t-----\n",
            prec_hflow, hflow_max0, prec_hflow, hflow_max02, prec_hflow, hflow_max03, prec_vol, vol_flow0/1000, prec_vol, vol_flow02/1000, prec_vol, vol_flow03/1000); // displaying release parameters


        #ifdef WITHGRASS


            printf("   0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t%.*f\t%.*f\t-----\n", prec_hflow, hflow_max, prec_hflow, hflow_max2,
                prec_hflow, hflow_max3, prec_vol, vol_flow/1000, prec_vol, vol_flow2/1000, prec_vol, vol_flow3/1000); // displaying initial parameters


        #else


            printf("   0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t%.*f\t%.*f\t-----", prec_hflow, hflow_max, prec_hflow, hflow_max2,
                prec_hflow, hflow_max3, prec_vol, vol_flow/1000, prec_vol, vol_flow2/1000, prec_vol, vol_flow3/1000); // displaying initial parameters


        #endif


        fflush(stdout); // forcing immediate display

        fprintf(f_summary, "nout\tnsum\tcfl\ttlength\ttsum\tdmax1\tvmax1\tdmax2\tvmax2\tdmax3\tvmax3\tvolume1\tvolume2\tvolume3\tekin\n");
        fprintf(f_summary, "R\t-----\t-----\t-----\t-----\t*%.*f\t-----\t*%.*f\t-----\t*%.*f\t-----\t*%.*f\t*%.*f\t*%.*f\t-----\n",
            prec_hflow, hflow_max0, prec_hflow, hflow_max02, prec_hflow, hflow_max03, prec_vol, vol_flow0, prec_vol, vol_flow02, prec_vol, vol_flow03 ); // writing release parameters to summary file
        fprintf(f_summary, "0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t%.*f\t%.*f\t-----\n",
            prec_hflow, hflow_max, prec_hflow, hflow_max2, prec_hflow, hflow_max3, prec_vol, vol_flow, prec_vol, vol_flow2, prec_vol, vol_flow3 ); // writing initial parameters to summary file
    }
    
    fprintf(f_volumes, "nout\tnsum"); // volumes file
    
    for ( z=0; z<nzones; z++ ) {
    
        fprintf(f_volumes, "\tz%i_vflow1\tz%ivchange1", z, z );
        if ( sico.MODEL == 7 ) fprintf(f_volumes, "\tz%i_vflow2\tz%ivchange2\tz%i_vflow3\tz%ivchange3", z, z, z, z );
    }

    fprintf(f_volumes, "\n0\t0"); // volumes file
    
    for ( z=0; z<nzones; z++ ) {
    
        fprintf(f_volumes, "\t%.3f\t0", vol_zone1[z] );
        if ( sico.MODEL == 7 ) fprintf(f_volumes, "\t%.3f\t0\t%.3f\t0", vol_zone2[z], vol_zone3[z] );
    }
    
    fprintf(f_volumes, "\n");
   
    
// -- STOP --- Summary of release and initial state of the flow -------------------------------------------------

    
// *** Start of loop over time steps until break criterion is met: ----------------------------------------------


    while ( (int)(round( 1000 * tsum )) <= (int)(round( 1000 * tmax )) && ccontinue == 1 ) {
    
        if ( (int)(round( 1000 * tint )) >= (int)( round( 1000 * tout ))) tint -= tout; // resetting time interval for output, if required


// *** Start of loop over two steps, each moving the system half of a cell (NOC scheme) -------------------------


        for ( p=0; p<2; p++ ) {
        
            cslide = 0; cflow = 0;
            for ( i=0; i<sico.IMAX; i++ ) {
                for ( k=0; k<sico.NVECTMIN; k++ ) {
                
                    awt[i][k] = 0; // resetting temporary state variables
                    wintd[i][k] = 0;
                    wintdtest[i][k] = 0;
                }
                
                pxslide[i] = 0; anctr = 0;

                if ( ptslide[i] > tsum || ptslide[i] == -777 ) cslide = 1; // resetting controls for relevance of sliding and flowing
                if ( ptslide[i] <= tsum || ptslide[i] == -777 ) cflow = 1;

                if ( ptslide[i] == -777 && sico.MODEL == 7 && cdomain[i] != 0 ) {
                
                    anctr = 1;
                    for ( j=0; j<9; j++ ) { if ( aw[in[i][j]][6] > sico.HFLOWMIN ) anctr = 0; }
                }
                
                if ( ptslide[i] > tsum || anctr == 1 ) pxslide[i] = 1;
            }

            if ( cflow == cflowpre ) tlength = tlength0; // time step length
            else tlength = 0.001;
            cflowpre = cflow;
            ib[1] = 0;


// -- START -- Preparing data for progressive collapse (constant volume) ----------------------------------------


            // *** Options trelease and trelstop, negative value of trelstop activates this release mode
            // *** Release of mass starts from the top of the release area
            // *** Similar volumes are released at each time step, until all release material has been released
            // *** The basal surface is updated at each time step (in contrast to all other release types, the initial terrain includes the release mass)
            // *** Applicable with the one-phase models and the multi-phase model

            if ( sico.COLLAPSE == 1 && sico.TRELEASE == 1 && sico.TRELSTOP == 1 ) {

                qtrelspan = 0;

                for ( i=0; i<sico.IMAX; i++ ) {

                    if ( tsum >= ptrelease[i] && tsum < ptrelstop[i] ) {

                        if ( ptrelstop[i] - ptrelease[i] > qtrelspan ) {

                            qtrelspan = ptrelstop[i] - ptrelease[i]; // start, stop, and duration of progressive collapse
                            qtrelstart = ptrelease[i];
                            qtrelstop = ptrelstop[i];
                        }
                    }
                }

                rhrem = tlengthpre / qtrelspan; // fraction of input hydrograph to be released at current time step

                if ( tsum >= qtrelstart && tsum < qtrelstop && qh_test > 0 ) {

                    qvol_test = 0; // initializing target release volume

                    if ( sico.MODEL <= 3 ) vol_flow0all = vol_flow0;
                    else vol_flow0all = vol_flow0 + vol_flow02 + vol_flow03;

                    while ( qvol_test < vol_flow0all * rhrem ) { // loop over test release heights until target release volume is met

                        qvol_test = 0;
                        qh_test -= vol_flow0all * pow(10, -8 ); // updating test release height
                        for ( i=0; i<sico.IMAX; i++ )  {

                            if ( tsum >= ptrelease[i] && tsum < ptrelstop[i] ) {
                            
                                if ( sico.MODEL <= 3 ) phreleaseall = phrelease[i]; // total release height
                                else phreleaseall = phrelease[i] + phrelease2[i] + phrelease3[i];
                            
                                qvol_test += sico.CSZ * sico.CSZ * ffmax( 0, phreleaseall - qh_test ); // updating target volume
                            }
                        }
                    }

                    for ( i=0; i<sico.IMAX; i++ ) {

                        if ( tsum >= ptrelease[i] && tsum < ptrelstop[i] && qh_test > 0 ) {

                            if ( sico.MODEL <= 3 ) { phreleaseall = phrelease[i]; alpha1 = 1; } // total release height and phase fractions
                            else { phreleaseall = phrelease[i] + phrelease2[i] + phrelease3[i];
                                if ( phreleaseall > 0 ) { alpha1 = phrelease[i] / phreleaseall; alpha2 = phrelease2[i] / phreleaseall; }
                            }

                            if ( phreleaseall > 0 ) {

                                qhreleaseall = ffmax( 0, phreleaseall - qh_test ); // release hydrograph height
                                qhrelease[i] = qhreleaseall * alpha1; phrelease[i] -= qhrelease[i]; // phase-specific release hydrograph heights and remaining release heights
                                
                                if ( sico.MODEL == 7 ) { 

                                    qhrelease2[i] = qhreleaseall * alpha2; phrelease2[i] -= qhrelease2[i];
                                    qhrelease3[i] = qhreleaseall * ( 1 - alpha1 - alpha2 ); phrelease3[i] -= qhrelease3[i];
                                }
                                
                                pelev[i] -= qhreleaseall; // updating terrain elevation

                                if ( sico.MODEL <= 3 ) aw[i][3] = aw[i][3] - qhrelease[i]; // updating changes of basal topography
                                else { aw[i][9] = aw[i][9] - qhrelease[i]; aw[i][10] = aw[i][10] - qhrelease2[i]; aw[i][11] = aw[i][11] - qhrelease3[i]; }
                            }
                        }
                    }

                } else {

                    for ( i=0; i<sico.IMAX; i++ ) {

                        qhrelease[i] = 0; // setting release hydrograph height to zero outside of release time span
                        if ( sico.MODEL == 7 ) { qhrelease2[i] = 0; qhrelease3[i] = 0; }
                    }
                }
                
                for ( ix=0; ix<ib[0]; ix++ ) { // updating slopes and topography-following cell sizes

                    i = ibasket[0][ix];

                    betax[i] = fbeta( pelev[in[i][5]], pelev[in[i][2]], 2.0, sico );
                    betay[i] = fbeta( pelev[in[i][4]], pelev[in[i][1]], 2.0, sico );
                    betaxy[i] = fbetaxy( betax[i], betay[i] );

                    if ( sico.CORRHEIGHT != 0 ) {

                        dx[i] = sico.CSZ / cos( betax[i] ); // topography-following cell spacing
                        dy[i] = sico.CSZ / cos( betay[i] );
                    }
                }

            } else if ( sico.TRELEASE == 1 ) {

                for ( i=0; i<sico.IMAX; i++ ) {

                    qhrelease[i] = phrelease[i]; // with release time but no progressive collapse, setting release hydrograph height to release height for further use
                    if ( sico.MODEL == 7 ) { qhrelease2[i] = phrelease2[i]; qhrelease3[i] = phrelease3[i]; }
                }
            }
            
            
// -- STOP --- Preparing data for progressive collapse (constant volume) ----------------------------------------


// -- START -- Preparing data for extrusive release -------------------------------------------------------------


            // *** Options trelease and trelstop, positive value of trelstop activates this release mode
            // *** Release height has to be provided as rate of extrusion (m/s)
            // *** Similar volumes are released at each time step, until all release material has been released
            // *** Applicable with the one-phase models and the multi-phase model

            if ( sico.TRELSTOP == 1 && sico.COLLAPSE == 0 ) tfact = tlengthpre; else tfact = 1;

            for ( i=0; i<sico.IMAX; i++ ) {

                if ( sico.TRELEASE == 1 ) ttrelease = ptrelease[i]; // release time
                else ttrelease = 0;

                if ( sico.TRELSTOP == 1 ) ttrelstop = ptrelstop[i]; // release hydrograph stop time
                else ttrelstop = ttrelease;

                if ( sico.TRELEASE == 1 && tsum >= ttrelease && ( ptrelease[i] != 1000 * sico.UNDEF || tsum < ttrelstop )) {

                    if ( sico.MODEL <= 3 ) {

                        aw[i][0] += fconvin( qhrelease[i], betaxy[i], sico ) * tfact; // release depth

                        aw[i][1] += fconvin( qhrelease[i], betaxy[i], sico ) * tfact * pviny[i] * -1;
                            // release momentum in x direction
                        aw[i][2] += fconvin( qhrelease[i], betaxy[i], sico ) * tfact * pvinx[i];
                            // release momentum in y direction
                            
                    } else if ( sico.MODEL == 7 ) {

                        aw[i][0] += fconvin( qhrelease[i], betaxy[i], sico ) * tfact; // PHASE 1 release depth
                        aw[i][3] += fconvin( qhrelease2[i], betaxy[i], sico ) * tfact; // PHASE 2 release depth
                        aw[i][6] += fconvin( qhrelease3[i], betaxy[i], sico ) * tfact; // PHASE 3 release depth

                        aw[i][1] += fconvin( qhrelease[i], betaxy[i], sico ) * tfact * pviny[i] * -1; // PHASE 1 release momenta in x and y directions
                        aw[i][2] += fconvin( qhrelease[i], betaxy[i], sico ) * tfact * pvinx[i];

                        aw[i][4] += fconvin( qhrelease2[i], betaxy[i], sico ) * tfact * pviny2[i] * -1; // PHASE 2 release momenta in x and y directions
                        aw[i][5] += fconvin( qhrelease2[i], betaxy[i], sico ) * tfact * pvinx2[i];

                        aw[i][7] += fconvin( qhrelease3[i], betaxy[i], sico ) * tfact * pviny3[i] * -1; // PHASE 3 release momenta in x and y directions
                        aw[i][8] += fconvin( qhrelease3[i], betaxy[i], sico ) * tfact * pvinx3[i];
                    }

                    if ( sico.COLLAPSE == 0 ) ptrelease[i] = 1000 * sico.UNDEF;
                }
            }


// -- STOP --- Preparing data for extrusive release -------------------------------------------------------------


// -- START -- Updating flow depths and velocities according to input hydrographs -------------------------------


            vol_hydbef = 0; vol_hydaft = 0;
            if ( sico.MODEL == 7 ) { vol_hydbef2 = 0; vol_hydaft2 = 0; }
            if ( sico.MODEL == 7 ) { vol_hydbef3 = 0; vol_hydaft3 = 0; }

            if ( hydrograph == 1 ) {

                if ( tsum == 0 || sico.HYDADD == 0 ) tlengthx = 1;
                else tlengthx = tlengthpre;

                for ( hydj = 0; hydj < hydnin; hydj++ ) { // loop over all input hydrographs

                    for ( hydk=0; hydk<=hydtmax[hydj]; hydk++ ) { // identifying relevant line of input hydrograph
                        if ( hydhyd[hydk][0][hydj] <= tsum ) hydt=hydk;
                        else break;
                    }

                    if ( sico.MODEL <= 3 ) hhyd0 = hydhyd[hydt][1][hydj]; // total flow height at centre of hydrograph
                    else if ( sico.MODEL == 7 ) hhyd0 = hydhyd[hydt][1][hydj] + hydhyd[hydt][3][hydj] + hydhyd[hydt][5][hydj];

                    if ( hhyd0 > sico.HFLOWMIN && sico.HYDADD < 2 ) {

                        for ( hydk=1; hydk<= hydp[0][hydj]; hydk++ ) { // loop over all cells of hydrograph profile

                            if ( sico.CORRHEIGHT == 0 ) carea = pow ( sico.CSZ, 2 );
                            else carea = pow ( sico.CSZ, 2 ) * pow( 1 - pow( sin( betax[hydp[hydk][hydj]] ) , 2 )
                                * pow ( sin( betay[hydp[hydk][hydj]] ) , 2 ) , 0.5 )
                                / ( cos( betax[hydp[hydk][hydj]] ) * cos( betay[hydp[hydk][hydj]] )); // topography-following area of cell

                            hhyd = hhyd0 - ( pelev[hydp[hydk][hydj]] - hydelev[hydj] ); // total flow height at cell

                            for ( k=0; k<sico.NVECTMIN; k++ ) {
                                if ( sico.HYDADD == 1 ) hydbef[k] = aw[hydp[hydk][hydj]][k]; else hydbef[k] = 0;
                            }

                            vol_hydbef += aw[hydp[hydk][hydj]][0] * carea;
                            if ( sico.MODEL == 7 ) vol_hydbef2 += aw[hydp[hydk][hydj]][3] * carea;
                            if ( sico.MODEL == 7 ) vol_hydbef3 += aw[hydp[hydk][hydj]][6] * carea;

                            if ( hhyd > sico.HFLOWMIN ) {

                                aw[hydp[hydk][hydj]][0] = hydbef[0] + tlengthx * hydhyd[hydt][1][hydj] * hhyd / hhyd0 * cos( betaxy[hydp[hydk][hydj]] ); // mixture or PHASE 1 flow depth
                                aw[hydp[hydk][hydj]][1] = hydbef[1] + tlengthx * hydhyd[hydt][1][hydj] * hhyd / hhyd0 * hydhyd[hydt][2][hydj]
                                    * cos( hydalpha[hydj] - sico.PI * 0.5 )  * cos( betaxy[hydp[hydk][hydj]] ); // mixture or PHASE 1 x momentum
                                aw[hydp[hydk][hydj]][2] = hydbef[2] + tlengthx * hydhyd[hydt][1][hydj] * hhyd / hhyd0 * hydhyd[hydt][2][hydj]
                                    * sin( hydalpha[hydj] - sico.PI * 0.5 )  * cos( betaxy[hydp[hydk][hydj]] ); // mixture or PHASE 1 y momentum
                            }

                            if ( hhyd > sico.HFLOWMIN && sico.MODEL == 7 ) {

                                aw[hydp[hydk][hydj]][3] = hydbef[3] + tlengthx * hydhyd[hydt][3][hydj] * hhyd / hhyd0  * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 2 flow depth
                                aw[hydp[hydk][hydj]][4] = hydbef[4] + tlengthx * hydhyd[hydt][3][hydj] * hhyd / hhyd0 * hydhyd[hydt][4][hydj]
                                    * cos( hydalpha[hydj] - sico.PI * 0.5 ) * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 2 x momentum
                                aw[hydp[hydk][hydj]][5] = hydbef[5] + tlengthx * hydhyd[hydt][3][hydj] * hhyd / hhyd0 * hydhyd[hydt][4][hydj]
                                    * sin( hydalpha[hydj] - sico.PI * 0.5 ) * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 2 y momentum
                            }

                            if ( hhyd > sico.HFLOWMIN && sico.MODEL == 7 ) {

                                aw[hydp[hydk][hydj]][6] = hydbef[6] + tlengthx * hydhyd[hydt][5][hydj] * hhyd / hhyd0  * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 3 flow depth
                                aw[hydp[hydk][hydj]][7] = hydbef[7] + tlengthx * hydhyd[hydt][5][hydj] * hhyd / hhyd0 * hydhyd[hydt][6][hydj]
                                    * cos( hydalpha[hydj] - sico.PI * 0.5 ) * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 3 x momentum
                                aw[hydp[hydk][hydj]][8] = hydbef[8] + tlengthx * hydhyd[hydt][5][hydj] * hhyd / hhyd0 * hydhyd[hydt][6][hydj]
                                    * sin( hydalpha[hydj] - sico.PI * 0.5 ) * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 3 y momentum
                            }

                            vol_hydaft += aw[hydp[hydk][hydj]][0] * carea;
                            if ( sico.MODEL == 7 ) vol_hydaft2 += aw[hydp[hydk][hydj]][3] * carea;
                            if ( sico.MODEL == 7 ) vol_hydaft3 += aw[hydp[hydk][hydj]][6] * carea;

                            if ( fvalid( aw[hydp[hydk][hydj]][0], aw[hydp[hydk][hydj]][3], aw[hydp[hydk][hydj]][6], sico ) == 1
                                && cstopped[hydp[hydk][hydj]] == 0 && cdomain[hydp[hydk][hydj]] > 0 ) { // flow cells

                                for ( j=0; j<9; j++ ) { // 1st row flow boundary cells

                                    if ( icheck[in[hydp[hydk][hydj]][j]][1] == 0 && cstopped[in[hydp[hydk][hydj]][j]] == 0
                                        && cdomain[in[hydp[hydk][hydj]][j]] > 0 ) {

                                        ibasket[1][ib[1]] = in[hydp[hydk][hydj]][j];
                                        ib[1] += 1;
                                        icheck[in[hydp[hydk][hydj]][j]][1] = 1;
                                    }
                                }
                            }
                        }

                    } else if ( hhyd0 / ( sico.CSZ * sico.CSZ ) > sico.HFLOWMIN && sico.HYDADD == 2 ) {

                        if ( sico.CORRHEIGHT == 0 ) carea = pow ( sico.CSZ, 2 );
                        else carea = pow ( sico.CSZ, 2 ) * pow( 1 - pow( sin( betax[hydi[hydj]] ) , 2 )
                            * pow ( sin( betay[hydi[hydj]] ) , 2 ) , 0.5 )
                            / ( cos( betax[hydi[hydj]] ) * cos( betay[hydi[hydj]] )); // topography-following area of cell

                        vol_hydbef += aw[hydi[hydj]][0] * carea;
                        if ( sico.MODEL == 7 ) vol_hydbef2 += aw[hydi[hydj]][3] * carea;
                        if ( sico.MODEL == 7 ) vol_hydbef3 += aw[hydi[hydj]][6] * carea;

                        aw[hydi[hydj]][0] += tlengthx * hydhyd[hydt][1][hydj] / carea;
                            // mixture or PHASE 1 flow depth
                        aw[hydi[hydj]][1] += tlengthx * hydhyd[hydt][1][hydj] / carea * hydhyd[hydt][2][hydj]
                            * cos( hydalpha[hydj] - sico.PI * 0.5 ); // mixture or PHASE 1 x momentum
                        aw[hydi[hydj]][2] += tlengthx * hydhyd[hydt][1][hydj] / carea * hydhyd[hydt][2][hydj]
                            * sin( hydalpha[hydj] - sico.PI * 0.5 ); // mixture or PHASE 1 y momentum

                        if ( sico.MODEL == 7 ) {

                            aw[hydi[hydj]][3] += tlengthx * hydhyd[hydt][3][hydj] / carea;
                                // PHASE 2 flow depth
                            aw[hydi[hydj]][4] += tlengthx * hydhyd[hydt][3][hydj] / carea * hydhyd[hydt][4][hydj]
                                * cos( hydalpha[hydj] - sico.PI * 0.5 ); // PHASE 2 x momentum
                            aw[hydi[hydj]][5] += tlengthx * hydhyd[hydt][3][hydj] / carea * hydhyd[hydt][4][hydj]
                                * sin( hydalpha[hydj] - sico.PI * 0.5 ); // PHASE 2 y momentum
                        }

                        if ( sico.MODEL == 7 ) {

                            aw[hydi[hydj]][6] += tlengthx * hydhyd[hydt][5][hydj] / carea;
                                // PHASE 3 flow depth
                            aw[hydi[hydj]][7] += tlengthx * hydhyd[hydt][5][hydj] / carea * hydhyd[hydt][6][hydj]
                                * cos( hydalpha[hydj] - sico.PI * 0.5 ); // PHASE 3 x momentum
                            aw[hydi[hydj]][8] += tlengthx * hydhyd[hydt][5][hydj] / carea * hydhyd[hydt][6][hydj]
                                * sin( hydalpha[hydj] - sico.PI * 0.5 ); // PHASE 3 y momentum
                        }

                        vol_hydaft += aw[hydi[hydj]][0] * carea;
                        if ( sico.MODEL == 7 ) vol_hydaft2 += aw[hydi[hydj]][3] * carea;
                        if ( sico.MODEL == 7 ) vol_hydaft3 += aw[hydi[hydj]][6] * carea;
                        
                        if ( fvalid( aw[hydi[hydj]][0], aw[hydi[hydj]][3], aw[hydi[hydj]][6], sico ) == 1
                            && cstopped[hydi[hydj]] == 0 && cdomain[hydi[hydj]] > 0 ) { // flow cells

                            for ( j=0; j<9; j++ ) { // 1st row flow boundary cells

                                if ( icheck[in[hydi[hydj]][j]][1] == 0 && cstopped[in[hydi[hydj]][j]] == 0
                                    && cdomain[in[hydi[hydj]][j]] > 0 ) {

                                    ibasket[1][ib[1]] = in[hydi[hydj]][j];
                                    ib[1] += 1;
                                    icheck[in[hydi[hydj]][j]][1] = 1;
                                }
                            }
                        }                        
                    }
                }

                vol_hyd += ( vol_hydaft - vol_hydbef ); // volume added through hydrograph
                if ( sico.MODEL == 7 ) vol_hyd2 += ( vol_hydaft2 - vol_hydbef2 );
                if ( sico.MODEL == 7 ) vol_hyd3 += ( vol_hydaft3 - vol_hydbef3 );
            }


// -- STOP --- Updating flow depths and velocities according to input hydrographs -------------------------------


// -- START -- Writing hydrograph infos to files ----------------------------------------------------------------


            if ( p == 1 && hydrograph == 1 ) {

                for ( hydj = 0; hydj < hydnin; hydj++ ) { // loop over all input hydrographs

                    if ( sico.MULT == 0 && ctrl_hydout == 1 ) { // for first time step after output

                        hydh = aw[hydi[hydj]][0] / cos( betaxy[hydi[hydj]] ); // mixture or PHASE 1 flow height
                        if ( hydh > sico.HFLOWMIN ) hydv = pow( pow( aw[hydi[hydj]][1], 2 ) + pow( aw[hydi[hydj]][2], 2 ), 0.5 ) / aw[hydi[hydj]][0];
                        else hydv= 0; // mixture or PHASE 1 flow velocity

                        if ( sico.MODEL <= 3 ) {

                            hyde = aw[hydi[hydj]][3]; // entrained height
                            hydh2 = 0; hydv2 = 0; hyde2 = 0; hydh3 = 0; hydv3 = 0; hyde3 = 0;
                            
                        } else if ( sico.MODEL == 7 ) {

                            hydh2 = aw[hydi[hydj]][3] / cos( betaxy[hydi[hydj]] ); // PHASE 2 flow height
                            if ( hydh2 > sico.HFLOWMIN ) hydv2 = pow( pow( aw[hydi[hydj]][4], 2 ) + pow( aw[hydi[hydj]][5], 2 ), 0.5 ) / aw[hydi[hydj]][3];
                            else hydv2= 0; // PHASE 2 flow velocity
                            hydh3 = aw[hydi[hydj]][6] / cos( betaxy[hydi[hydj]] ); // PHASE 3 flow height
                            if ( hydh3 > sico.HFLOWMIN ) hydv3 = pow( pow( aw[hydi[hydj]][7], 2 ) + pow( aw[hydi[hydj]][8], 2 ), 0.5 ) / aw[hydi[hydj]][6];
                            else hydv3= 0; // PHASE 3 flow velocity
                            hyde = aw[hydi[hydj]][9];
                            hyde2 = aw[hydi[hydj]][10];
                            hyde3 = aw[hydi[hydj]][11]; // changes of basal topography
                        }

                        hydq = 0; hydq2 = 0; hydq3 = 0; // resetting discharges

                        if ( hydalpha[hydj] == 0 || hydalpha[hydj] == sico.PI * 0.5 || hydalpha[hydj] == sico.PI || hydalpha[hydj] == 3 * sico.PI * 0.5 ) hydfalpha = 1;
                        else if ( fabs( 1 / sin( hydalpha[hydj] )) < fabs( 1 / cos( hydalpha[hydj] ))) hydfalpha = fabs( 1 / sin( hydalpha[hydj] ));
                        else hydfalpha = fabs( 1 / cos( hydalpha[hydj] )); // correction factor for profile direction

                        for ( hydk = 1; hydk <= hydp[0][hydj]; hydk++ ) { // loop over all cells of profile

                            alpha = falpha( betax[hydp[hydk][hydj]], betay[hydp[hydk][hydj]], sico ); // aspect
                            hydbeta = atan ( tan( betaxy[hydp[hydk][hydj]] ) * cos ( alpha - hydalpha[hydj] + sico.PI * 0.5 )); // corrected slope

                            hydfcorr = sico.CSZ * hydfalpha / cos( hydbeta ); // reference length for discharge
                            hydm0 = pow( pow(aw[hydp[hydk][hydj]][1], 2 ) + pow( aw[hydp[hydk][hydj]][2], 2), 0.5 );

                            if ( aw[hydp[hydk][hydj]][0] > sico.HFLOWMIN && hydm0 > 0 ) {

                                hydmx = aw[hydp[hydk][hydj]][1];
                                hydmy = aw[hydp[hydk][hydj]][2];

                                if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 );
                                else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                hydq += hydm * hydfcorr; // updating mixture or PHASE 1 discharge
                            }

                            if ( sico.MODEL == 7 ) {

                                hydm0 = pow( pow(aw[hydp[hydk][hydj]][4], 2 ) + pow( aw[hydp[hydk][hydj]][5], 2), 0.5 );

                                if ( aw[hydp[hydk][hydj]][3] > sico.HFLOWMIN && hydm0 > 0 ) {

                                    hydmx = aw[hydp[hydk][hydj]][4];
                                    hydmy = aw[hydp[hydk][hydj]][5];

                                    if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 );
                                    else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                    hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                    hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                    hydq2 += hydm * hydfcorr; // updating PHASE 2 discharge
                                }

                                hydm0 = pow( pow(aw[hydp[hydk][hydj]][7], 2 ) + pow( aw[hydp[hydk][hydj]][8], 2), 0.5 );

                                if ( aw[hydp[hydk][hydj]][6] > sico.HFLOWMIN && hydm0 > 0 ) {

                                    hydmx = aw[hydp[hydk][hydj]][7];
                                    hydmy = aw[hydp[hydk][hydj]][8];

                                    if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 );
                                    else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                    hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                    hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                    hydq3 += hydm * hydfcorr; // updating PHASE 3 discharge
                                }
                            }
                        }

                        fprintf(f_hydinfo[hydj], "%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                            tsum, hydh, hydv, hyde, hydq, hydh2, hydv2, hyde2, hydq2, hydh3, hydv3, hyde3, hydq3); // writing hydrograph info to file
                    }
                }
                
                ctrl_hydout = 0; // resetting control for hydrograph output
            }


// -- STOP --- Writing hydrograph infos to files ----------------------------------------------------------------


// -- START -- Updating computational domains -------------------------------------------------------------------


            // *** cdomain[i] = 0: 1st row edge cell
            // *** cdomain[i] = 1: no-flux cell (suppression of oscillations)
            // *** cdomain[i] = 3: included cells, flow depth < minimum at & around cell
            // *** cdomain[i] = 4: included cells, flow depth >= minimum around cell (flow boundary)
            // *** cdomain[i] = 5: included cells, flow depth >= minimum (flow)

            if ( cslide == 0 ) iloop = ib[0]; else iloop = sico.IMAX;
            for ( ix=0; ix<iloop; ix++ ) {

                if ( cslide == 0 ) i = ibasket[0][ix]; else i = ix;

                if (( fvalid( aw[i][0], aw[i][3], aw[i][6], sico ) == 1 || ( sico.TRELEASE == 1 && ( ptrelease[i] > 0 || sico.COLLAPSE != 0 )))
                    && cstopped[i] == 0 && cdomain[i] > 0 ) { // flow cells

                    for ( j=0; j<9; j++ ) { // 1st row flow boundary cells

                        if ( icheck[in[i][j]][1] == 0 && cstopped[in[i][j]] == 0 && cdomain[in[i][j]] > 0 ) {

                            ibasket[1][ib[1]] = in[i][j];
                            ib[1] += 1;
                            icheck[in[i][j]][1] = 1;
                        }
                    }
                }
            }

            ib[0] = 0;

            for ( ix=0; ix<ib[1]; ix++ ) { // 2nd row flow boundary cells

                i = ibasket[1][ix];

                for ( j=0; j<9; j++ ) {

                    if ( icheck[in[i][j]][0] == 0 && cstopped[in[i][j]] == 0 && cdomain[in[i][j]] > 0 ) {

                        ibasket[0][ib[0]] = in[i][j];
                        ib[0] += 1;
                        icheck[in[i][j]][0] = 1;
                    }
                }
            }

            for ( ix=0; ix<ib[0]; ix++ ) { // updating domains and resetting checks

                i = ibasket[0][ix];

                if ( sico.MODEL <= 3 ) hflow = aw[i][0];
                else if ( sico.MODEL == 7 ) hflow = aw[i][0] + aw[i][3] + aw[i][6];

                if ( cdomain[i] != 0 ) { // for all non-edge cells

                    ctrlv = 0; ctrlr = 1; // resetting controls
                    for ( j=1; j<9; j++ ) { // loop over neighbourhood

                        if ( sico.MODEL <= 3 ) hflown = aw[in[i][j]][0];
                        else if ( sico.MODEL == 7 ) hflown = aw[in[i][j]][0] + aw[in[i][j]][3] + aw[in[i][j]][6];

                        if ( cflux[i] == 1 && ( cflux[in[i][j]] == 1 && ( hflown > sico.HFLOWMIN || hflow > sico.HFLOWMIN ))) ctrlv = 1;
                            // setting control for cell to positive if not at or surrounded by no-flux area

                        if ( cstopped[i] == 1 || cstopped[in[i][j]] == 1 ) ctrlr = 0;
                            // setting control for cell to negative if at or adjacent to stopped area
                    }
                    if ( ctrlv == 1 && ctrlr == 1 ) cdomain[i] = 3; // updating computational domain
                }

                if ( cdomain[i] == 3 ) { // if cell is part of computational domain

                    ctrlr = 0; // resetting control

                    for ( j=1; j<9; j++ ) { // loop over neighbourhood

                        if ( sico.MODEL <= 3 ) hflown = aw[in[i][j]][0];
                        else if ( sico.MODEL == 7 ) hflown = aw[in[i][j]][0] + aw[in[i][j]][3] + aw[in[i][j]][6];
                        if ( hflown > sico.HFLOWMIN ) ctrlr = 1; // if neighbour cell is flow cell, setting control to positive
                    }
                    if ( ctrlr == 1 ) cdomain[i] = 4; // flow depth > minimum around cell
                    if ( hflow > sico.HFLOWMIN ) cdomain[i] = 5; // flow depth > minimum
                }

                icheck[i][0] = 0;
                icheck[i][1] = 0;
            }


// -- STOP --- Updating computational domains -------------------------------------------------------------------


// *** Start of loop over various time step lengths until CFL criterion is met ----------------------------------


            ccfl = 0; // resetting control for fulfilment of CFL criterion
            while ( ccfl == 0 ) {


// -- START -- Applying initial block sliding -------------------------------------------------------------------


                // *** Modified mass point model, limited deformation, search radius and distance-dependent weights for mass point characteristics defined by option slidepar
                // *** Voellmy-type mixture model, one-parameter friction model (basal friction) with one-phase or multi-phase model, Pudasaini model (downslope acceleration only, deactivated)
                // *** Duration before slide changes to flow locally defined through parameter tslide (raster map) or, with tslide=-777, through presence of fluid

                if ( cslide == 1 ) {

                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];
                        if ( pxslide[i] == 1 ) { awt[i][0] = 0; awt[i][1] = 0; awt[i][2] = 0; } // initializing state variables
                    }

                    if ( sico.SLIDERAD <= 0 ) { // global centre coordinates and averages 

                        ansumh = 0; ansumslopex = 0; ansumslopey = 0; // initializing cumulative flow height and slopes
                   
                        for ( ix=0; ix<ib[0]; ix++ ) {

                            i = ibasket[0][ix];
                            if ( pxslide[i] == 1 && cdomain[i] == 5 ) {
                            
                                if ( sico.MODEL <= 3 ) hflow = aw[i][0];
                                else if ( sico.MODEL == 7 ) hflow = aw[i][0] + aw[i][3] + aw[i][6];

                                if ( sico.MODEL <= 3 ) {
                                                    
                                    anslopex = fbeta( pelev[in[i][5]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[i][5]][0] ), 
                                        pelev[in[i][2]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[i][2]][0] ), 2.0, sico ); // slopes
                                    anslopey = fbeta( pelev[in[i][4]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[i][4]][0] ), 
                                        pelev[in[i][1]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[i][1]][0] ), 2.0, sico );
                                                    
                                } else {

                                    anslopex = fbeta( pelev[in[i][5]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[i][5]][0] + aw[in[i][5]][3] + aw[in[i][5]][6] ), 
                                        pelev[in[i][2]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[i][2]][0] + aw[in[i][2]][3] + aw[in[i][2]][6] ), 2.0, sico );
                                    anslopey = fbeta( pelev[in[i][4]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[i][4]][0] + aw[in[i][4]][3] + aw[in[i][4]][6] ), 
                                        pelev[in[i][1]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[i][1]][0] + aw[in[i][1]][3] + aw[in[i][1]][6] ), 2.0, sico );
                                }

                                ansumh += hflow; // cumulative flow height
                                ansumslopex += hflow * tan(anslopex); // weighted cumulative x and y slopes
                                ansumslopey += hflow * tan(anslopey);
                            }
                        }
                    }

                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];
                        if ( pxslide[i] == 1 && cdomain[i] == 5 ) {

                            if ( sico.MODEL <= 3 ) { hflow = aw[i][0]; /*momx = aw[i][1]; momy = aw[i][2];*/ // flow heights and momenta
                            } else if ( sico.MODEL == 7 ) { 
                            
                                hflow = aw[i][0] + aw[i][3] + aw[i][6];
                                //momx = aw[i][1] + aw[i][4] + aw[i][7]; momy = aw[i][2] + aw[i][5] + aw[i][8];
                            }

                            if ( sico.SLIDERAD > 0 ) { // local centre coordinates and averages 

                                if ( sico.MODEL <= 3 || sico.PHASES[1] != 1 ) {

                                    andhdx = fbeta( aw[in[i][2]][0], aw[in[i][5]][0], 2.0, sico ); // gradients of flow height in x and y direction
                                    andhdy = fbeta( aw[in[i][1]][0], aw[in[i][4]][0], 2.0, sico );
                                
                                } else if ( sico.MODEL == 7 ) {
                                
                                    andhdx = fbeta( aw[in[i][2]][0] + aw[in[i][2]][3] + aw[in[i][2]][6], aw[in[i][5]][0] + aw[in[i][5]][3] + aw[in[i][5]][6], 2.0, sico );
                                    andhdy = fbeta( aw[in[i][1]][0] + aw[in[i][1]][3] + aw[in[i][1]][6], aw[in[i][4]][0] + aw[in[i][4]][3] + aw[in[i][4]][6], 2.0, sico );
                                }                         

                                ansumh = 0; ansumslopex = 0; ansumslopey = 0; // initializing cumulative flow height and slopes

                                iy = ffmax( 1, ( px[i] - anwin0 ) * sico.N ) + ffmax( 1, py[i] - anwin0 ); // loop over all cells within search window
                                for ( x=ffmax(1,px[i]-anwin0); x<ffmin((sico.M-1),px[i]+anwin0+1); x++ ) { 
                                    for ( y=ffmax(1,py[i]-anwin0); y<ffmin((sico.N-1),py[i]+anwin0+1); y++ ) {

                                        anrad = pow( pow( sico.CSZ * ( px[iy] - px[i] ), 2 ) + pow( sico.CSZ * ( py[iy] - py[i] ), 2 ), 0.5 ); // distance to target pixel
                                    
                                        if ( anrad < sico.SLIDERAD && pxslide[iy] == 1 && cdomain[iy] == 5 ) {
                            
                                            if ( sico.MODEL <= 3 ) hflown = aw[iy][0];
                                            else if ( sico.MODEL == 7 ) hflown = aw[iy][0] + aw[iy][3] + aw[iy][6];

                                            anwhtd = pow( 1 - anrad / sico.SLIDERAD, sico.SLIDEEXP );

                                            if ( sico.MODEL <= 3 ) {
                                                    
                                                anslopex = fbeta( pelev[in[iy][5]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[iy][5]][0] ), 
                                                    pelev[in[iy][2]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[iy][2]][0] ), 2.0, sico ); // slopes
                                                anslopey = fbeta( pelev[in[iy][4]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[iy][4]][0] ), 
                                                    pelev[in[iy][1]] +ffmax( 0, sico.SLIDEDEF ) * ( aw[in[iy][1]][0] ), 2.0, sico );                    
                                                    
                                            } else {

                                                anslopex = fbeta( pelev[in[iy][5]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[iy][5]][0] + aw[in[iy][5]][3] + aw[in[iy][5]][6] ), 
                                                    pelev[in[iy][2]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[iy][2]][0] + aw[in[iy][2]][3] + aw[in[iy][2]][6] ), 2.0, sico );
                                                anslopey = fbeta( pelev[in[iy][4]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[iy][4]][0] + aw[in[iy][4]][3] + aw[in[iy][4]][6] ), 
                                                    pelev[in[iy][1]] + ffmax( 0, sico.SLIDEDEF ) * ( aw[in[iy][1]][0] + aw[in[iy][1]][3] + aw[in[iy][1]][6] ), 2.0, sico );
                                            }

                                            ansumh += ( hflown * anwhtd ); // cumulative flow height
                                            ansumslopex += ( hflown * tan(anslopex) * anwhtd ); // weighted cumulative x and y slopes
                                            ansumslopey += ( hflown * tan(anslopey) * anwhtd );

                                        }
                                        iy += 1;
                                    }
                                    iy += ffmax( 2, sico.N - 2 * anwin0 - 1 );
                                }
                            }

                            anavgslopex = atan( ansumslopex / ansumh ) + ffmin( 0, sico.SLIDEDEF ) * andhdx; // average x and y slopes
                            anavgslopey = atan( ansumslopey / ansumh ) + ffmin( 0, sico.SLIDEDEF ) * andhdy;

                            anslope = fbetaxy( anavgslopex, anavgslopey ); // average slope, gravities, and aspect
                            angx = sico.GRAVITY * sin( anavgslopex );
                            angy = sico.GRAVITY * sin( anavgslopey );
                            // angxy = sico.GRAVITY * sin( anslope );
                            angz = sico.GRAVITY * cos( anslope );
                            // anaspect = falpha( anavgslopex, anavgslopey, sico );
                            
                            if ( sico.MODEL <= 3 ) { alpha1 = aw[i][0] / hflow; alpha2 = aw[i][3] / hflow; } else { alpha1 = 1; alpha2 = 0; }
                            // ansy0 = pow( pow( fdiv( momx, hflow, sico.HFLOWMIN ), 2 ) + pow( fdiv( momy, hflow, sico.HFLOWMIN ), 2 ), 0.5 );

                            if ( sico.DELTA == 1 ) { // mixture or PHASE 1 basal friction angle
                                if ( pdelta[i] != sico.UNDEF ) sflow.DELTA[0] = pdelta[i];
                                else sflow.DELTA[0] = sflow.DELTA0[0];
                            } else sflow.DELTA[0] = sflow.DELTA0[0];
                                
                            if ( sico.DELTA2 == 1 ) { // mixture or PHASE 2 basal friction angle
                                if ( pdelta2[i] != sico.UNDEF ) sflow.DELTA[1] = pdelta2[i];
                                else sflow.DELTA[1] = sflow.DELTA0[1];
                            } else sflow.DELTA[1] = sflow.DELTA0[1];

                            if ( sico.TUFRI == 1 ) { // turbulent friction coefficient
                                if ( ptufri[i] != sico.UNDEF ) sflow.TUFRI = ptufri[i];
                                else sflow.TUFRI = sflow.TUFRI0;
                            } else sflow.TUFRI = sflow.TUFRI0;
                               
                            if ( sico.PHASES[0] != 3 ) andelta = alpha1 * sflow.DELTA[0] + alpha2 * sflow.DELTA[1]; // basal friction angle of mixture
                            else andelta = 0;

                            //if ( sico.MODEL <= 3 || aw[i][0] + aw[i][3] <= sico.HFLOWMIN || aw[i][6] <= sico.HFLOWMIN ) { analphas = 1; angamma = 0; // solid and fluid fractions and density ratios
                            //} else if ( sico.MODEL == 7 ) { analphas = alpha1 + alpha2; angamma = ( alpha1 * sflow.RHO1 + alpha2 * sflow.RHO2 ) / ( analphas * sflow.RHO3 ); }

                            for ( l=0; l<1; l++ ) { // 9 for Pudasaini model

                                if ( sico.MODEL <= 3 || sico.PHASES[1] != 1 ) {

                                    wu[l] = fdiv( aw[in[i][l]][1], aw[in[i][l]][0], sico.HFLOWMIN ); // velocities of mixture in x and y direction
                                    wv[l] = fdiv( aw[in[i][l]][2], aw[in[i][l]][0], sico.HFLOWMIN );
                                    
                                } else if ( sico.MODEL == 7 ) {
                                    
                                    wu[l] = fdiv( aw[in[i][l]][1] + aw[in[i][l]][4] + aw[in[i][l]][7], aw[in[i][l]][0] + aw[in[i][l]][3] + aw[in[i][l]][6], sico.HFLOWMIN );
                                    wv[l] = fdiv( aw[in[i][l]][2] + aw[in[i][l]][5] + aw[in[i][l]][8], aw[in[i][l]][0] + aw[in[i][l]][3] + aw[in[i][l]][6], sico.HFLOWMIN );
                                }
                            }

                            /*wdu[0] = ( wu[2] - wu[5] ) / ( 2 * dx[i] ); // velocity gradients in x and y direction
                            wdv[0] = ( wv[1] - wv[4] ) / ( 2 * dy[i] );
                                
                            if ( sico.MODEL <= 3 || sico.PHASES[1] != 1 ) {
                                
                                hekin = 0.5 * hflow * sflow.RHO1 * ( pow( wu[0], 2 ) + pow( wv[0], 2 )); // kinetic energy, solid flow height, and earth pressure coefficients
                                hflows = aw[i][0];
                                ank = fk( hflow, hflows, wu[0], wv[0], wdu[0]+wdv[0], hekin, 0, sico, sflow );
                                    
                            } else if ( sico.PHASES[0] != 3 ) {
                                
                                hekin = 0.5 * hflow * ( alpha1 * sflow.RHO1 + alpha2 * sflow.RHO2 + ( 1 - alpha1 - alpha2 ) * sflow.RHO3 ) * ( pow( wu[0], 2 ) + pow( wv[0], 2 ));
                                hflows = aw[i][0] + aw[i][3];
                                ank = fk( hflow, hflows, alpha1/analphas, alpha2/analphas, wdv[0]+wdv[0], hekin, 9, sico, sflow );
                                    
                            } else { hflows = 0; ank = 0; }

                            andhdxy = fbetaxy( andhdx, andhdy );
                            anaspectdh = falpha( andhdx, andhdy, sico );
                            andh = tan( andhdxy * cos( anaspectdh - anaspect )); // flow height gradient in downslope direction
                            analpha = angxy - ( 1 - angamma ) * analphas * angz * tan( andelta ) - angz * ((( 1 - angamma ) * ank + angamma ) * analphas + (1 - analphas )) * andh;
                            anbeta = 0.001; // alpha and beta parameters

                            if ( analpha > 0 && pow( anbeta / analpha, 0.5 ) * ansy0 <= 1 ) {

                                anu[i] = pow( analpha / anbeta, 0.5 ) * tanh( pow( analpha * anbeta, 0.5 ) * tlength + atanh( pow( anbeta / analpha, 0.5 ) * ansy0 )); // sliding velocity
                                anx[i] = 1 / anbeta * log( cosh( pow( analpha * anbeta, 0.5 ) * tlength + atanh( pow( anbeta / analpha, 0.5 ) * ansy0 ))) 
                                    - 1 / anbeta * log( cosh( atanh( pow( anbeta / analpha, 0.5 ) * ansy0 ))); // sliding distance in downslope direction
                            }
                                
                            anpx = anx[i] * cos( anaspect ); 
                            anpy = anx[i] * sin( anaspect );
                            anupx = anu[i] * cos( anaspect ); 
                            anupy = anu[i] * sin( anaspect );*/

                            vflow = pow( pow( wu[0], 2 ) + pow( wv[0], 2 ), 0.5 ); // flow velocity and velocity ratios

                            if ( vflow != 0 ) {
                                vflowxratio = wu[0] / vflow;
                                vflowyratio = wv[0] / vflow;
                                    
                            } else { vflowxratio = 0; vflowyratio = 0; }                               

                            if ( sico.MODEL == 0 ) { // momentum production with Voellmy-type model
               
                                anupx = wu[0] + ( angx * hflow - vflowxratio * ( tan( andelta ) * angz * hflow
                                    + sico.GRAVITY * pow( vflow, 2 ) / sflow.TUFRI )) / hflow * tlength;
                                anupy = wv[0] + ( angy * hflow - vflowyratio * ( tan( andelta ) * angz * cos( anslope ) * hflow
                                    + sico.GRAVITY * pow( vflow, 2 ) / sflow.TUFRI )) / hflow * tlength;
                                
                            } else { // momentum production with one-parameter friction model
                                
                                anupx = wu[0] + ( angx - vflowxratio * tan( andelta ) * angz ) * tlength;
                                anupy = wv[0] + ( angy - vflowyratio * tan( andelta ) * angz ) * tlength;
                            }

                            anpx = px[i] + anupx * tlength / sico.CSZ; // sliding distances in raster cells in x and y direction
                            anpy = py[i] + anupy * tlength / sico.CSZ;

                            anpx1 = (int)( anpx ); anpx2 = (int)( anpx ) + 1; anpy1 = (int)( anpy ); anpy2 = (int)( anpy ) + 1; // nearest neighbour cells
                            anwht[0] = ( anpx2 - anpx ) * ( anpy2 - anpy ); anwht[1] = ( anpx2 - anpx ) * ( anpy - anpy1 );
                            anwht[2] = ( anpx - anpx1 ) * ( anpy2 - anpy ); anwht[3] = ( anpx - anpx1 ) * ( anpy - anpy1 );

                            andist = (int)( tlength * ffmax(anupx, anupy) / sico.CSZ + 2 ); anctr = 0;

                            iz = ffmax( 1, ( px[i] - andist ) * sico.N ) + ffmax( 1, py[i] - andist ); // loop over all relevant cells
                            for ( x=ffmax(1,px[i]-andist); x<ffmin((sico.M-1),px[i]+andist+1); x++ ) { 
                                for ( y=ffmax(1,py[i]-andist); y<ffmin((sico.N-1),py[i]+andist+1); y++ ) {

                                    anid = -1;

                                    if ( px[iz] == anpx1 && py[iz] == anpy1 ) { anid = 0; anctr += 1; } // searching for target cells
                                    if ( px[iz] == anpx1 && py[iz] == anpy2 ) { anid = 1; anctr += 1; }
                                    if ( px[iz] == anpx2 && py[iz] == anpy1 ) { anid = 2; anctr += 1; }
                                    if ( px[iz] == anpx2 && py[iz] == anpy2 ) { anid = 3; anctr += 1; }

                                    if ( anid >= 0 ) {

                                        awt[iz][0] += aw[i][0] * anwht[anid]; // assigning flow heights and momenta to target cells
                                        awt[iz][1] += aw[i][0] * anwht[anid] * anupx; 
                                        awt[iz][2] += aw[i][0] * anwht[anid] * anupy;
                                        
                                        if ( sico.MODEL == 7 ) {
                                        
                                            awt[iz][3] += aw[i][3] * anwht[anid]; 
                                            awt[iz][4] += aw[i][3] * anwht[anid] * anupx; 
                                            awt[iz][5] += aw[i][3] * anwht[anid] * anupy;
                                            
                                            awt[iz][6] += aw[i][6] * anwht[anid]; 
                                            awt[iz][7] += aw[i][6] * anwht[anid] * anupx; 
                                            awt[iz][8] += aw[i][6] * anwht[anid] * anupy;
                                        }
                                    }
                                    iz += 1;
                                }
                                iz += ffmax( 2, sico.N - 2 * andist - 1 );
                            }
                        }
                    }
                }

                
// -- STOP --- Applying initial block sliding -------------------------------------------------------------------


// -- START -- Flow propagation with NOC scheme -----------------------------------------------------------------


                // *** High-resolution scheme, system is moved half a cell size each time step
                // *** Technically, the system is shifted one cell each second time step
                // *** At the end, the new state variables are written to a temporary vector (awt[i][k])
                // *** The values from the temporary vector are later written to the permanent vector (aw[i][k])
                // *** as soon as it has been verified that the time step length is sufficiently short to fulfil the CFL criterion
                // *** Otherwise, the procedure is repeated with shorter time step length


                // Fluxes, source terms, and gradients (original coordinate system)

                if ( cflow == 1 ) {

                  for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] <= 3 ) {

                        for ( k=0; k<sico.NVECTMIN; k++ ) { af[i][k]=0; ag[i][k]=0; as[i][k]=0; asigma_x[i][k]=0; asigma_y[i][k]=0; }

                    } else {

                        if ( sico.PHI == 1 ) {
                            if ( pphi[i] != sico.UNDEF ) sflow.PHI[0] = pphi[i]; // internal friction angle of mixture or PHASE 1
                            else sflow.PHI[0] = sflow.PHI0[0];
                        } else sflow.PHI[0] = sflow.PHI0[0];

                        if ( sico.PHI2 == 1 ) {
                            if ( pphi2[i] != sico.UNDEF ) sflow.PHI[1] = pphi2[i]; // internal friction angle of PHASE 2
                            else sflow.PHI[1] = sflow.PHI0[1];
                        } else sflow.PHI[1] = sflow.PHI0[1];

                        if ( sico.PHI3 == 1 ) {
                            if ( pphi3[i] != sico.UNDEF ) sflow.PHI[2] = pphi3[i]; // internal friction angle of PHASE 3
                            else sflow.PHI[2] = sflow.PHI0[2];
                        } else sflow.PHI[2] = sflow.PHI0[2];

                        if ( sico.DELTA == 1 ) {
                            if ( pdelta[i] != sico.UNDEF ) sflow.DELTA[0] = pdelta[i]; // basal friction angle of mixture or PHASE 1
                            else sflow.DELTA[0] = sflow.DELTA0[0];
                        } else sflow.DELTA[0] = sflow.DELTA0[0];

                        if ( sico.DELTA2 == 1 ) {
                            if ( pdelta2[i] != sico.UNDEF ) sflow.DELTA[1] = pdelta2[i]; // basal friction angle of PHASE 2
                            else sflow.DELTA[1] = sflow.DELTA0[1];
                        } else sflow.DELTA[1] = sflow.DELTA0[1];

                        if ( sico.DELTA3 == 1 ) {
                            if ( pdelta3[i] != sico.UNDEF ) sflow.DELTA[2] = pdelta3[i]; // basal friction angle of PHASE 3
                            else sflow.DELTA[2] = sflow.DELTA0[2];
                        } else sflow.DELTA[2] = sflow.DELTA0[2];
                        
                        if ( sico.TUFRI == 1 ) {
                            if ( ptufri[i] != sico.UNDEF ) sflow.TUFRI = ptufri[i]; // turbulent friction coefficient
                            else sflow.TUFRI = sflow.TUFRI0;
                        } else sflow.TUFRI = sflow.TUFRI0;

                        if ( frictiograph == 1 ) {

                            for ( frik=0; frik<=fritmax; frik++ ) { // identifying relevant line of input frictiograph
                                if ( frifri[frik][0] <= tsum ) frit=frik;
                                else break;
                            }

                            sflow.PHI[0] = frifri[frit][1] * sico.PI / 180; sflow.DELTA[0] = frifri[frit][2] * sico.PI / 180;
                            if ( sico.MODEL == 7 ) { sflow.PHI[1] = frifri[frit][3] * sico.PI / 180; sflow.DELTA[1] = frifri[frit][4] * sico.PI / 180; }
                            if ( sico.MODEL == 7 ) { sflow.PHI[2] = frifri[frit][5] * sico.PI / 180; sflow.DELTA[2] = frifri[frit][6] * sico.PI / 180; }
                        }

                        if ( sico.NYSS == 1 ) {
                            if ( pnyss[i] != sico.UNDEF ) sflow.NY[0] = pnyss[i]; // viscosity of mixture or PHASE 1
                            else sflow.NY[0] = sflow.NY0[0];
                        } else sflow.NY[0] = sflow.NY0[0];

                        if ( sico.NYFS == 1 ) {
                            if ( pnyfs[i] != sico.UNDEF ) sflow.NY[1] = pnyfs[i]; // viscosity of PHASE 2
                            else sflow.NY[1] = sflow.NY0[1];
                        } else sflow.NY[1] = sflow.NY0[1];

                        if ( sico.NYFF == 1 ) {
                            if ( pnyff[i] != sico.UNDEF ) sflow.NY[2] = pnyff[i]; // viscosity of PHASE 3
                            else sflow.NY[2] = sflow.NY0[2];
                        } else sflow.NY[2] = sflow.NY0[2];

                        if ( sico.AMBDRAG == 1 ) {
                            if ( pambdrag[i] != sico.UNDEF ) sflow.AMBDRAG = pambdrag[i]; // ambient drag coefficient
                            else sflow.AMBDRAG = sflow.AMBDRAG0;
                        } else sflow.AMBDRAG = sflow.AMBDRAG0;

                        if ( sico.FLUFRI == 1 ) {
                            if ( pflufri[i] != sico.UNDEF ) sflow.FLUFRI = pflufri[i]; // fluid friction coefficient
                            else sflow.FLUFRI = sflow.FLUFRI0;
                        } else sflow.FLUFRI = sflow.FLUFRI0;


                        // Components of gravity

                        grav[0] = sico.GRAVITY * cos( betax[i] );
                        grav[1] = sico.GRAVITY * cos( betay[i] );
                        grav[2] = sico.GRAVITY * sin( betax[i] );
                        grav[3] = sico.GRAVITY * sin( betay[i] );
                        grav[4] = sico.GRAVITY * cos( betaxy[i] );


                        // Flow velocities in x and y directions

                        vflowx = fdiv( aw[i][1], aw[i][0], sico.HFLOWMIN );
                        vflowy = fdiv( aw[i][2], aw[i][0], sico.HFLOWMIN );

                        if ( sico.MODEL == 7 ) vflowx2 = fdiv( aw[i][4], aw[i][3], sico.HFLOWMIN ); else vflowx2 = 0;
                        if ( sico.MODEL == 7 ) vflowy2 = fdiv( aw[i][5], aw[i][3], sico.HFLOWMIN ); else vflowy2 = 0;

                        if ( sico.MODEL == 7 ) vflowx3 = fdiv( aw[i][7], aw[i][6], sico.HFLOWMIN ); else vflowx3 = 0;
                        if ( sico.MODEL == 7 ) vflowy3 = fdiv( aw[i][8], aw[i][6], sico.HFLOWMIN ); else vflowy3 = 0;
                        
                        
                        // Total flow depths, velocities, fractions, and slopes of neighbouring cells

                        welev[0] = pelev[i];

                        for ( l=0; l<9; l++ ) {

                            welev[l] = pelev[in[i][l]];

                            nbetax[l] = betax[in[i][l]];
                            nbetay[l] = betay[in[i][l]];

                            if ( sico.MODEL <= 3 ) wh[l] = aw[in[i][l]][0];
                            else if ( sico.MODEL == 7 ) wh[l] = aw[in[i][l]][0] + aw[in[i][l]][3] + aw[in[i][l]][6];

                            wu[l] = fdiv( aw[in[i][l]][1], aw[in[i][l]][0], sico.HFLOWMIN );
                            wv[l] = fdiv( aw[in[i][l]][2], aw[in[i][l]][0], sico.HFLOWMIN );

                            if ( sico.MODEL == 7 ) wu2[l] = fdiv( aw[in[i][l]][4], aw[in[i][l]][3], sico.HFLOWMIN ); else wu2[l] = 0;
                            if ( sico.MODEL == 7 ) wv2[l] = fdiv( aw[in[i][l]][5], aw[in[i][l]][3], sico.HFLOWMIN ); else wv2[l] = 0;

                            if ( sico.MODEL == 7 ) wu3[l] = fdiv( aw[in[i][l]][7], aw[in[i][l]][6], sico.HFLOWMIN ); else wu3[l] = 0;
                            if ( sico.MODEL == 7 ) wv3[l] = fdiv( aw[in[i][l]][8], aw[in[i][l]][6], sico.HFLOWMIN ); else wv3[l] = 0;

                            walpha[l] = fdiv( aw[in[i][l]][0], wh[l], sico.HFLOWMIN );
                            if ( sico.MODEL == 7 ) walpha2[l] = fdiv( aw[in[i][l]][3], wh[l], sico.HFLOWMIN ); else walpha2[l] = 0;
                            if ( sico.MODEL == 7 ) walpha3[l] = fdiv( aw[in[i][l]][6], wh[l], sico.HFLOWMIN ); else walpha3[l] = 0;
                        }


                        // Gradients of total flow depth, flow velocities, and fractions

                        whx = ( wh[2] - wh[5] ) / ( 2 * dx[i] );
                        if ( sico.NOOSC > 0 && sico.NOOSC < 5 && ( wh[2] == 0 || wh[5] == 0 )) whx = tan( betax[i] ); // balance of forces at reservoir edge

                        why = ( wh[1] - wh[4] ) / ( 2 * dy[i] );
                        if ( sico.NOOSC > 0 && sico.NOOSC < 5 && ( wh[1] == 0 || wh[4] == 0 )) why = tan( betay[i] ); 

                        wdu[0] = ( wu[2] - wu[5] ) / ( 2 * dx[i] );
                        wdv[0] = ( wv[1] - wv[4] ) / ( 2 * dy[i] );

                        if ( sico.MODEL == 7 ) wdu[1] = ( wu2[2] - wu2[5] ) / ( 2 * dx[i] ); else wdu[1] = 0;
                        if ( sico.MODEL == 7 ) wdv[1] = ( wv2[1] - wv2[4] ) / ( 2 * dy[i] ); else wdv[1] = 0;

                        if ( sico.MODEL == 7 ) wdu[2] = ( wu3[2] - wu3[5] ) / ( 2 * dx[i] ); else wdu[2] = 0;
                        if ( sico.MODEL == 7 ) wdv[2] = ( wv3[1] - wv3[4] ) / ( 2 * dy[i] ); else wdv[2] = 0;

                        if ( walpha[2] ) xwdv[0] = ( wv[2] - wv[5] ) / ( 2 * dx[i] );
                        if ( walpha[1] ) xwdu[0] = ( wu[1] - wu[4] ) / ( 2 * dy[i] );

                        if ( sico.MODEL == 7 ) xwdv[1] = ( wv2[2] - wv2[5] ) / ( 2 * dx[i] ); else xwdv[1] = 0;
                        if ( sico.MODEL == 7 ) xwdu[1] = ( wu2[1] - wu2[4] ) / ( 2 * dy[i] ); else xwdu[1] = 0;

                        if ( sico.MODEL == 7 ) xwdv[2] = ( wv3[2] - wv3[5] ) / ( 2 * dx[i] ); else xwdv[2] = 0;
                        if ( sico.MODEL == 7 ) xwdu[2] = ( wu3[1] - wu3[4] ) / ( 2 * dy[i] ); else xwdu[2] = 0;

                        walphax = ( walpha[2] - walpha[5] ) / ( 2 * dx[i] );
                        walphay = ( walpha[1] - walpha[4] ) / ( 2 * dy[i] );

                        if ( sico.MODEL == 7 ) walphax2 = ( walpha2[2] - walpha2[5] ) / ( 2 * dx[i] ); else walphax2 = 0;
                        if ( sico.MODEL == 7 ) walphay2 = ( walpha2[1] - walpha2[4] ) / ( 2 * dy[i] ); else walphay2 = 0;

                        if ( sico.MODEL == 7 ) walphax3 = ( walpha3[2] - walpha3[5] ) / ( 2 * dx[i] ); else walphax3 = 0;
                        if ( sico.MODEL == 7 ) walphay3 = ( walpha3[1] - walpha3[4] ) / ( 2 * dy[i] ); else walphay3 = 0;


                        // Earth pressure coefficients

                        if ( sico.MODEL <= 3 ) { whflow = aw[i][0], whflow2 = 0, whflow3 = 0; }
                        else if ( sico.MODEL == 7 ) { whflow = aw[i][0], whflow2 = aw[i][3], whflow3 = aw[i][6]; }

                        hekin = 0.5 * whflow * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 ));
                        if ( sico.MODEL == 7 ) hekin += 0.5 * whflow2 * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                        if ( sico.MODEL == 7 ) hekin += 0.5 * whflow3 * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));

                        for ( l=0; l<sico.PMAX; l++ ) {

                            if ( sico.PHASES[l] < 2 ) {

                                gkx[l] = fk( whflow+whflow2+whflow3, aw[i][3*l], vflowx, vflowy, wdu[l], hekin, l, sico, sflow );
                                gky[l] = fk( whflow+whflow2+whflow3, aw[i][3*l], vflowx, vflowy, wdv[l], hekin, l, sico, sflow );
                            }
                        }

                        for ( k=0; k<18; k++ ) wwd[k] = ad[i][k+sico.NVECTMIN];


                        // Curvature, flux, source, and deceleration terms

                        kappau = fcurv( vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, welev, grav, sico );
                        vm = fvm( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, sico, sflow );
                        cdrag = fdrag( whflow, whflow2, whflow3, sico, sflow );
                        gze = fgze( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wh, whx, why,
                            wdu, wdv, xwdu, xwdv, wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, dx[i], dy[i],
                            grav, betax[i], betay[i], cdrag, tsum, kappau, sico, sflow );
                        disp = fdisp( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wdu, wdv, xwdu, xwdv, 
                            betax[i], betay[i], gkx, gky, vm, cdrag, tsum, sico, sflow );

                        gf = ff( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                            gkx, gze, kappau, vm, whx, walphax, cdrag, grav, disp, betax[i], sflow, sico );
                        gg = fg( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                            gky, gze, kappau, vm, why, walphay, cdrag, grav, disp, betay[i], sflow, sico );
                        gs = fs( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                            whx, why, walphax, walphay, walphax2, walphay2, walphax3, walphay3, grav, gze, betax[i], betay[i], kappau, cdrag, sico, sflow );
                        gdecel = fd( wh, whx, why, wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, nbetax, nbetay, gze, dx[i], dy[i],
                            wwd, kappau, hekin, grav, sico, sflow );

                        for ( k=0; k<sico.NVECTMIN; k++ ) { 
                            af[i][k] = gf[k]; ag[i][k] = gg[k]; as[i][k] = gs[k];
                        }
                        for ( k=0; k<sico.NVECTMIN+18; k++ ) ad[i][k] = gdecel[k];

                        free( kappau ); free( vm ); free( cdrag ); free( gze ); free( disp ); free( gf ); free( gg ); free( gs ); free( gdecel );


                        // Slopes of the vector components                        

                        asigma_xelev[i] = fsigma( pelev[i], pelev[in[i][5]], pelev[in[i][2]], 1, dx, dy, i, sico );
                        asigma_yelev[i] = fsigma( pelev[i], pelev[in[i][4]], pelev[in[i][1]], 2, dx, dy, i, sico );

                        for (k=0; k<sico.NVECTMIN; k++) {

                            asigma_x[i][k] = fsigma( aw[i][k], aw[in[i][5]][k], aw[in[i][2]][k], 1, dx, dy, i, sico );
                            asigma_y[i][k] = fsigma( aw[i][k], aw[in[i][4]][k], aw[in[i][1]][k], 2, dx, dy, i, sico );
                        }
                    }
                }


                // Slopes of the fluxes

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] <= 3 ) {

                        for ( k=0; k<sico.NVECTMIN; k++) { asigma_f[i][k]=0; asigma_g[i][k]=0; }
		
                    } else {

                        for ( k=0; k<sico.NVECTMIN; k++ ) {

                            asigma_f[i][k] = fsigma ( af[i][k], af[in[i][5]][k], af[in[i][2]][k], 1, dx, dy, i, sico );
                            asigma_g[i][k] = fsigma ( ag[i][k], ag[in[i][4]][k], ag[in[i][1]][k], 2, dx, dy, i, sico );
                        }
                    }
                }


                // Values of vector at quarter of cell after half time step

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] > 3 ) {

                        fj[0][0]=1; fj[0][1]=1; fj[1][0]=1; fj[1][1]=-1; fj[2][0]=-1; fj[2][1]=1; fj[3][0]=-1; fj[3][1]=-1; // factors for positive or negative gradients
                        for ( j=0; j<4; j++ ) wintelev[i][j]=pelev[in[i][j]]+fj[j][0]*0.25*dx[i]*asigma_xelev[in[i][j]]+fj[j][1]*0.25*dy[i]*asigma_yelev[in[i][j]];
                    }
                }

                ctrl_noosc = 0;
                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] > 3 ) {

                        fj[0][0]=1; fj[0][1]=1; fj[1][0]=1; fj[1][1]=-1; fj[2][0]=-1; fj[2][1]=1; fj[3][0]=-1; fj[3][1]=-1;
                            // factors for positive or negative gradients

                        for ( j=0; j<6; j++ ) { // loop over adjacent cells

                            for ( k=0; k<sico.NVECTMIN; k++ ) { // loop over all components of the vector

                                dw_dt[k]=-asigma_f[in[i][j]][k]-asigma_g[in[i][j]][k]+as[in[i][j]][k]; // gradient of vector
                                dw_dttest[k]=-asigma_f[in[i][j]][k]-asigma_g[in[i][j]][k]+(as[in[i][j]][k]-ad[in[i][j]][k]); // gradient of vector
                            }

                            for ( k=0; k<sico.NVECTMIN; k++ ) { // loop over all components of the vector

                                if ( j < 4 ) winta[i][j][k]=aw[in[i][j]][k]+fj[j][0]*0.25*dx[i]*asigma_x[in[i][j]][k]+fj[j][1]
                                    *0.25*dy[i]*asigma_y[in[i][j]][k]; // value at quarter of cell

                                wintb[i][j][k]=aw[in[i][j]][k]+0.5*tlength*dw_dt[k];
                                wintbtest=aw[in[i][j]][k]+0.5*tlength*dw_dttest[k];

                                if (( fsign(wintb[i][j][k]) == fsign(wintbtest) && fabs(wintb[i][j][k]) > fabs(wintbtest)) || sico.ORIGINAL == 1 ) wintb[i][j][k] = wintbtest;
                                else if ( fsign(wintb[i][j][k]) == fsign(wintbtest)) wintb[i][j][k] = wintb[i][j][k];
                                else if ( k==0 || k==3 || k==6 ) wintb[i][j][k] = aw[in[i][j]][k];
                                else wintb[i][j][k] = 0;

                                if ( j < 4 ) {

                                    wintc[i][j][k]=aw[in[i][j]][k]+0.5*tlength*dw_dt[k]+fj[j][0]*0.25*dx[i]*asigma_x[in[i][j]][k]
                                    +fj[j][1]*0.25*dy[i]*asigma_y[in[i][j]][k];

                                    wintctest=aw[in[i][j]][k]+0.5*tlength*dw_dttest[k]+fj[j][0]*0.25*dx[i]*asigma_x[in[i][j]][k]
                                    +fj[j][1]*0.25*dy[i]*asigma_y[in[i][j]][k];

                                    if (( fsign(wintc[i][j][k]) == fsign(wintctest) && fabs(wintc[i][j][k]) > fabs(wintctest)) || sico.ORIGINAL==1 ) wintc[i][j][k] = wintctest;
                                    else if ( fsign(wintc[i][j][k]) == fsign(wintctest)) wintc[i][j][k] = wintc[i][j][k];
                                    else if ( k==0 || k==3 || k==6 ) wintc[i][j][k] = aw[in[i][j]][k];
                                    else wintc[i][j][k] = 0;

                                }
                            }
                        }

                        if ( sico.NOOSC > 0 && sico.NOOSC < 5 ) { // balance of forces at reservoir edge

                            ctrl_lake = 0;
                            minlevel = 99999;

                            for ( j=0; j<4; j++ ) {

                                if ( sico.MODEL <= 3 ) {

                                    wh[j] = aw[in[i][j]][0];
                                    wmx = aw[in[i][j]][1];
                                    wmy = aw[in[i][j]][2];

                                } else if ( sico.MODEL == 7 ) {

                                    wh[j] = aw[in[i][j]][0] + aw[in[i][j]][3] + aw[in[i][j]][6];
                                    wmx = aw[in[i][j]][1] + aw[in[i][j]][4] + aw[in[i][j]][7];
                                    wmy = aw[in[i][j]][2] + aw[in[i][j]][5] + aw[in[i][j]][8];
                                }

                                if ( wh[j] > 0 && wmx == 0 && wmy == 0 ) {

                                    ctrl_lake = 1;
                                    if ( wh[j] + pelev[in[i][j]] < minlevel ) minlevel = wh[j] + pelev[in[i][j]];
                                }
                            }

                            if ( ctrl_lake == 1 ) {

                                ctrl_noosc = 1;
                                for ( j=0; j<4; j++ ) {

                                    if ( sico.MODEL <= 3 ) {

                                        wh[j] = aw[i][0] - ( wintelev[i][j] - pelev[i] );
                                        walpha[j] = 1.0;

                                    } else if ( sico.MODEL == 7 ) {

                                        wh[j] = aw[i][0] + aw[i][3] + aw[i][6] - ( wintelev[i][j] - pelev[i] );
                                        walpha[j] = fdiv( aw[i][0], aw[i][0] + aw[i][3] + aw[i][6], 0 );
                                        walpha2[j] = fdiv( aw[i][3], aw[i][0] + aw[i][3] + aw[i][6], 0 );
                                        walpha3[j] = fdiv( aw[i][6], aw[i][0] + aw[i][3] + aw[i][6], 0 );
                                    }

                                    wintc[i][j][0] = wh[j] * walpha[j];
                                    if ( wintelev[i][j] + wh[j] > minlevel ) wintc[i][j][0] = ( minlevel - wintelev[i][j] ) * walpha[j];
                                    if ( wintc[i][j][0] < 0 ) wintc[i][j][0] = 0;

                                    if ( sico.MODEL == 7 ) {

                                        wintc[i][j][3] = wh[j] * walpha2[j];
                                        if ( wintelev[i][j] + wh[j] > minlevel ) wintc[i][j][3] = ( minlevel - wintelev[i][j] ) * walpha2[j];
                                        if ( wintc[i][j][3] < 0 ) wintc[i][j][3] = 0;

                                        wintc[i][j][6] = wh[j] * walpha3[j];
                                        if ( wintelev[i][j] + wh[j] > minlevel ) wintc[i][j][6] = ( minlevel - wintelev[i][j] ) * walpha3[j];
                                        if ( wintc[i][j][6] < 0 ) wintc[i][j][6] = 0;
                                    }
                                }
                            }
                        }

                    } else {

                        for ( j=0; j<4; j++ ) {
                            for ( k=0; k<sico.NVECTMIN; k++ ) {

                                if ( j < 4 ) winta[i][j][k]=0; //!!!CHECK j
                                wintb[i][j][k]=0; 
                                if ( j < 4 ) wintc[i][j][k]=0;
                            }
                        }
                    }
                }


                // Fluxes and source terms (shifted coordinate system)

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] > 3 ) {

                        if ( sico.PHI == 1 ) {
                            if ( pphi[i] != sico.UNDEF ) sflow.PHI[0] = pphi[i]; // internal friction angle of mixture or PHASE 1
                            else sflow.PHI[0] = sflow.PHI0[0];
                        } else sflow.PHI[0] = sflow.PHI0[0];

                        if ( sico.PHI2 == 1 ) {
                            if ( pphi2[i] != sico.UNDEF ) sflow.PHI[1] = pphi2[i]; // internal friction angle of PHASE 2
                            else sflow.PHI[1] = sflow.PHI0[1];
                        } else sflow.PHI[1] = sflow.PHI0[1];

                        if ( sico.PHI3 == 1 ) {
                            if ( pphi3[i] != sico.UNDEF ) sflow.PHI[2] = pphi3[i]; // internal friction angle of PHASE 3
                            else sflow.PHI[2] = sflow.PHI0[2];
                        } else sflow.PHI[2] = sflow.PHI0[2];

                        if ( sico.DELTA == 1 ) {
                            if ( pdelta[i] != sico.UNDEF ) sflow.DELTA[0] = pdelta[i]; // basal friction angle of mixture or PHASE 1
                            else sflow.DELTA[0] = sflow.DELTA0[0];
                        } else sflow.DELTA[0] = sflow.DELTA0[0];

                        if ( sico.DELTA2 == 1 ) {
                            if ( pdelta2[i] != sico.UNDEF ) sflow.DELTA[1] = pdelta2[i]; // basal friction angle of PHASE 2
                            else sflow.DELTA[1] = sflow.DELTA0[1];
                        } else sflow.DELTA[1] = sflow.DELTA0[1];

                        if ( sico.DELTA3 == 1 ) {
                            if ( pdelta3[i] != sico.UNDEF ) sflow.DELTA[2] = pdelta3[i]; // basal friction angle of PHASE 3
                            else sflow.DELTA[2] = sflow.DELTA0[2];
                        } else sflow.DELTA[2] = sflow.DELTA0[2];

                        if ( sico.TUFRI == 1 ) {
                            if ( ptufri[i] != sico.UNDEF ) sflow.TUFRI = ptufri[i]; // turbulent friction coefficient
                            else sflow.TUFRI = sflow.TUFRI0;
                        } else sflow.TUFRI = sflow.TUFRI0;

                        if ( frictiograph == 1 ) {

                            for ( frik=0; frik<=fritmax; frik++ ) { // identifying relevant line of input frictiograph
                                if ( frifri[frik][0] <= tsum ) frit=frik;
                                else break;
                            }

                            sflow.PHI[0] = frifri[frit][1] * sico.PI / 180; sflow.DELTA[0] = frifri[frit][2] * sico.PI / 180;
                            if ( sico.MODEL == 7 ) { sflow.PHI[1] = frifri[frit][3] * sico.PI / 180; sflow.DELTA[1] = frifri[frit][4] * sico.PI / 180; }
                            if ( sico.MODEL == 7 ) { sflow.PHI[2] = frifri[frit][5] * sico.PI / 180; sflow.DELTA[2] = frifri[frit][6] * sico.PI / 180; }
                        }

                        if ( sico.NYSS == 1 ) {
                            if ( pnyss[i] != sico.UNDEF ) sflow.NY[0] = pnyss[i]; // viscosity of PHASE 1
                            else sflow.NY[0] = sflow.NY0[0];
                        } else sflow.NY[0] = sflow.NY0[0];

                        if ( sico.NYFS == 1 ) {
                            if ( pnyfs[i] != sico.UNDEF ) sflow.NY[1] = pnyfs[i]; // viscosity of PHASE 2
                            else sflow.NY[1] = sflow.NY0[1];
                        } else sflow.NY[1] = sflow.NY0[1];

                        if ( sico.NYFF == 1 ) {
                            if ( pnyff[i] != sico.UNDEF ) sflow.NY[2] = pnyff[i]; // viscosity of PHASE 3
                            else sflow.NY[2] = sflow.NY0[2];
                        } else sflow.NY[2] = sflow.NY0[2];

                        if ( sico.AMBDRAG == 1 ) {
                            if ( pambdrag[i] != sico.UNDEF ) sflow.AMBDRAG = pambdrag[i]; // ambient drag coefficient
                            else sflow.AMBDRAG = sflow.AMBDRAG0;
                        } else sflow.AMBDRAG = sflow.AMBDRAG0;

                        if ( sico.FLUFRI == 1 ) {
                            if ( pflufri[i] != sico.UNDEF ) sflow.FLUFRI = pflufri[i]; // fluid friction coefficient
                            else sflow.FLUFRI = sflow.FLUFRI0;
                        } else sflow.FLUFRI = sflow.FLUFRI0;

                        for ( j=0; j<4; j++ ) {


                            // Slopes, topography-following cell sizes, and gravity components at quarter of cell

                            if ( cdomain2[i] == 1 ) {

                                wbetax = betax[i];
                                wbetay = betay[i];
                                wbetaxy = betaxy[i];

                            } else {

                                wbetax = fbeta( wintelev[in[i][5]][j], wintelev[in[i][2]][j], 2.0, sico );
                                wbetay = fbeta( wintelev[in[i][4]][j], wintelev[in[i][1]][j], 2.0, sico );
                                wbetaxy = fbetaxy( wbetax, wbetay );

                            }

                            if ( sico.CORRHEIGHT != 0 ) {

                                wdx = sico.CSZ / cos( wbetax );
                                wdy = sico.CSZ / cos( wbetay );

                            } else {

                                wdx = sico.CSZ;
                                wdy = sico.CSZ;
                            }

                            wgrav[0] = sico.GRAVITY * cos( wbetax );
                            wgrav[1] = sico.GRAVITY * cos( wbetay );
                            wgrav[2] = sico.GRAVITY * sin( wbetax );
                            wgrav[3] = sico.GRAVITY * sin( wbetay );
                            wgrav[4] = sico.GRAVITY * cos( wbetaxy );


                            // Total flow depths

                            for ( l=1; l<9; l++ ) {

                                if ( sico.MODEL <= 3 ) wh[l] = wintb[in[i][l]][j][0];
                                else if ( sico.MODEL == 7 ) wh[l] = wintb[in[i][l]][j][0] + wintb[in[i][l]][j][3] + wintb[in[i][l]][j][6];
                            }


                            // Gradients of total flow depth

                            whx = ( wh[2] - wh[5] ) / ( 2 * dx[i] );
                            if ( sico.NOOSC > 0 && sico.NOOSC < 5 && ( wh[2] == 0 || wh[5] == 0 )) whx = tan( betax[i] ); // balance of forces at reservoir edge

                            why = ( wh[1] - wh[4] ) / ( 2 * dy[i] );
                            if ( sico.NOOSC > 0 && sico.NOOSC < 5 && ( wh[1] == 0 || wh[4] == 0 )) why = tan( betay[i] ); 


                            // Flow velocities in x and y directions

                            vflowx = fdiv( wintb[i][j][1], wintb[i][j][0], sico.HFLOWMIN );
                            vflowy = fdiv( wintb[i][j][2], wintb[i][j][0], sico.HFLOWMIN );

                            if ( sico.MODEL == 7 ) vflowx2 = fdiv( wintb[i][j][4], wintb[i][j][3], sico.HFLOWMIN ); else vflowx2 = 0;
                            if ( sico.MODEL == 7 ) vflowy2 = fdiv( wintb[i][j][5], wintb[i][j][3], sico.HFLOWMIN ); else vflowy2 = 0;

                            if ( sico.MODEL == 7 ) vflowx3 = fdiv( wintb[i][j][7], wintb[i][j][6], sico.HFLOWMIN ); else vflowx3 = 0;
                            if ( sico.MODEL == 7 ) vflowy3 = fdiv( wintb[i][j][8], wintb[i][j][6], sico.HFLOWMIN ); else vflowy3 = 0;


                            // Gradients of flow velocities

                            for ( l=1; l<9; l++ ) {

                                wu[l] = fdiv( wintb[in[i][l]][j][1], wintb[in[i][l]][j][0], sico.HFLOWMIN );
                                wv[l] = fdiv( wintb[in[i][l]][j][2], wintb[in[i][l]][j][0], sico.HFLOWMIN );

                                if ( sico.MODEL == 7 ) wu2[l] = fdiv( wintb[in[i][l]][j][4], wintb[in[i][l]][j][3], sico.HFLOWMIN ); else wu2[l] = 0;
                                if ( sico.MODEL == 7 ) wv2[l] = fdiv( wintb[in[i][l]][j][5], wintb[in[i][l]][j][3], sico.HFLOWMIN ); else wv2[l] = 0;

                                if ( sico.MODEL == 7 ) wu3[l] = fdiv( wintb[in[i][l]][j][7], wintb[in[i][l]][j][6], sico.HFLOWMIN ); else wu3[l] = 0;
                                if ( sico.MODEL == 7 ) wv3[l] = fdiv( wintb[in[i][l]][j][8], wintb[in[i][l]][j][6], sico.HFLOWMIN ); else wv3[l] = 0;

                                walpha[l] = fdiv( wintb[in[i][l]][j][0], wh[l], sico.HFLOWMIN );
                                if ( sico.MODEL == 7 ) walpha2[l] = fdiv( wintb[in[i][l]][j][3], wh[l], sico.HFLOWMIN ); else walpha2[l] = 0;
                                if ( sico.MODEL == 7 ) walpha3[l] = fdiv( wintb[in[i][l]][j][6], wh[l], sico.HFLOWMIN ); else walpha3[l] = 0;
                            }

                            wdu[0] = ( wu[2] - wu[5] ) / ( 2 * wdx );
                            wdv[0] = ( wv[1] - wv[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7  ) wdu[1] = ( wu2[2] - wu2[5] ) / ( 2 * wdx ); else wdu[1] = 0;
                            if ( sico.MODEL == 7  ) wdv[1] = ( wv2[1] - wv2[4] ) / ( 2 * wdy ); else wdv[1] = 0;

                            if ( sico.MODEL == 7 ) wdu[2] = ( wu3[2] - wu3[5] ) / ( 2 * wdx ); else wdu[2] = 0;
                            if ( sico.MODEL == 7 ) wdv[2] = ( wv3[1] - wv3[4] ) / ( 2 * wdy ); else wdv[2] = 0;

                            xwdv[0] = ( wv[2] - wv[5] ) / ( 2 * wdx );
                            xwdu[0] = ( wu[1] - wu[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) xwdv[1] = ( wv2[2] - wv2[5] ) / ( 2 * dx[i] ); else xwdv[1] = 0;
                            if ( sico.MODEL == 7 ) xwdu[1] = ( wu2[1] - wu2[4] ) / ( 2 * dy[i] ); else xwdu[1] = 0;

                            if ( sico.MODEL == 7 ) xwdv[2] = ( wv3[2] - wv3[5] ) / ( 2 * dx[i] ); else xwdv[2] = 0;
                            if ( sico.MODEL == 7 ) xwdu[2] = ( wu3[1] - wu3[4] ) / ( 2 * dy[i] ); else xwdu[2] = 0;

                            walphax = ( walpha[2] - walpha[5] ) / ( 2 * wdx );
                            walphay = ( walpha[1] - walpha[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) walphax2 = ( walpha2[2] - walpha2[5] ) / ( 2 * wdx ); else walphax2 = 0;
                            if ( sico.MODEL == 7 ) walphay2 = ( walpha2[1] - walpha2[4] ) / ( 2 * wdy ); else walphay2 = 0;

                            if ( sico.MODEL == 7 ) walphax3 = ( walpha3[2] - walpha3[5] ) / ( 2 * wdx ); else walphax3 = 0;
                            if ( sico.MODEL == 7 ) walphay3 = ( walpha3[1] - walpha3[4] ) / ( 2 * wdy ); else walphay3 = 0;

 
                            // Earth pressure coefficients

                            if ( sico.MODEL <= 3 ) { whflow = wintb[i][j][0], whflow2 = 0, whflow3 = 0; }
                            else if ( sico.MODEL == 7 ) { whflow = wintb[i][j][0], whflow2 = wintb[i][j][3], whflow3 = wintb[i][j][6]; }

                            hekin = 0.5 * whflow * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * whflow2 * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * whflow3 * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));

                            for ( ll=0; ll<sico.PMAX; ll++ ) {

                                if ( sico.PHASES[ll] < 2 ) {

                                    gkx[ll] = fk( whflow+whflow2+whflow3, wintb[i][j][3*ll], vflowx, vflowy, wdu[ll], hekin, ll, sico, sflow );
                                    gky[ll] = fk( whflow+whflow2+whflow3, wintb[i][j][3*ll], vflowx, vflowy, wdv[ll], hekin, ll, sico, sflow );
                                }
                            }


                            // Flux terms

                            kappau = fcurv( vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, welev, wgrav, sico );
                            vm = fvm( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, sico, sflow );
                            cdrag = fdrag( whflow, whflow2, whflow3, sico, sflow );
                            gze = fgze( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wh, whx, why, wdu, wdv, xwdu, xwdv, 
                                wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, wdx, wdy, wgrav, wbetax, wbetay, cdrag, tsum, kappau, sico, sflow );
                            disp = fdisp( whflow,whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wdu, wdv, xwdu, xwdv, 
                                wbetax, wbetay, gkx, gky, vm, cdrag, tsum, sico, sflow );

                            gf = ff( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                                gkx, gze, kappau, vm, whx, walphax, cdrag, wgrav, disp, wbetax, sflow, sico );
                            gg = fg( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                                gky, gze, kappau, vm, why, walphay, cdrag, wgrav, disp, wbetay, sflow, sico );

                            for ( k=0; k<sico.NVECTMIN; k++ ) { f[j][k] = gf[k]; g[j][k] = gg[k]; }
                            free( kappau ); free( vm ); free( cdrag ); free( gze ); free( disp ); free( gf ); free( gg );


                            // Flow velocities in x and y directions

                            vflowx = fdiv( wintc[i][j][1], wintc[i][j][0], sico.HFLOWMIN );
                            vflowy = fdiv( wintc[i][j][2], wintc[i][j][0], sico.HFLOWMIN );

                            if ( sico.MODEL == 7 ) vflowx2 = fdiv( wintc[i][j][4], wintc[i][j][3], sico.HFLOWMIN ); else vflowx2 = 0;
                            if ( sico.MODEL == 7 ) vflowy2 = fdiv( wintc[i][j][5], wintc[i][j][3], sico.HFLOWMIN ); else vflowy2 = 0;

                            if ( sico.MODEL == 7 ) vflowx3 = fdiv( wintc[i][j][7], wintc[i][j][6], sico.HFLOWMIN ); else vflowx3 = 0;
                            if ( sico.MODEL == 7 ) vflowy3 = fdiv( wintc[i][j][8], wintc[i][j][6], sico.HFLOWMIN ); else vflowy3 = 0;


                            // Total flow depths, velocities, fractions, and slopes of neighbouring cells

                            welev[0] = wintelev[i][0];

                            for ( l=0; l<9; l++ ) {

                                welev[l] = wintelev[in[i][l]][j];
                                
                                if ( cdomain2[in[i][l]] == 1 ) {

                                    nbetax[l] = betax[in[i][l]];
                                    nbetay[l] = betay[in[i][l]];

                                } else {

                                    nbetax[l] = fbeta( wintelev[in[in[i][l]][5]][j], wintelev[in[in[i][l]][2]][j], 2.0, sico );
                                    nbetay[l] = fbeta( wintelev[in[in[i][l]][4]][j], wintelev[in[in[i][l]][1]][j], 2.0, sico );

                                }

                                if ( sico.MODEL <= 3 ) wh[l] = wintc[in[i][l]][j][0];
                                else if ( sico.MODEL == 7 ) wh[l] = wintc[in[i][l]][j][0] + wintc[in[i][l]][j][3] + wintc[in[i][l]][j][6];

                                wu[l] = fdiv( wintc[in[i][l]][j][1], wintc[in[i][l]][j][0], sico.HFLOWMIN );
                                wv[l] = fdiv( wintc[in[i][l]][j][2], wintc[in[i][l]][j][0], sico.HFLOWMIN );

                                if ( sico.MODEL == 7 ) wu2[l] = fdiv( wintc[in[i][l]][j][4], wintc[in[i][l]][j][3], sico.HFLOWMIN ); else wu2[l] = 0;
                                if ( sico.MODEL == 7 ) wv2[l] = fdiv( wintc[in[i][l]][j][5], wintc[in[i][l]][j][3], sico.HFLOWMIN ); else wv2[l] = 0;

                                if ( sico.MODEL == 7 ) wu3[l] = fdiv( wintc[in[i][l]][j][7], wintc[in[i][l]][j][6], sico.HFLOWMIN ); else wu3[l] = 0;
                                if ( sico.MODEL == 7 ) wv3[l] = fdiv( wintc[in[i][l]][j][8], wintc[in[i][l]][j][6], sico.HFLOWMIN ); else wv3[l] = 0;

                                walpha[l] = fdiv( wintc[in[i][l]][j][0], wh[l], sico.HFLOWMIN );
                                if ( sico.MODEL == 7 ) walpha2[l] = fdiv( wintc[in[i][l]][j][3], wh[l], sico.HFLOWMIN ); else walpha2[l] = 0;
                                if ( sico.MODEL == 7 ) walpha3[l] = fdiv( wintc[in[i][l]][j][6], wh[l], sico.HFLOWMIN ); else walpha3[l] = 0;
                            }


                            // Gradients of total flow depth, flow velocities, and fractions

                            whx = ( wh[2] - wh[5] ) / ( 2 * wdx );
                            if ( sico.NOOSC > 0 && sico.NOOSC < 5 && ( wh[2] == 0 || wh[5] == 0 )) whx = tan( wbetax ); // balance of forces at reservoir edge

                            why = ( wh[1] - wh[4] ) / ( 2 * wdy );
                            if ( sico.NOOSC > 0 && sico.NOOSC < 5 && ( wh[1] == 0 || wh[4] == 0 )) why = tan( wbetay );

                            wdu[0] = ( wu[2] - wu[5] ) / ( 2 * wdx );
                            wdv[0] = ( wv[1] - wv[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) wdu[1] = ( wu2[2] - wu2[5] ) / ( 2 * wdx ); else wdu[1] = 0;
                            if ( sico.MODEL == 7 ) wdv[1] = ( wv2[1] - wv2[4] ) / ( 2 * wdy ); else wdv[1] = 0;

                            if ( sico.MODEL == 7 ) wdu[2] = ( wu3[2] - wu3[5] ) / ( 2 * wdx ); else wdu[2] = 0;
                            if ( sico.MODEL == 7 ) wdv[2] = ( wv3[1] - wv3[4] ) / ( 2 * wdy ); else wdv[2] = 0;

                            xwdv[0] = ( wv[2] - wv[5] ) / ( 2 * wdx );
                            xwdu[0] = ( wu[1] - wu[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) xwdv[1] = ( wv2[2] - wv2[5] ) / ( 2 * dx[i] ); else xwdv[1] = 0;
                            if ( sico.MODEL == 7 ) xwdu[1] = ( wu2[1] - wu2[4] ) / ( 2 * dy[i] ); else xwdu[1] = 0;

                            if ( sico.MODEL == 7 ) xwdv[2] = ( wv3[2] - wv3[5] ) / ( 2 * dx[i] ); else xwdv[2] = 0;
                            if ( sico.MODEL == 7 ) xwdu[2] = ( wu3[1] - wu3[4] ) / ( 2 * dy[i] ); else xwdu[2] = 0;

                            walphax = ( walpha[2] - walpha[5] ) / ( 2 * wdx );
                            walphay = ( walpha[1] - walpha[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) walphax2 = ( walpha2[2] - walpha2[5] ) / ( 2 * wdx ); else walphax2 = 0;
                            if ( sico.MODEL == 7 ) walphay2 = ( walpha2[1] - walpha2[4] ) / ( 2 * wdy ); else walphay2 = 0;

                            if ( sico.MODEL == 7 ) walphax3 = ( walpha3[2] - walpha3[5] ) / ( 2 * wdx ); else walphax3 = 0;
                            if ( sico.MODEL == 7 ) walphay3 = ( walpha3[1] - walpha3[4] ) / ( 2 * wdy ); else walphay3 = 0;


                            // Source terms

                            if ( sico.MODEL <= 3 ) { whflow = wintc[i][j][0], whflow2 = 0, whflow3 = 0; }
                            else if ( sico.MODEL == 7 ) { whflow = wintc[i][j][0], whflow2 = wintc[i][j][3], whflow3 = wintc[i][j][6]; }

                            hekin = 0.5 * whflow * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * whflow2 * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * whflow3 * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));

                            kappau = fcurv( vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, welev, wgrav, sico );
                            vm = fvm( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, sico, sflow );
                            cdrag = fdrag( whflow, whflow2, whflow3, sico, sflow );

                            gze = fgze( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wh, whx, why, wdu, wdv, xwdu, xwdv, 
                                wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, wdx, wdy, wgrav, wbetax, wbetay, cdrag, tsum, kappau, sico, sflow );

                            gs = fs( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                                whx, why, walphax, walphay, walphax2, walphay2, walphax3, walphay3, wgrav, gze, wbetax, wbetay, kappau, cdrag, sico, sflow );

                            for ( k=0; k<sico.NVECTMIN; k++ ) s[j][k] = gs[k];
                            free( gs );


                            // Deceleration terms

                            for ( k=0; k<18; k++ ) wwd[k] = d[i][j][k+sico.NVECTMIN];

                            gdecel = fd( wh, whx, why, wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, nbetax, nbetay, gze, wdx, wdy, 
                                wwd, kappau, hekin, wgrav, sico, sflow );

                            for ( k=0; k<sico.NVECTMIN+18; k++ ) d[i][j][k] = gdecel[k];
                            free( kappau ); free( cdrag ); free( vm ); free( gze ); free( gdecel );
                        }


                        // Value of vector at top right corner of the cell

                        for ( k=0; k<sico.NVECTMIN; k++ ) {

                            wintd[i][k]=0.25*(winta[i][0][k]+winta[i][1][k]+winta[i][2][k]+winta[i][3][k])
                                -tlength/dx[i]*(0.5*(f[2][k]+f[3][k])-0.5*(f[0][k]+f[1][k]))
                                -tlength/dy[i]*(0.5*(g[1][k]+g[3][k])-0.5*(g[0][k]+g[2][k]))
                                +0.25*tlength*(s[0][k]+s[1][k]+s[2][k]+s[3][k]);

                            if ( sico.NOOSC > 0 && sico.NOOSC < 5 && fabs(-tlength/dx[i]*(0.5*(f[2][k]+f[3][k])-0.5*(f[0][k]+f[1][k]))
                                -tlength/dy[i]*(0.5*(g[1][k]+g[3][k])-0.5*(g[0][k]+g[2][k]))
                                 +0.25*tlength*(s[0][k]+s[1][k]+s[2][k]+s[3][k])) < pow( 10, -7 )) { // balance of forces at reservoir edge

                                 if ( p == 0 ) wintd[i][k] = aw[in[i][0]][k];
                                 else wintd[i][k] = aw[in[i][3]][k];
                            }

                            wintdtest[i][k]=0.25*(winta[i][0][k]+winta[i][1][k]+winta[i][2][k]+winta[i][3][k])
                                -tlength/dx[i]*(0.5*(f[2][k]+f[3][k])-0.5*(f[0][k]+f[1][k]))
                                -tlength/dy[i]*(0.5*(g[1][k]+g[3][k])-0.5*(g[0][k]+g[2][k]))
                                +0.25*tlength*(s[0][k]+s[1][k]+s[2][k]+s[3][k]-d[i][0][k]-d[i][1][k]-d[i][2][k]-d[i][3][k]);

                            if ( sico.NOOSC > 0 && sico.NOOSC < 5 && fabs(-tlength/dx[i]*(0.5*(f[2][k]+f[3][k])-0.5*(f[0][k]+f[1][k]))
                                -tlength/dy[i]*(0.5*(g[1][k]+g[3][k])-0.5*(g[0][k]+g[2][k]))
                                +0.25*tlength*(s[0][k]+s[1][k]+s[2][k]+s[3][k]-d[i][0][k]-d[i][1][k]-d[i][2][k]-d[i][3][k])) < pow( 10, -7 )) {

                                if ( p == 0 ) wintdtest[i][k] = aw[in[i][0]][k];
                                else wintdtest[i][k] = aw[in[i][3]][k];
                            }
                        }

                    } else for ( k=0; k<sico.NVECTMIN; k++ ) { wintd[i][k] = aw[i][k]; wintdtest[i][k] = aw[i][k]; }
                }


                // Moving vector if second sub-timestep and writing values to temporary vector

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( pxslide[i] == 0 ) {

                        ctrlv = 1; if ( in[i][6] == 0 ) ctrlv = 0;
                        if ( ctrlv == 1 ) {
                            ctrlr = 0;
                            if ( sico.MODEL <= 3 ) {

                                if ( wintd[i][0] > 0  ) ctrlr = 1;
                                if ( wintd[in[i][6]][0] > 0 ) ctrlr = 1;

                            } else if ( sico.MODEL == 7 ) {

                                if ( wintd[i][0] > 0 || wintd[i][3] > 0 || wintd[i][6] > 0 ) ctrlr = 1;
                                if ( wintd[in[i][6]][0] > 0 || wintd[in[i][6]][3] > 0 || wintd[in[i][6]][6] > 0 ) ctrlr = 1;
                            }
                        }

                        if ( ctrlv == 1 && ctrlr == 1 && cdomain[i] > 3 ) {

                            for ( k=0; k<sico.NVECTMIN; k++ ) {
                                if ( p == 1 && ( cdomain[in[i][6]] > 3 || cflux[in[i][6]] == 0 )) {

                                    awt[i][k]=wintd[in[i][6]][k];
                                    awttest[k]=wintdtest[in[i][6]][k];
                                }

                                else if ( p != 1 && ( cdomain[i] > 3 || cflux[i] == 0 )) {

                                    awt[i][k]=wintd[i][k];
                                    awttest[k]=wintdtest[i][k];
                                }
                                else if ( k<sico.NVECTMIN ) { awt[i][k] = 0; awttest[k] = 0; }
                                else { awt[i][k] = aw[i][k]; awttest[k] = aw[i][k]; }
                                if (( fsign(awt[i][k]) == fsign(awttest[k]) && fabs(awt[i][k]) > fabs(awttest[k])) || sico.ORIGINAL == 1 ) awt[i][k] = awttest[k];
                                else if ( fsign(awt[i][k]) == fsign(awttest[k])) awt[i][k] = awt[i][k];
                                else if ( k==0 || k==3 || k==6 ) awt[i][k] = aw[i][k];
                                else awt[i][k] = 0;
                            }
                        }
                        else if ( cflux[i] != 0 && cstopped[i] != 1 ) { for ( k=0; k<sico.NVECTMIN; k++ ) awt[i][k] = 0; }
                        else { for ( k=0; k<sico.NVECTMIN; k++ ) awt[i][k] = aw[i][k]; }

                        if ( awt[i][0] < 0 ) { awt[i][0] = 0; awt[i][1] = 0; awt[i][2] = 0; }
                        if ( sico.MODEL == 7 && awt[i][3] < 0 ) { awt[i][3] = 0; awt[i][4] = 0; awt[i][5] = 0; }
                        if ( sico.MODEL == 7 && awt[i][6] < 0 ) { awt[i][6] = 0; awt[i][7] = 0; awt[i][8] = 0; }
                        for ( k=0; k<sico.NVECTMIN; k++ ) { if ( isnan( awt[i][k] ) != 0 ) awt[i][k] = 0; }

                        if ( cdomain2[i] == 1 ) { for ( k=0; k<3; k++ ) awt[i][k] = aw[i][k]; }
                        if ( sico.MODEL == 7 && cdomain2[i] == 1 ) { for ( k=3; k<6; k++ ) awt[i][k] = aw[i][k]; }
                        if ( sico.MODEL == 7 && cdomain2[i] == 1 ) { for ( k=6; k<9; k++ ) awt[i][k] = aw[i][k]; }
                    }
                  }
                }
                

// -- STOP --- Flow propagation with NOC scheme -----------------------------------------------------------------


// -- START -- Diffusion control --------------------------------------------------------------------------------


                // *** Experimental, may yield unplausible results under certain conditions
                // *** cedge[i][j] (cedge2[i][j], cedge3[i][j]): edge cell with regard to mixture or PHASE 1 (PHASE 2, 3) in direction j - 1 = yes, 2 = no
                // *** cready[i][j] (cready2[i][j], cready3[i][j]): degree of mixture or PHASE 1 (PHASE 2, 3) fill of neighbour cells in direction j - ratio
                // *** cedge0[i] (cedge02[i], cedge03[i]): number of cells which may provide mixture or PHASE 1 (PHASE 2, 3) inflow
                // *** cneighbours[i] (cneighbours2[i], cneighbours3[i]): number of neighbour cells with non-zero flow depth


                if ( sico.CURVCTRL > 2 && ctrl_noosc == 0 ) {

                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];

                        cedge0[i] = 0; cedge02[i] = 0; cedge03[i] = 0;
                    }


                // Updating degree of fill of neighbour cells

                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];
                        if ( p == 1 ) i2 = i;
                        else i2 = i;

                        if ( fvalid( aw[i2][0], aw[i2][3], aw[i2][6], sico ) == 1 ) { // if cell was flow cell at the beginning of the time step

                            if ( aw[i2][0] > 0 ) {
                                // identifying directions for which cell is edge cell with regard to mixture or PHASE 1

                                for ( j=1; j<9; j++ ) {
                                    if ( aw[in[i2][j]][0] == 0  ) cedge[i][j] = 1;
                                    else cedge[i][j] = 0;
                                }
                            }
                            if ( sico.MODEL == 7 && aw[i2][3] > 0 ) {
                                // identifying directions for which cell is edge cell with regard to PHASE 2

                                for ( j=1; j<9; j++ ) {
                                    if ( aw[in[i2][j]][3] == 0 ) cedge2[i][j] = 1;
                                    else cedge2[i][j] = 0;
                                }
                            }
                            if ( sico.MODEL == 7 && aw[i2][6] > 0 ) {
                                // identifying directions for which cell is edge cell with regard to PHASE 3

                                for ( j=1; j<9; j++ ) {
                                    if ( aw[in[i2][j]][6] == 0 ) cedge3[i][j] = 1;
                                    else cedge3[i][j] = 0;
                                }
                            }

                            if ( nsum == 1 ) { // for first time step:

                                if ( cedge[i][2] == 1 && aw[i2][0] > 0 ) cready[i][2] = 2; // degree of fill for mixture or PHASE 1
                                if ( cedge[i][5] == 1 && aw[i2][0] > 0 ) cready[i][5] = 2;
                                if ( cedge[i][1] == 1 && aw[i2][0] > 0 ) cready[i][1] = 2;
                                if ( cedge[i][4] == 1 && aw[i2][0] > 0 ) cready[i][4] = 2;

                                if ( sico.MODEL == 7 && cedge2[i][2] == 1 && aw[i2][3] > 0 ) // degree of fill for PHASE 2
                                    cready2[i][2] = 2;
                                if ( sico.MODEL == 7 && cedge2[i][5] == 1 && aw[i2][3] > 0 )
                                    cready2[i][5] = 2;
                                if ( sico.MODEL == 7 && cedge2[i][1] == 1 && aw[i2][3] > 0 )
                                    cready2[i][1] = 2;
                                if ( sico.MODEL == 7 && cedge2[i][4] == 1 && aw[i2][3] > 0 )
                                    cready2[i][4] = 2;

                                if ( sico.MODEL == 7 && cedge3[i][2] == 1 && aw[i2][6] > 0 ) // degree of fill for PHASE 3
                                    cready3[i][2] = 2;
                                if ( sico.MODEL == 7 && cedge3[i][5] == 1 && aw[i2][6] > 0 )
                                    cready3[i][5] = 2;
                                if ( sico.MODEL == 7 && cedge3[i][1] == 1 && aw[i2][6] > 0 )
                                    cready3[i][1] = 2;
                                if ( sico.MODEL == 7 && cedge3[i][4] == 1 && aw[i2][6] > 0 )
                                    cready3[i][4] = 2;

                            } else { // for all other time steps:

                                if ( cedge[i][2] == 1 && aw[i2][0] > 0 )
                                    cready[i][2] = cready[i2][2] + aw[i2][1] / aw[i2][0] * 2*tlength / sico.CSZ; // degree of fill for mixture or PHASE 1
                                if ( cedge[i][5] == 1 && aw[i2][0] > 0 )
                                    cready[i][5] = cready[i2][5] - aw[i2][1] / aw[i2][0] * 2*tlength / sico.CSZ;
                                if ( cedge[i][1] == 1 && aw[i2][0] > 0 )
                                    cready[i][1] = cready[i2][1] + aw[i2][2] / aw[i2][0] * 2*tlength / sico.CSZ;
                                if ( cedge[i][4] == 1 && aw[i2][0] > 0 )
                                    cready[i][4] = cready[i2][4] - aw[i2][2] / aw[i2][0] * 2*tlength / sico.CSZ;

                                if ( sico.MODEL == 7 && cedge2[i][2] == 1 && aw[i2][3] > 0 ) // degree of fill for PHASE 2
                                    cready2[i][2] = cready2[i2][2] + aw[i2][4] / aw[i2][3] * 2*tlength / sico.CSZ;
                                if ( sico.MODEL == 7 && cedge2[i][5] == 1 && aw[i2][3] > 0 )
                                    cready2[i][5] = cready2[i2][5] - aw[i2][4] / aw[i2][3] * 2*tlength / sico.CSZ;
                                if ( sico.MODEL == 7 && cedge2[i][1] == 1 && aw[i2][3] > 0 )
                                    cready2[i][1] = cready2[i2][1] + aw[i2][5] / aw[i2][3] * 2*tlength / sico.CSZ;
                                if ( sico.MODEL == 7 && cedge2[i][4] == 1 && aw[i2][3] > 0 )
                                    cready2[i][4] = cready2[i2][4] - aw[i2][5] / aw[i2][3] * 2*tlength / sico.CSZ;

                                if ( sico.MODEL == 7 && cedge3[i][2] == 1 && aw[i2][6] > 0 ) // degree of fill for PHASE 3
                                    cready3[i][2] = cready3[i2][2] + aw[i2][7] / aw[i2][6] * 2*tlength / sico.CSZ;
                                if ( sico.MODEL == 7 && cedge3[i][5] == 1 && aw[i2][6] > 0 )
                                    cready3[i][5] = cready3[i2][5] - aw[i2][7] / aw[i2][6] * 2*tlength / sico.CSZ;
                                if ( sico.MODEL == 7 && cedge3[i][1] == 1 && aw[i2][6] > 0 )
                                    cready3[i][1] = cready3[i2][1] + aw[i2][8] / aw[i2][6] * 2*tlength / sico.CSZ;
                                if ( sico.MODEL == 7 && cedge3[i][4] == 1 && aw[i2][6] > 0 )
                                    cready3[i][4] = cready3[i2][4] - aw[i2][8] / aw[i2][6] * 2*tlength / sico.CSZ;
                            }
                        }
                    }


                    // Number of cells which may provide inflow to neighbour cells (cedge0, cedge02, cedge03 )

                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];
                        if ( p == 1 ) i2 = i;
                        else i2 = i;

                        cneighbours[i] = 0; cneighbours2[i] = 0; cneighbours3[i] = 0;

                        if ( fvalid( aw[i2][0], aw[i2][3], aw[i2][6], sico ) == 1 ) {

                            for ( j=1; j<3; j++ ) {

                                if (( cedge [i][j] == 0 || cready [i][j] > 1 )
                                    && aw[i2][0] > 0 ) cedge0 [in[i][j]] += 1;
                                if ( sico.MODEL == 7 && ( cedge2[i][j] == 0 || cready2[i][j] > 1 )
                                    && aw[i2][3] > 0 ) cedge02[in[i][j]] += 1;
                                if ( sico.MODEL == 7 && ( cedge3[i][j] == 0 || cready3[i][j] > 1 )
                                    && aw[i2][6] > 0 ) cedge03[in[i][j]] += 1;
                            }
                            for ( j=4; j<6; j++ ) {

                                if (( cedge [i][j] == 0 || cready [i][j] > 1 )
                                    && aw[i2][0] > 0 ) cedge0 [in[i][j]] += 1;
                                if ( sico.MODEL == 7 && ( cedge2[i][j] == 0 || cready2[i][j] > 1 )
                                    && aw[i2][3] > 0 ) cedge02[in[i][j]] += 1;
                                if ( sico.MODEL == 7 && ( cedge3[i][j] == 0 || cready3[i][j] > 1 )
                                    && aw[i2][6] > 0 ) cedge03[in[i][j]] += 1;
                            }
                            if (( cedge[i][3] == 0 || cready[i][1] > 1 || cready[i][2] > 1 )
                                && aw[i2][0] > 0 ) cedge0 [in[i][3]] += 1;
                            if (( cedge[i][6] == 0 || cready[i][4] > 1 || cready[i][5] > 1 )
                                && aw[i2][0] > 0 ) cedge0 [in[i][6]] += 1;
                            if (( cedge[i][7] == 0 || cready[i][2] > 1 || cready[i][4] > 1 )
                                && aw[i2][0] > 0 ) cedge0 [in[i][7]] += 1;
                            if (( cedge[i][8] == 0 || cready[i][1] > 1 || cready[i][5] > 1 )
                                && aw[i2][0] > 0 ) cedge0 [in[i][8]] += 1;

                            if ( sico.MODEL == 7 && ( cedge2[i][3] == 0 || cready2[i][1] > 1
                                || cready2[i][2] > 1 ) && aw[i2][3] > 0 ) cedge02[in[i][3]] += 1;
                            if ( sico.MODEL == 7 && ( cedge2[i][6] == 0 || cready2[i][4] > 1
                                || cready2[i][5] > 1 ) && aw[i2][3] > 0 ) cedge02[in[i][6]] += 1;
                            if ( sico.MODEL == 7 && ( cedge2[i][7] == 0 || cready2[i][2] > 1
                                || cready2[i][4] > 1 ) && aw[i2][3] > 0 ) cedge02[in[i][7]] += 1;
                            if ( sico.MODEL == 7 && ( cedge2[i][8] == 0 || cready2[i][1] > 1
                                || cready2[i][5] > 1 ) && aw[i2][3] > 0 ) cedge02[in[i][8]] += 1;

                            if ( sico.MODEL == 7 && ( cedge3[i][3] == 0 || cready3[i][1] > 1
                                || cready3[i][2] > 1 ) && aw[i2][6] > 0 ) cedge03[in[i][3]] += 1;
                            if ( sico.MODEL == 7 && ( cedge3[i][6] == 0 || cready3[i][4] > 1
                                || cready3[i][5] > 1 ) && aw[i2][6] > 0 ) cedge03[in[i][6]] += 1;
                            if ( sico.MODEL == 7 && ( cedge3[i][7] == 0 || cready3[i][2] > 1
                                || cready3[i][4] > 1 ) && aw[i2][6] > 0 ) cedge03[in[i][7]] += 1;
                            if ( sico.MODEL == 7 && ( cedge3[i][8] == 0 || cready3[i][1] > 1
                                || cready3[i][5] > 1 ) && aw[i2][6] > 0 ) cedge03[in[i][8]] += 1;
                        }
                    }


                    // Number of neighbour cells with non-zero flow depth (cneighbours, cneighbours2, cneighbours3)

                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];

                        if ( cdomain[i] >= 3 ) {

                            for ( j=1; j<9; j++ ) { 

                                if ( cedge [i][j] == 1 && cedge0[in[i][j]] == 0 ) {

                                    if ( awt[in[i][j]][0] > 0 ) cneighbours[in[i][j]] += 1;
                                }
                            }
                            for ( j=1; j<9; j++ ) {

                                if ( cedge2 [i][j] == 1 && cedge02[in[i][j]] == 0 ) { 

                                    if ( sico.MODEL == 7 && awt[in[i][j]][3] > 0 ) cneighbours2[in[i][j]] += 1;
                                }
                            }
                            for ( j=1; j<9; j++ ) {

                                if ( cedge3 [i][j] == 1 && cedge03[in[i][j]] == 0 ) { 

                                    if ( sico.MODEL == 7 && awt[in[i][j]][6] > 0 ) cneighbours3[in[i][j]] += 1;
                                }
                            }
                        }
                    }


                    // Reallocating flow depth and momentum from outside cells to edge cells

                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];

                        if ( cdomain[i] >= 3 && pxslide[i] == 0 ) {

                            for ( j=1; j<9; j++ ) {

                                if ( cedge [i][j] == 1 && cedge0[in[i][j]] == 0 ) {

                                    for ( k=0; k<3; k++ ) {
                                        if ( cneighbours[in[i][j]] > 0 )
                                            awt[i][k] += awt[in[i][j]][k] / cneighbours[in[i][j]];
                                    }
                                }
                                if ( sico.MODEL == 7 && cedge2[i][j] == 1 && cedge02[in[i][j]] == 0 ) {

                                    for ( k=3; k<6; k++ ) {
                                        if ( cneighbours2[in[i][j]] > 0 )
                                            awt[i][k] += awt[in[i][j]][k] / cneighbours2[in[i][j]];
                                    }
                                }
                                if ( sico.MODEL == 7 && cedge3[i][j] == 1 && cedge03[in[i][j]] == 0 ) {

                                    for ( k=6; k<9; k++ ) {
                                        if ( cneighbours3[in[i][j]] > 0 )
                                            awt[i][k] += awt[in[i][j]][k] / cneighbours3[in[i][j]];
                                    }
                                }
                            }
                        }
                    }


                    // Setting flow depth and momentum at outside cells to zero

                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];

                        if ( cdomain[i] >= 3 && pxslide[i] == 0 ) {

                            for ( j=1; j<9; j++ ) {

                                if ( cedge [i][j] == 1 && cedge0 [in[i][j]] == 0 ) {
                                    for ( k=0; k<3; k++ ) awt[in[i][j]][k] = 0;
                                }
                                if ( sico.MODEL == 7 && cedge2[i][j] == 1 && cedge02[in[i][j]] == 0 ) {
                                    for ( k=3; k<6; k++ ) awt[in[i][j]][k] = 0;
                                }
                                if ( sico.MODEL == 7 && cedge3[i][j] == 1 && cedge03[in[i][j]] == 0 ) {
                                    for ( k=6; k<9; k++ ) awt[in[i][j]][k] = 0;
                                }
                            }
                        }
                    }
                }

 
// -- STOP --- Diffusion control --------------------------------------------------------------------------------


// -- START -- Evaluating time step length and validity of time step --------------------------------------------


                vcelr = 0;

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( fvalid( awt[i][0], awt[i][3], awt[i][6], sico ) == 1 && pxslide[i] == 0 ) {

                        for ( l=1; l<9; l++ ) {

                            wu[l] = fdiv( awt[in[i][l]][1], awt[in[i][l]][0], sico.HFLOWMIN ); // flow velocities of neighbour cells
                            wv[l] = fdiv( awt[in[i][l]][2], awt[in[i][l]][0], sico.HFLOWMIN );

                            if ( sico.MODEL == 7 ) wu2[l] = fdiv( awt[in[i][l]][4], awt[in[i][l]][3], sico.HFLOWMIN ); else wu2[l] = 0;
                            if ( sico.MODEL == 7 ) wv2[l] = fdiv( awt[in[i][l]][5], awt[in[i][l]][3], sico.HFLOWMIN ); else wv2[l] = 0;

                            if ( sico.MODEL == 7 ) wu3[l] = fdiv( awt[in[i][l]][7], awt[in[i][l]][6], sico.HFLOWMIN ); else wu3[l] = 0;
                            if ( sico.MODEL == 7 ) wv3[l] = fdiv( awt[in[i][l]][8], awt[in[i][l]][6], sico.HFLOWMIN ); else wv3[l] = 0;
                        }

                        wdu[0] = wu[2] - wu[5]; // spatial changes of flow velocities
                        wdv[0] = wv[1] - wv[4];

                        if ( sico.MODEL == 7 ) { wdu[1] = wu2[2] - wu2[5]; wdv[1] = wv2[1] - wv2[4]; wdu[2] = wu3[2] - wu3[5]; wdv[2] = wv3[1] - wv3[4]; }

                        vflowx = fdiv( awt[i][1], awt[i][0], sico.HFLOWMIN ); // flow velocities at cell
                        vflowy = fdiv( awt[i][2], awt[i][0], sico.HFLOWMIN );

                        if ( sico.MODEL == 7 ) { 
                        
                            vflowx2 = fdiv( awt[i][4], awt[i][3], sico.HFLOWMIN );
                            vflowy2 = fdiv( awt[i][5], awt[i][3], sico.HFLOWMIN );
                            vflowx3 = fdiv( awt[i][7], awt[i][6], sico.HFLOWMIN );
                            vflowy3 = fdiv( awt[i][8], awt[i][6], sico.HFLOWMIN );
                            
                        } else { vflowx2 = 0; vflowy2 = 0; vflowx3 = 0;  vflowy3 = 0; }

                        if ( sico.MODEL <= 3 ) { whflow = aw[i][0], whflow2 = 0, whflow3 = 0; } // flow heights
                        else if ( sico.MODEL == 7 ) { whflow = aw[i][0], whflow2 = aw[i][3], whflow3 = aw[i][6]; }

                        hekin = 0.5 * whflow * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 )); // flow kinetic energy
                        if ( sico.MODEL == 7 ) hekin += 0.5 * whflow2 * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                        if ( sico.MODEL == 7 ) hekin += 0.5 * whflow3 * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));

                        for ( l=0; l<sico.PMAX; l++ ) {

                            if ( sico.PHASES[l] < 2 ) {

                                gkx[l] = fk( whflow+whflow2+whflow3, aw[i][3*l], vflowx, vflowy, wdu[l], hekin, l, sico, sflow ); // earth pressure coefficients
                                gky[l] = fk( whflow+whflow2+whflow3, aw[i][3*l], vflowx, vflowy, wdv[l], hekin, l, sico, sflow );
                            }
                        }

                        vcelr0 = fcvelr2( whflow, whflow2, whflow3, vflowx, vflowy, vflowx2, vflowy2, vflowx3, vflowy3, sflow, gkx, gky, betaxy, i, sico ); // flow velocity plus wave speed
                        if ( vcelr0 > vcelr ) vcelr = vcelr0;
                    }
                }

                cfl = vcelr * tlength / sico.CSZ; // updating CFL value

                if ( cfl <= sico.CFL[0] || vcelr == 0 ) { // if time step length meets CFL criterion

                    ccfl = 1; // setting control for fulfilment of CFL criterion to positive
                    if ( cfl > cflmax ) cflmax = cfl; // updating maximum CFL value, if necessary

                    tlengthpre = tlength; // updating length of previous time step

                    tsum += tlength; // updating time since release
                    tint += tlength; // updating time since last output

                    if ( tlength <= pow( 10, -7 )) { ccontinue = 0; csuccess = 0; } // defining numerical failure

                    if ( vcelr > 0 ) tlength0 = sico.CFL[0] * 0.90 * sico.CSZ / vcelr; // updating length of time step according to CFL criterion
                    else  tlength0 = sico.CFL[1]; 

                    vcelr = 0; // resetting flow velocity plus wave speed
                    
                } else { // if time step is too long too meet CFL criterion, repeating time step with reduced length

                    tlength = 0.75 * tlength * sico.CFL[0] / cfl; // defining new time step length
                    vcelr = 0; // resetting flow velocity plus wave speed
                }


// -- STOP --- Evaluating time step length and validity of time step --------------------------------------------


                #ifdef WITHGRASS


                    printf("   %i\t%i\t%.3f\t%.1f  \t%.1f\t...\r", nout, nsum, cflmax, 100*tlength, tsum ); // progress
                    fflush(stdout); // forcing immediate display


                #endif


// *** End of loop over time step lengths -----------------------------------------------------------------------


            }

            if ( cslide == 0 ) iloop = ib[0]; else iloop = sico.IMAX;
            for ( ix=0; ix<iloop; ix++ ) {

                if ( cslide == 0 ) i = ibasket[0][ix]; else i = ix;
                for ( k=0; k<sico.NVECTMIN; k++ ) aw[i][k] = awt[i][k]; // applying values of temporary vectors


// -- START -- Phase transformations ----------------------------------------------------------------------------
 

                if ( sico.MODEL == 7 ) {

                    if ( sico.TRANSSSFS == 1 ) { // PHASE 1 - PHASE 2 transformation coefficient
                        if ( ptransssfs[i] != sico.UNDEF ) sflow.TRANSSSFS = ptransssfs[i];
                        else sflow.TRANSSSFS = sflow.TRANSSSFS0;
                    } else sflow.TRANSSSFS = sflow.TRANSSSFS0;

                    if ( sico.TRANSSSFF == 1 ) { // PHASE 1 - PHASE 3 transformation coefficient
                        if ( ptransssff[i] != sico.UNDEF ) sflow.TRANSSSFF = ptransssff[i];
                        else sflow.TRANSSSFF = sflow.TRANSSSFF0;
                    } else sflow.TRANSSSFF = sflow.TRANSSSFF0;

                    if ( sico.TRANSFSFF == 1 ) { // PHASE 2 - PHASE 3 transformation coefficient
                        if ( ptransfsff[i] != sico.UNDEF ) sflow.TRANSFSFF = ptransfsff[i];
                        else sflow.TRANSFSFF = sflow.TRANSFSFF0;
                    } else sflow.TRANSFSFF = sflow.TRANSFSFF0;

                    if ( transformograph == 1 ) { // reading transformograph, if available

                        for ( trak=0; trak<=tratmax; trak++ ) {
                            if ( tratra[trak][0] <= tsum ) tratx=trak;
                            else break;
                        }

                        sflow.TRANSSSFS = tratra[tratx][1];
                        sflow.TRANSSSFF = tratra[tratx][2];
                        sflow.TRANSFSFF = tratra[tratx][3];
                    }

                    if ( aw[i][0] > 0 ) { vflowx = aw[i][1] / aw[i][0]; vflowy = aw[i][2] / aw[i][0]; } else { vflowx = 0; vflowy = 0; } // flow velocities
                    if ( aw[i][3] > 0 ) { vflowx2 = aw[i][4] / aw[i][3]; vflowy2 = aw[i][5] / aw[i][3]; } else { vflowx2 = 0; vflowy2 = 0; }
                    if ( aw[i][6] > 0 ) { vflowx3 = aw[i][7] / aw[i][6]; vflowy3 = aw[i][8] / aw[i][6]; } else { vflowx3 = 0; vflowy3 = 0; }

                    if ( transformograph == 1 ) ekin = 1; // transformograph does not use flow kinetic energy
                    else ekin = ( aw[i][0] * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 ))
                        + aw[i][3] * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ))
                        + aw[i][6] * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ))) * 0.5; // flow kinetic energy

                    if ( sflow.TRANSSSFF > 0 ) { // PHASE 1 to PHASE 3 transformation

                        trans = ffmin( aw[i][0], tlength * pow( 10, -fabs( sflow.TRANSSSFF )) * ekin );

                        aw[i][0] -= trans;
                        aw[i][1] -= trans * vflowx;
                        aw[i][2] -= trans * vflowy;

                        aw[i][6] += trans * sflow.RHO1 / sflow.RHO3;
                        aw[i][7] += trans * sflow.RHO1 / sflow.RHO3 * vflowx;
                        aw[i][8] += trans * sflow.RHO1 / sflow.RHO3 * vflowy;
                    }

                    if ( sflow.TRANSSSFF < 0 ) { // PHASE 3 to PHASE 1 transformation

                        trans = ffmax( -aw[i][6], -tlength * pow( 10, -fabs( sflow.TRANSSSFF )) * ekin );

                        aw[i][0] -= trans * sflow.RHO3 / sflow.RHO1;
                        aw[i][1] -= trans * vflowx3 * sflow.RHO3 / sflow.RHO1;
                        aw[i][2] -= trans * vflowy3 * sflow.RHO3 / sflow.RHO1;

                        aw[i][6] += trans;
                        aw[i][7] += trans * vflowx3;
                        aw[i][8] += trans * vflowy3;
                    }

                    if ( sflow.TRANSSSFS > 0 ) { // PHASE 1 to PHASE 2 transformation

                        trans = ffmin( aw[i][0], tlength * pow( 10, -fabs( sflow.TRANSSSFS )) * ekin );

                        aw[i][0] -= trans;
                        aw[i][1] -= trans * vflowx;
                        aw[i][2] -= trans * vflowy;

                        aw[i][3] += trans * sflow.RHO1 / sflow.RHO2;
                        aw[i][4] += trans * sflow.RHO1 / sflow.RHO2 * vflowx;
                        aw[i][5] += trans * sflow.RHO1 / sflow.RHO2 * vflowy;
                    }

                    if ( sflow.TRANSSSFS < 0 ) { // PHASE 2 to PHASE 1 transformation

                        trans = ffmax( -aw[i][3], -tlength * pow( 10, -fabs( sflow.TRANSSSFS )) * ekin );

                        aw[i][0] -= trans * sflow.RHO2 / sflow.RHO1;
                        aw[i][1] -= trans * vflowx2 * sflow.RHO2 / sflow.RHO1;
                        aw[i][2] -= trans * vflowy2 * sflow.RHO2 / sflow.RHO1;

                        aw[i][3] += trans;
                        aw[i][4] += trans * vflowx2;
                        aw[i][5] += trans * vflowy2;
                    }         

                    if ( sflow.TRANSFSFF > 0 ) { // PHASE 2 to PHASE 3 transformation

                        trans = ffmin( aw[i][3], tlength * pow( 10, -fabs( sflow.TRANSFSFF )) * ekin );

                        aw[i][3] -= trans;
                        aw[i][4] -= trans * vflowx2;
                        aw[i][5] -= trans * vflowy2;

                        aw[i][6] += trans * sflow.RHO2 / sflow.RHO3;
                        aw[i][7] += trans * sflow.RHO2 / sflow.RHO3 * vflowx2;
                        aw[i][8] += trans * sflow.RHO2 / sflow.RHO3 * vflowy2;
                    }

                    if ( sflow.TRANSFSFF < 0 ) { // PHASE 3 to PHASE 2 transformation

                        trans = ffmax( -aw[i][6], -tlength * pow( 10, -fabs( sflow.TRANSFSFF )) * ekin );

                        aw[i][3] -= trans * sflow.RHO3 / sflow.RHO2;
                        aw[i][4] -= trans * vflowx3 * sflow.RHO3 / sflow.RHO2;
                        aw[i][5] -= trans * vflowy3 * sflow.RHO3 / sflow.RHO2;

                        aw[i][6] += trans;
                        aw[i][7] += trans * vflowx3;
                        aw[i][8] += trans * vflowy3;
                    }
                }
            }


// -- STOP --- Phase transformations ----------------------------------------------------------------------------


// -- START -- Entrainment --------------------------------------------------------------------------------------


            // *** 1 = Entrainment coefficient multiplied with flow momentum
            // *** 2 = Pudasaini and Fischer (2020) and Pudasaini and Krautblatter (2021) erosion-deposition models
            // *** 3 = Combined approach: 1 for entrainment, 2 for deposition


            if ( sico.ENTRAINMENT != 0 ) {

                if ( cslide == 0 ) iloop = ib[0]; else iloop = sico.IMAX;
                for ( ix=0; ix<iloop; ix++ ) {

                    if ( cslide == 0 ) i = ibasket[0][ix]; else i = ix;
                    qentr = 0;

                    if ( sico.CENTR == 1 ) {
                        if ( pcentr[i] != sico.UNDEF ) sflow.CENTR = pcentr[i]; // entrainment coefficient
                        else sflow.CENTR = sflow.CENTR0;
                    } else sflow.CENTR = sflow.CENTR0;

                    if ( sico.MODEL <= 3 ) { hflow = aw[i][0]; alpha1 = 1; } // flow heights and phase fractions
                    else { hflow = aw[i][0] + aw[i][3] + aw[i][6]; alpha1 = aw[i][0] / hflow; alpha2 = aw[i][3] / hflow; alpha3 = aw[i][6] / hflow; }
                    
                    if ( sico.MODEL <= 3 ) { hentrmaxx = phentrmax[i]; alphab1 = 1; } // entrainable heights and phase fractions
                    else {
                     
                        hentrmaxx = phentrmax[i] + phentrmax2[i] + phentrmax3[i]; 
                        if ( hentrmaxx > 0 ) { alphab1 = phentrmax[i] / hentrmaxx; alphab2 = phentrmax2[i] / hentrmaxx; alphab3 = phentrmax3[i] / hentrmaxx; }
                        else { alphab1 = 1; alphab2 = 0; alphab3 = 0; }
                    }

                    if ( aw[i][0] > 0 ) { vflowx1 = aw[i][1]/aw[i][0]; vflowy1 = aw[i][2]/aw[i][0]; } else { vflowx1 = 0; vflowy1 = 0; } // flow velocities of individual phases
                    if ( sico.MODEL == 7 && aw[i][3] > 0 ) { vflowx2 = aw[i][4]/aw[i][3]; vflowy2 = aw[i][5]/aw[i][3]; } else { vflowx2 = 0; vflowy2 = 0; }
                    if ( sico.MODEL == 7 && aw[i][6] > 0 ) { vflowx3 = aw[i][7]/aw[i][6]; vflowy3 = aw[i][8]/aw[i][6]; } else { vflowx3 = 0; vflowy3 = 0; }

                    for ( j=0; j<9; j++ ) {

                        if ( sico.MODEL <= 3 && hflow > sico.HFLOWMIN ) hflowj[j] = aw[in[i][j]][0];
                        else if ( hflow > sico.HFLOWMIN ) hflowj[j] = aw[in[i][j]][0] + aw[in[i][j]][3] + aw[in[i][j]][6];
                        elevtot[j] = pelev[in[i][j]]; // elevation of neighbour cell
                        
                        if ( sico.MODEL <= 3 && hflowj[j] > sico.HFLOWMIN && hflow > sico.HFLOWMIN ) {
                            
                            vflowxj[j] = aw[in[i][j]][1] / hflowj[j]; // flow velocities of neighbour cell
                            vflowyj[j] = aw[in[i][j]][2] / hflowj[j];
                            
                        } else if ( sico.MODEL == 7 && hflowj[j] > sico.HFLOWMIN && hflow > sico.HFLOWMIN ) {
                            
                            vflowxj[j] = ( aw[in[i][j]][1] + aw[in[i][j]][4] + aw[in[i][j]][5] ) / hflowj[j];
                            vflowyj[j] = ( aw[in[i][j]][2] + aw[in[i][j]][5] + aw[in[i][j]][8] ) / hflowj[j]; 

                        } else { vflowxj[j] = 0; vflowyj[j] = 0; }
                    }

                    wbetax = fbeta( elevtot[5], elevtot[2], 2.0, sico );
                    wbetay = fbeta( elevtot[4], elevtot[1], 2.0, sico );
                    wbetaxy = fbetaxy( wbetax, wbetay ); // slope, including flow height

                    alphav = falphav( vflowxj[0], vflowyj[0], sico ); // direction of movement
                    alpha = falpha( wbetax, wbetay, sico ); // aspect

                    if (( sico.ENTRAINMENT == 1 || sico.ENTRAINMENT == 3 ) && hflow > sico.HFLOWMIN && hentrmaxx > 0 ) { // first model

                        mom = aw[i][0] * sflow.RHO1 * pow( pow( vflowx1, 2 ) + pow( vflowy1, 2 ), 0.5 ); // flow momenta
                        if ( sico.MODEL == 7 ) mom += aw[i][3] * sflow.RHO2 * pow( pow( vflowx2, 2 ) + pow( vflowy2, 2 ), 0.5 );
                        if ( sico.MODEL == 7 ) mom += aw[i][6] * sflow.RHO3 * pow( pow( vflowx3, 2 ) + pow( vflowy3, 2 ), 0.5 );

                        if ( sflow.CENTR > 0 ) betav = ffmax(0, tan( wbetaxy * cos( alpha - alphav ))); // movement-following slope
                        else betav = 1;
                
                        qentr = tlength * pow( 10, -fabs( sflow.CENTR )) * mom * betav; // entrainment rate
                        momfact = 0; // mobility generator
                  
                    } else if (( sico.ENTRAINMENT == 2 || sico.ENTRAINMENT == 3 ) && hflow > sico.HFLOWMIN ) { // second model

                        if ( sico.PHASES[0] != 3 ) { // if mixture, solid, or fine solid material is involved

                            if ( sico.DELTA == 1 ) {
                                if ( pdelta[i] != sico.UNDEF ) sflow.DELTA[0] = pdelta[i]; // basal friction angle of PHASE 1
                                else sflow.DELTA[0] = sflow.DELTA0[0];
                            } else sflow.DELTA[0] = sflow.DELTA0[0];

                            if ( sico.DELTA2 == 1 ) {
                                if ( pdelta2[i] != sico.UNDEF ) sflow.DELTA[1] = pdelta2[i]; // basal friction angle of PHASE 2
                                else sflow.DELTA[1] = sflow.DELTA0[1];
                            } else sflow.DELTA[1] = sflow.DELTA0[1];

                            if ( sico.CVSHEAR == 1 ) {
                                if ( pcvshear[i] != sico.UNDEF ) sflow.CVSHEAR = pcvshear[i]; // shear velocity ratio
                                else sflow.CVSHEAR = sflow.CVSHEAR0;
                            } else sflow.CVSHEAR = sflow.CVSHEAR0;
                            sflow.CVSHEAR = pow( 1 / sflow.CVSHEAR, 2 );

                            if ( sico.DELTAB == 1 ) {
                                if ( pdeltab[i] != sico.UNDEF ) sflow.DELTAB = pdeltab[i]; // basal friction difference
                                else sflow.DELTAB = sflow.DELTAB0;
                            } else sflow.DELTAB = sflow.DELTAB0;

                            if ( sflow.CVSHEAR > 0 ) betav = wbetaxy;
                            else betav = ffmax( 0, wbetaxy * cos( alpha - alphav ));
                            gz = sico.GRAVITY * cos( betav ) * hflow; // normal force

                            if ( sico.MODEL <= 3 ) {

                                alphasfs = 1; // mixture characteristics (one-phase models)
                                rhom = sflow.RHO1;
                                mym = tan( sflow.DELTA[0] );
                                gammam = 0;
                            
                                alphabsfs = 1; // characteristics of basal surface (one-phase models)
                                rhob = sflow.RHOB1;
                                myb = tan( sflow.DELTAB ) + mym;
                                gammab = 0;
        
                            } else {
                                                       
                                alphasfs = ( aw[i][0] + aw[i][3] ) / hflow; // mixture characteristics (multi-phase model)
                                rhom = sflow.RHO1 * alpha1 + sflow.RHO2 * alpha2 + sflow.RHO3 * alpha3;
                                mym = tan( sflow.DELTA[0] * alpha1 + sflow.DELTA[1] * alpha2 );
                                if ( alphasfs > 0 ) gammam = sflow.RHO3 / (( aw[i][0] * sflow.RHO1 + aw[i][3] * sflow.RHO2 ) / ( aw[i][0] + aw[i][3] )); else gammam = 0;
                            
                                if ( phentrmax[i] + phentrmax2[i] > 0 ) {
                            
                                    alphabsfs = ( phentrmax[i] + phentrmax2[i] ) / ( hentrmaxx ); // characteristics of basal surface (multi-phase model)
                                    rhob = sflow.RHOB1 * alphab1 + sflow.RHOB2 * alphab2 + sflow.RHOB3 * alphab3;
                                    if ( alphabsfs > 0 ) gammab = sflow.RHOB3 / (( phentrmax[i] * sflow.RHOB1 + phentrmax2[i] * sflow.RHOB2 ) / ( phentrmax[i] + phentrmax2[i] )); else gammab = 0;
                                
                                } else { alphabsfs = 1; rhob = sflow.RHO1; gammab = 0; }
                                
                                myb = tan( sflow.DELTAB ) + mym;
                            }
                            
                            lambdam = 1.0; // drift factors
                            lambdab = lambdam / ( 1 + ( rhob * alphabsfs / rhom / alphasfs ));
                            
                            qentrup = gz * (( 1 - gammam ) * rhom * mym * alphasfs - ( 1 - gammab ) * rhob * myb * alphabsfs ); // components of entrainment rate equation
                            qentrdown = fabs( sflow.CVSHEAR ) * ( rhom * lambdam * alphasfs - rhob * lambdab * alphabsfs );
                            if ( sico.MODEL == 7 ) qentrtest = tlength * pow( 10, -fabs( sflow.CENTR )) * aw[i][6] * sflow.RHO3 * pow( pow( vflowx3, 2 ) + pow( vflowy3, 2 ), 0.5 );
                            else qentrtest = 0;

                            if ( qentrdown != 0 ) { // entrainment rate and mobility generator
                            
                                if ( qentrup > 0 ) {
                                
                                    qentrtest += tlength * pow ( fabs( qentrup ), 0.5 ) / pow( fabs( qentrdown ), 0.5 );
                                    momfacttest = 2 * lambdab - 1;
                                    
                                } else {
                                
                                    qentrtest -= tlength * pow ( fabs( qentrup ), 0.5 ) / pow( fabs( qentrdown ), 0.5 );
                                    momfacttest = 1;
                                }
                            }

                        } else { // if PHASE 1 is fluid
                            
                            qentrtest = tlength * pow( 10, -fabs( sflow.CENTR )) * aw[i][0] * sflow.RHO1 * pow( pow( vflowx, 2 ) + pow( vflowy, 2 ), 0.5 );
                            momfact = 0;
                        }
                        
                        if ( sico.ENTRAINMENT == 2 ) { qentr = qentrtest; momfact = momfacttest; } // managing combined approach (third model)
                        else if ( qentrtest < 0 ) { qentr += qentrtest; if ( qentr > 0 ) momfact = 0; else momfact = 1; }
                        
                    } else if ( sico.ENTRAINMENT == 4 && hflow > sico.HFLOWMIN ) { // fourth model
                    
                        dux = ( vflowxj[2] * hflowj[2] - vflowxj[5] * hflowj[5] ) / ( 2 * sico.CSZ ); duy = ( vflowyj[1] * hflowj[1] - vflowyj[4] * hflowj[4] ) / ( 2 * sico.CSZ );
                            // velocity gradients in x and y direction
                        alpha = falphav( dux, duy, sico ); // direction of velocity gradient
                        duxy = pow( pow( dux, 2 ) + pow( duy, 2 ), 0.5 ); // absolute value of velocity gradient
                        dumain = duxy * cos( alpha - alphav ); // movement-following velocity gradient
                        //dumain = duxy * cos( alpha - alphav ) / pow( vflowxj[0] * vflowxj[0] + vflowyj[0] * vflowyj[0], 0.5 ); // movement-following velocity gradient
                        
                        if ( sico.DELTAB == 1 ) {
                            if ( pdeltab[i] != sico.UNDEF ) sflow.DELTAB = pdeltab[i]; // basal friction difference
                            else sflow.DELTAB = sflow.DELTAB0;
                        } else sflow.DELTAB = sflow.DELTAB0;

                        for ( l=0; l<9; l++ ) welev[l] = pelev[in[i][l]]; // elevation of adjacent cells
                        kappau = fcurv( vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, welev, grav, sico ); // curvature
                        betav = cos( wbetaxy * cos( alpha - alphav )); // cosine of movement-following slope
                                            
                        qentr = tlength * betav * fsign( dumain ) * fabs( pow( dumain, sflow.DELTAB ) * pow( 10, -fabs( sflow.CENTR ))); // entrainment rate
                        if ( qentr < 0 && kappau < 0 ) qentr = 0; // avoiding deposition in convex topography
                        
                        momfact = 0; // mobility generator
                    }

                    if ( qentr != 0 ) {

                        if ( qentr > 0 ) qentr1 = ffmin( phentrmax[i], alphab1 * qentr ); // constraining entrainment rate
                        else qentr1 = ffmax( -aw[i][0], alpha1 * qentr );
                        phentrmax[i] -= ffmax( 0, qentr1 ); // updating entrainable height
                        
                        momaddsx = momfact * vflowx1 * qentr1; // contribution of entrainment to momentum
                        momaddsy = momfact * vflowy1 * qentr1;

                        aw[i][0] = aw[i][0] + qentr1; // updating flow depths and momenta (mixture or PHASE 1)
                        aw[i][1] = aw[i][1] + momaddsx;
                        aw[i][2] = aw[i][2] + momaddsy;
                        
                        if ( sico.MODEL == 7 ) {
                        
                            if ( qentr > 0 ) qentr2 = ffmin( phentrmax2[i], alphab2 * qentr ); // constraining entrainment rates
                            else qentr2 = ffmax( -aw[i][3], alpha2 * qentr );
                            if ( qentr > 0 ) qentr3 = ffmin( phentrmax3[i], alphab3 * qentr );
                            else {
                            
                                qentr3 = ffmax( -aw[i][6], alpha3 * qentr );
                                qentr3 = ffmax( qentr3, sflow.THETAS * ( qentr1 + qentr2 ));
                            }
                            phentrmax2[i] -= ffmax( 0, qentr2 ); // updating entrainable heights
                            phentrmax3[i] -= ffmax( 0, qentr3 );
                            
                            momaddfsx = momfact * vflowx2 * qentr2;  // contribution of entrainment to momenta
                            momaddfx = momfact * vflowx3 * qentr3;
                            momaddfsy = momfact * vflowy2 * qentr2;
                            momaddfy = momfact * vflowy3 * qentr3;

                            aw[i][9] = aw[i][9] - qentr1; // updating change of basal topography (PHASE 1)

                            aw[i][3] = aw[i][3] + qentr2; // updating flow depths, momenta, and change of basal topography (PHASE 2)
                            aw[i][4] = aw[i][4] + momaddfsx;
                            aw[i][5] = aw[i][5] + momaddfsy;
                            aw[i][10] = aw[i][10] - qentr2; 

                            aw[i][6] = aw[i][6] + qentr3; // updating flow depths, momenta, and change of basal topography (PHASE 3)
                            aw[i][7] = aw[i][7] + momaddfx;
                            aw[i][8] = aw[i][8] + momaddfy;
                            aw[i][11] = aw[i][11] - qentr3;

                        } else aw[i][3] = aw[i][3] - qentr1; // updating change of basal topography
                    }
                }
            }


// -- STOP --- Entrainment --------------------------------------------------------------------------------------


// -- START -- Stopping of flow ---------------------------------------------------------------------------------


            // *** 1 = Stopping based on fraction of maximum kinetic energy, material is deposited and simulation is terminated
            // ***     when stopping occurs
            // *** 2 = Stopping based on fraction of maximum momentum, material is deposited and simulation is terminated
            // ***     when stopping occurs
            // *** 3 = Stopping based on flow pressure threshold, material is deposited 
            // ***     and simulation is terminated when stopping occurs
            // *** Negative numbers allow for the release of the fluid material along with an equal amount of solid material (multi-phase model only) 


            if ( tsum > 0 && sico.STOPPING != 0 ) { // stopping through fraction of maximum kinetic energy or momentum, or minimum pressure criterion

                hekin_sum = 0; // resetting sum of momentum or kinetic energy or momentum
                ctrl_pressthr = 0; // resetting minimum pressure criterion

                if ( cslide == 0 ) iloop = ib[0]; else iloop = sico.IMAX;
                for ( ix=0; ix<iloop; ix++ ) {

                    if ( cslide == 0 ) i = ibasket[0][ix]; else i = ix;

                    if ( cdomain2[i] == 0 && tsum >= fabs( pstoptime[i] )) {

                        if ( aw[i][0] > sico.HFLOWMIN ) { vflowx = aw[i][1]/aw[i][0]; vflowy = aw[i][2]/aw[i][0]; } else { vflowx = 0; vflowy = 0; } // flow velocities
                        if ( sico.MODEL == 7 && aw[i][3] > sico.HFLOWMIN ) { vflowx2 = aw[i][4]/aw[i][3]; vflowy2 = aw[i][5]/aw[i][3]; } else { vflowx2 = 0; vflowy2 = 0; }
                        if ( sico.MODEL == 7 && aw[i][6] > sico.HFLOWMIN ) { vflowx3 = aw[i][7]/aw[i][6]; vflowy3 = aw[i][8]/aw[i][6]; } else { vflowx3 = 0; vflowy3 = 0; }

                        if ( fabs( sico.STOPPING ) == 1 ) { // fraction of kinetic energy

                            hekin = 0.5 * aw[i][0] * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 )); // flow kinetic energies
                            if ( sico.MODEL == 7 ) hekin += 0.5 * aw[i][3] * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * aw[i][6] * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));
                            hekin_sum += hekin; // updating sum of kinetic energy

                        } else if ( fabs( sico.STOPPING ) == 2 ) { // fraction of momentum

                            hekin = aw[i][0] * sflow.RHO1 * pow( pow( vflowx, 2 ) + pow( vflowy, 2 ), 0.5 ); // flow momenta
                            if ( sico.MODEL == 7 ) hekin += aw[i][3] * sflow.RHO2 * pow( pow( vflowx2, 2 ) + pow( vflowy2, 2 ), 0.5 );
                            if ( sico.MODEL == 7 ) hekin += aw[i][6] * sflow.RHO3 * pow( pow( vflowx3, 2 ) + pow( vflowy3, 2 ), 0.5 );
                            hekin_sum += hekin; // updating sum of momentum
                            
                        } else if ( fabs( sico.STOPPING ) == 3 ) { // minimum pressure criterion                        
                        
                            if ( sico.MODEL <= 3 ) hflow = aw[i][0]; // flow heights
                            else if ( sico.MODEL == 7 ) hflow = aw[i][0] + aw[i][3] + aw[i][6];

                            hekin = 0.5 * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 )); // flow pressures
                            if ( sico.MODEL == 7 ) hekin += 0.5 * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));
                            if ( hflow > sico.IMPTHR[0] ) { if ( hekin > sflow.CSTOP ) ctrl_pressthr = 1; } // updating minimum pressure criterion
                        }
                    }
                }

                if ( hekin_sum > hekin_max ) hekin_max = hekin_sum; // updating maximum momentum or kinetic energy

                if ( cslide == 0 ) iloop = ib[0]; else iloop = sico.IMAX;
                for ( ix=0; ix<iloop; ix++ ) {

                    if ( cslide == 0 ) i = ibasket[0][ix]; else i = ix;

                    if ( tsum >= fabs( pstoptime[i] )) {

                        if ( pstoptime[i] < 0 || ( fabs( sico.STOPPING ) == 3 && ctrl_pressthr == 0 ) || ( fabs( sico.STOPPING ) != 3 && hekin_sum < hekin_max * sflow.CSTOP )) { 
                            // applying stopping criterion

                            cstopped[i] = 1; // updating control for stopping

                            if ( sico.STOPPING < 0 && sico.MODEL == 7 ) { // stopping with material release
                            
                                for ( j=0; j<imax; j++ ) {

                                    aw[in[i][j]][0] = aw[in[i][j]][6]; aw[in[i][j]][1] = aw[in[i][j]][7]; aw[in[i][j]][2] = aw[in[i][j]][8]; // constraining PHASE 1
                                }
                            }
                        }
                    }
                }
            }


// -- STOP --- Stopping of flow ---------------------------------------------------------------------------------


// -- START -- Updating vectors of composite, maximum and cumulative values, and basal topography ---------------


            if ( cslide == 0 ) iloop = ib[0]; else iloop = sico.IMAX;
            for ( ix=0; ix<iloop; ix++ ) {

                if ( cslide == 0 ) i = ibasket[0][ix]; else i = ix;

                if ( sico.CORRHEIGHT == 0 ) corrfact = 1; // correction factor for depth to height conversion
                else corrfact = 1 / cos( betaxy[i] );

                if ( sico.MODEL <= 3 ) { // one-phase models

                    pelev[i] = pelev0[i] + aw[i][3] * corrfact; // correcting elevation for entrainment and deposition

                    if ( aw[i][0] > sico.HFLOWMIN ) aw[i][4] = pow( pow( aw[i][1]/aw[i][0], 2 ) + pow( aw[i][2]/aw[i][0], 2 ) , 0.5 ); else aw[i][4] = 0; // flow velocity
                    aw[i][5] = 0.5 * sflow.RHO1 * aw[i][0] * pow( aw[i][4], 2 ); // flow kinetic energy
                    aw[i][6] = 0.5 * sflow.RHO1 * pow( aw[i][4], 2 ); // flow pressure

                    if ( aw[i][0] > aw[i][7] ) aw[i][7] = aw[i][0]; // maximum flow depth
                    if ( aw[i][4] > aw[i][8] && aw[i][0] >= sico.IMPTHR[0] ) aw[i][8] = aw[i][4]; // maximum flow velocity
                    if ( aw[i][5] > aw[i][9] ) aw[i][9] = aw[i][5]; //maximum flow kinetic energy
                    if ( aw[i][6] > aw[i][10] ) aw[i][10] = aw[i][6]; // maximum flow pressure

                    if ( aw[i][8] > vmax ) vmax = aw[i][8]; // overall maximum flow velocity
                
                } else if ( sico.MODEL == 7 ) { // multi-phase model

                    pelev[i] = pelev0[i] + ( aw[i][9] + aw[i][10] + aw[i][11] ) * corrfact; // correcting elevation for entrainment and deposition

                    if ( aw[i][0] > sico.HFLOWMIN ) aw[i][12] = pow( pow( aw[i][1]/aw[i][0], 2 ) + pow( aw[i][2]/aw[i][0], 2 ) , 0.5 ); else aw[i][12] = 0;
                    if ( aw[i][3] > sico.HFLOWMIN ) aw[i][13] = pow( pow( aw[i][4]/aw[i][3], 2 ) + pow( aw[i][5]/aw[i][3], 2 ) , 0.5 ); else aw[i][13] = 0;
                    if ( aw[i][6] > sico.HFLOWMIN ) aw[i][14] = pow( pow( aw[i][7]/aw[i][6], 2 ) + pow( aw[i][8]/aw[i][6], 2 ) , 0.5 ); else aw[i][14] = 0; // flow velocities
                    
                    aw[i][15] = aw[i][0] + aw[i][3] + aw[i][6]; // total flow depth
                    
                    aw[i][16] = 0.5 * sflow.RHO1 * aw[i][0] * pow( aw[i][12], 2 );
                    aw[i][17] = 0.5 * sflow.RHO2 * aw[i][3] * pow( aw[i][13], 2 );
                    aw[i][18] = 0.5 * sflow.RHO3 * aw[i][6] * pow( aw[i][14], 2 );
                    aw[i][19] = aw[i][16] + aw[i][17] + aw[i][18]; // flow kinetic energies
                    
                    if ( aw[i][0] > sico.HFLOWMIN ) aw[i][20] = aw[i][16] / aw[i][0]; else aw[i][20] = 0;
                    if ( aw[i][3] > sico.HFLOWMIN ) aw[i][21] = aw[i][17] / aw[i][3]; else aw[i][21] = 0;
                    if ( aw[i][6] > sico.HFLOWMIN ) aw[i][22] = aw[i][18] / aw[i][6]; else aw[i][22] = 0;
                    if ( aw[i][15] > sico.HFLOWMIN ) aw[i][23] = aw[i][19] / aw[i][15]; else aw[i][23] = 0; // flow pressures
                    
                    aw[i][24] = aw[i][9] + aw[i][10] + aw[i][11]; // total entrained or deposited depth

                    if ( aw[i][0] > aw[i][25] ) aw[i][25] = aw[i][0];
                    if ( aw[i][12] > aw[i][26] && aw[i][0] >= sico.IMPTHR[0] ) aw[i][26] = aw[i][12];
                    if ( aw[i][3] > aw[i][27] ) aw[i][27] = aw[i][3];
                    if ( aw[i][13] > aw[i][28] && aw[i][3] >= sico.IMPTHR[0] ) aw[i][28] = aw[i][13];
                    if ( aw[i][6] > aw[i][29] ) aw[i][29] = aw[i][6];
                    if ( aw[i][14] > aw[i][30] && aw[i][6] >= sico.IMPTHR[0] ) aw[i][30] = aw[i][14]; // maximum flow depths and velocities

                    if ( pow( pow( aw[i][1] + aw[i][4] + aw[i][7], 2 ) + pow( aw[i][2] + aw[i][5] + aw[i][8], 2 ), 0.5 ) > aw[i][47] ) {
                    
                        aw[i][47] = pow( pow( aw[i][1] + aw[i][4] + aw[i][7], 2 ) + pow( aw[i][2] + aw[i][5] + aw[i][8], 2 ), 0.5 ); // maximum total momentum
                        aw[i][45] = pow( pow(( aw[i][1] + aw[i][4] + aw[i][7] ) / ( aw[i][0] + aw[i][3] + aw[i][6] ), 2 )
                            + pow(( aw[i][2] + aw[i][5] + aw[i][8] ) / ( aw[i][0] + aw[i][3] + aw[i][6] ), 2 ), 0.5 ); // flow velocity at maximum total momentum
                        aw[i][43] = aw[i][0] / aw[i][15];
                        aw[i][44] = aw[i][6] / aw[i][15]; // phase fractions at maximum total momentum
                    }

                    if ( aw[i][15] > aw[i][31] ) aw[i][31] = aw[i][15]; // maximum total flow depth

                    if ( aw[i][16] > aw[i][32] ) aw[i][32] = aw[i][16];
                    if ( aw[i][20] > aw[i][33] ) aw[i][33] = aw[i][20];
                    if ( aw[i][17] > aw[i][34] ) aw[i][34] = aw[i][17];
                    if ( aw[i][21] > aw[i][35] ) aw[i][35] = aw[i][21];
                    if ( aw[i][18] > aw[i][36] ) aw[i][36] = aw[i][18];
                    if ( aw[i][22] > aw[i][37] ) aw[i][37] = aw[i][22];
                    if ( aw[i][19] > aw[i][38] ) aw[i][38] = aw[i][19];
                    if ( aw[i][23] > aw[i][39] ) aw[i][39] = aw[i][23]; // maximum kinetic energies and flow pressures

                    if ( aw[i][26] > vmax ) vmax = aw[i][26];
                    if ( aw[i][28] > vmax ) vmax = aw[i][28];
                    if ( aw[i][30] > vmax ) vmax = aw[i][30]; // overall maximum velocity
                }
            }


// -- STOP --- Updating vectors of composite, maximum and cumulative values, and basal topography ---------------


// *** End of loop over two steps, each moving the system half of a cell (NOC scheme) --------------------------


        }

        if ( cslide == 0 ) iloop = ib[0]; else iloop = sico.IMAX;
        for ( ix=0; ix<iloop; ix++ ) {

            if ( cslide == 0 ) i = ibasket[0][ix]; else i = ix;

            betax[i] = fbeta( pelev[in[i][5]], pelev[in[i][2]], 2.0, sico ); // slopes
            betay[i] = fbeta( pelev[in[i][4]], pelev[in[i][1]], 2.0, sico );
            betaxy[i] = fbetaxy( betax[i], betay[i] );

            if ( sico.CORRHEIGHT != 0 ) {

                dx[i] = sico.CSZ / cos( betax[i] ); // topography-following cell spacing
                dy[i] = sico.CSZ / cos( betay[i] );
            }
        }


// -- START -- Updating maximum values of flow parameters and volumes, time of reach ----------------------------


        hflow_max = 0; vflow_max = 0; vol_flow = 0; vol_edge = 0; vol_entr = 0; vol_noflux = 0; ekin_flow = 0;
        hflow_max2 = 0; vflow_max2 = 0; vol_flow2 = 0; vol_edge2 = 0; vol_entr2 = 0; vol_noflux2 = 0;
        hflow_max3 = 0; vflow_max3 = 0; vol_flow3 = 0; vol_edge3 = 0; vol_entr3 = 0; vol_noflux3 = 0;

        for ( z=0; z<nzones; z++ ) { vol_zone1[z] = 0; vol_zone2[z] = 0; vol_zone3[z] = 0; vol_czone1[z] = 0; vol_czone2[z] = 0; vol_czone3[z] = 0; } // zone-specific volumes

        for ( i=0; i<sico.IMAX; i++ ) {

            if ( cdomain2[i] == 0 ) {

                if ( sico.CORRHEIGHT == 0 ) carea = pow ( sico.CSZ, 2 );
                else carea = pow ( sico.CSZ, 2 ) * pow( 1 - pow( sin( betax[i] ) , 2 ) * pow ( sin( betay[i] ) , 2 ) , 0.5 )
                    / ( cos( betax[i] ) * cos( betay[i] ) ); // topography-following area of cell

                if ( aw[i][0] > hflow_max ) hflow_max = aw[i][0]; // maximum mixture or PHASE 1 flow depth
                if ( cstopped[i] == 0 && cflux[i] == 1 && aw[i][0] > sico.HFLOWMIN ) { vol_flow += aw[i][0] * carea; vol_zone1[pzones[i]] += aw[i][0] * carea; } // mixture or PHASE 1 flow volume
                if (( cflux[i] == 0 || cstopped[i] == 1 ) && aw[i][0] > sico.HFLOWMIN ) vol_noflux += aw[i][0] * carea; // mixture or PHASE 1 volume without fluxes

                if ( sico.MODEL <= 3 ) { // one-phase models

                    if ( hflow_max > hflow_maxmax ) hflow_maxmax = hflow_max; // absolute maximum flow depth
                    if ( aw[i][0] >= sico.IMPTHR[0] && aw[i][4] > vflow_max ) vflow_max = aw[i][4]; // maximum flow velocity
                    if ( cdomain[i] != 0 && ( sico.ENTRAINMENT != 0 || sico.STOPPING != 0 )) { vol_entr += aw[i][3] * carea; vol_czone1[pzones[i]] += aw[i][3] * carea; } // entrained or deposited volume

                    ekin_flow += aw[i][5] * carea; // kinetic energy of flow
                    
                } else if ( sico.MODEL == 7 ) { // multi-phase model

                    if ( aw[i][3] > hflow_max2 ) hflow_max2 = aw[i][3];
                    if ( aw[i][6] > hflow_max3 ) hflow_max3 = aw[i][6];
                    if ( hflow_max + hflow_max2 + hflow_max3 > hflow_maxmax ) hflow_maxmax = hflow_max + hflow_max2 + hflow_max3; // absolute maximum flow depths

                    if ( cdomain[i] != 0 && aw[i][0] >= sico.IMPTHR[0] && aw[i][12] > vflow_max ) vflow_max = aw[i][12];
                    if ( cdomain[i] != 0 && aw[i][3] >= sico.IMPTHR[0] && aw[i][13] > vflow_max2 ) vflow_max2 = aw[i][13];
                    if ( cdomain[i] != 0 && aw[i][6] >= sico.IMPTHR[0] && aw[i][14] > vflow_max3 ) vflow_max3 = aw[i][14]; // maximum velocities

                    if ( cstopped[i] == 0 && cflux[i] == 1 && aw[i][3] > sico.HFLOWMIN ) { vol_flow2 += aw[i][3] * carea; vol_zone2[pzones[i]] += aw[i][3] * carea; }
                    if ( cstopped[i] == 0 && cflux[i] == 1 && aw[i][6] > sico.HFLOWMIN ) { vol_flow3 += aw[i][6] * carea; vol_zone3[pzones[i]] += aw[i][6] * carea; } // flow volumes

                    if (( cflux[i] == 0 || cstopped[i] == 1 ) && aw[i][3] > sico.HFLOWMIN ) vol_noflux2 += aw[i][3] * carea;
                    if (( cflux[i] == 0 || cstopped[i] == 1 ) && aw[i][6] > sico.HFLOWMIN ) vol_noflux3 += aw[i][6] * carea; // volumes without fluxes

                    if ( cdomain[i] != 0 && ( sico.ENTRAINMENT != 0 || sico.STOPPING != 0 )) { vol_entr += aw[i][9] * carea; vol_czone1[pzones[i]] += aw[i][9] * carea; }
                    if ( cdomain[i] != 0 && ( sico.ENTRAINMENT != 0 || sico.STOPPING != 0 )) { vol_entr2 += aw[i][10] * carea; vol_czone2[pzones[i]] += aw[i][10] * carea; }
                    if ( cdomain[i] != 0 && ( sico.ENTRAINMENT != 0 || sico.STOPPING != 0 )) { vol_entr3 += aw[i][11] * carea; vol_czone3[pzones[i]] += aw[i][11] * carea; }
                        // entrained or deposited volumes

                    ekin_flow += aw[i][19] * carea; // kinetic energy of flow
                }
            }

            if ( cdomain2[i] == 1 && aw[i][0] > sico.HFLOWMIN ) vol_edge += aw[i][0] * carea;
            if ( sico.MODEL == 7 && cdomain2[i] == 1 && aw[i][3] > sico.HFLOWMIN ) vol_edge2 += aw[i][3] * carea;
            if ( sico.MODEL == 7 && cdomain2[i] == 1 && aw[i][6] > sico.HFLOWMIN ) vol_edge3 += aw[i][6] * carea; // flow volumes leaving area of interest
        }

        for ( i=0; i<sico.IMAX; i++ ) { // updating time of reach

            if ( cdomain2[i] == 0 ) {

                if ( sico.MODEL <= 3 ) hflow = aw[i][0];
                else if ( sico.MODEL == 7 ) hflow = aw[i][0] + aw[i][3] + aw[i][6];
                if ( aw[i][nvect_all-2] == sico.UNDEF && hflow >= sico.IMPTHR[0] ) {

                    aw[i][nvect_all-8] = pow( pow(( aw[i][1] + aw[i][4] + aw[i][7] ) / ( aw[i][0] + aw[i][3] + aw[i][6] ), 2 )
                        + pow(( aw[i][2] + aw[i][5] + aw[i][8] ) / ( aw[i][0] + aw[i][3] + aw[i][6] ), 2 ), 0.5 ); 
                    aw[i][nvect_all-7] = aw[i][0] / hflow;
                    aw[i][nvect_all-6] = aw[i][3] / hflow; 
                    aw[i][nvect_all-2] = tsum;
                }
            }
        }


// -- STOP --- Updating maximum values of flow parameters and volumes, time of reach ----------------------------


// -- START -- Managing stopped flows and numerical failures ----------------------------------------------------


        if (( sico.MODEL <= 3 && vol_flow == 0 ) || ( sico.MODEL == 7 && vol_flow == 0 && vol_flow2 == 0 && vol_flow3 == 0 )) {
            // stopping simulation if the entire flow has stopped moving or left the area of interest (successful simulation)

            printf("\n                       \n                       ------------------ FLOW STOPPED ------------------\n\n");
            fflush(stdout); // forcing immediate display
            fprintf(f_summary, "\n                       \n                       ------------------ FLOW STOPPED ------------------\n\n");
            ccontinue = 0; // setting control for continuation to negative
        }

        if ( ccontinue == 0 && csuccess == 0 ) {
            // stopping simulation if time step length indicates numerical failure (unsuccessful simulation)

            printf("\n                       \n                       --------------- NUMERICAL FAILURE ----------------\n\n");
            fflush(stdout); // forcing immediate display
            fprintf(f_summary, "\n                       \n                       --------------- NUMERICAL FAILURE ----------------\n\n");
        }

        if ( sico.STOPPING != 0 && ccontinue == 0 && csuccess == 1 ) { // correcting deposited depth and basal topography for stopping

            if ( cslide == 0 ) iloop = ib[0]; else iloop = sico.IMAX;
            for ( ix=0; ix<iloop; ix++ ) {

                if ( cslide == 0 ) i = ibasket[0][ix]; else i = ix;

                if ( sico.CORRHEIGHT == 0 ) corrfact = 1; // correction factor for depth to height conversion
                else corrfact = 1 / cos( betaxy[i] );

                if ( sico.MODEL <= 3 ) { // one-phase models

                    aw[i][3] += cstopped[i] * aw[i][0]; // change of basal surface
                    pelev[i] += cstopped[i] * aw[i][0] * corrfact; // elevation
                    aw[i][0] -= cstopped[i] * aw[i][0]; // flow depth
                    
                } else if ( sico.MODEL == 7 ) { // multi-phase model

                    aw[i][9] += cstopped[i] * aw[i][0];
                    aw[i][10] += cstopped[i] * aw[i][3];
                    aw[i][11] += cstopped[i] * aw[i][6];
                    aw[i][17] = aw[i][9] + aw[i][10] + aw[i][11]; // change of basal surface
                    
                    pelev[i] += ( cstopped[i] * aw[i][0] + cstopped[i] * aw[i][3] + cstopped[i] * aw[i][6] ) * corrfact; // elevation
                    aw[i][0] -= cstopped[i] * aw[i][0];
                    aw[i][3] -= cstopped[i] * aw[i][3];
                    aw[i][6] -= cstopped[i] * aw[i][6]; // flow depths
                }
            }
        }


// -- STOP --- Managing stopped flows and numerical failures ----------------------------------------------------


        if ( (int)(round( 1000 * tint )) >= (int)( round( 1000 * tout )) || ccontinue == 0 ) { // if defined interval or last time step is reached


// -- START -- Writing hydrograph infos to files ----------------------------------------------------------------


            if ( sico.MULT == 0 && hydrograph == 1 ) {

                ctrl_hydout = 1; // control for hydrograph output

                for ( hydj = hydnin; hydj < hydnin + hydnout; hydj++ ) {  // loop over all output hydrographs

                    hydh = aw[hydi[hydj]][0]  / cos( betaxy[hydi[hydj]] ); // mixture or PHASE 1 flow height
                    if ( hydh > sico.HFLOWMIN ) hydv = pow( pow( aw[hydi[hydj]][1], 2 ) + pow( aw[hydi[hydj]][2], 2 ), 0.5 ) / aw[hydi[hydj]][0];
                    else hydv= 0; // mixture or PHASE 1 flow velocity

                    if ( sico.MODEL <= 3 ) {

                        hyde = aw[hydi[hydj]][3]; // entrained height
                        hydh2 = 0; hydv2 = 0; hyde2 = 0; hydh3 = 0; hydv3 = 0; hyde3 = 0;
                        
                    } else if ( sico.MODEL == 7 ) {

                        hydh2 = aw[hydi[hydj]][3] / cos( betaxy[hydi[hydj]] ); // PHASE 2 flow height
                        hydh3 = aw[hydi[hydj]][6] / cos( betaxy[hydi[hydj]] ); // PHASE 3 flow height
                        if ( hydh2 > sico.HFLOWMIN ) hydv2 = pow( pow( aw[hydi[hydj]][4], 2 ) + pow( aw[hydi[hydj]][5], 2 ), 0.5 ) / aw[hydi[hydj]][3];
                        else hydv2= 0; // PHASE 2 flow velocity
                        if ( hydh3 > sico.HFLOWMIN ) hydv3 = pow( pow( aw[hydi[hydj]][7], 2 ) + pow( aw[hydi[hydj]][8], 2 ), 0.5 ) / aw[hydi[hydj]][6];
                        else hydv3= 0; // PHASE 3 flow velocity
                        hyde = aw[hydi[hydj]][9]; // change of basal surface mixture or PHASE 1
                        hyde2 = aw[hydi[hydj]][10]; // change of basal surface PHASE 2
                        hyde3 = aw[hydi[hydj]][11]; // change of basal surface PHASE 3
                    }

                    hydq = 0; hydq2 = 0; hydq3 = 0; // resetting discharges

                    if ( hydalpha[hydj] == 0 || hydalpha[hydj] == sico.PI * 0.5 || hydalpha[hydj] == sico.PI || hydalpha[hydj] == 3 * sico.PI * 0.5 ) hydfalpha = 1;
                    else if ( fabs( 1 / sin( hydalpha[hydj] )) < fabs( 1 / cos( hydalpha[hydj] ))) hydfalpha = fabs( 1 / sin( hydalpha[hydj] ));
                    else hydfalpha = fabs( 1 / cos( hydalpha[hydj] )); // correction factor for profile direction

                    for ( hydk = 1; hydk <= hydp[0][hydj]; hydk++ ) { // loop over all cells of profile

                        alpha = falpha( betax[hydp[hydk][hydj]], betay[hydp[hydk][hydj]], sico ); // aspect
                        hydbeta = atan ( tan( betaxy[hydp[hydk][hydj]] ) * cos ( alpha - hydalpha[hydj] + sico.PI * 0.5 )); // corrected slope

                        hydfcorr = sico.CSZ * hydfalpha / cos( hydbeta ); // reference length for discharge

                        hydm0 = pow( pow(aw[hydp[hydk][hydj]][1], 2 ) + pow( aw[hydp[hydk][hydj]][2], 2), 0.5 );

                        if ( aw[hydp[hydk][hydj]][0] > sico.HFLOWMIN && hydm0 > 0 ) {

                            hydmx = aw[hydp[hydk][hydj]][1];
                            hydmy = aw[hydp[hydk][hydj]][2];

                            if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 ); else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                            hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                            hydm = hydm0 * cos( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                            hydq += hydm * hydfcorr; // updating mixture or PHASE 1 discharge
                        }

                        if ( sico.MODEL == 7 ) {

                            hydm0 = pow( pow(aw[hydp[hydk][hydj]][4], 2 ) + pow( aw[hydp[hydk][hydj]][5], 2), 0.5 );

                            if ( aw[hydp[hydk][hydj]][3] > sico.HFLOWMIN && hydm0 > 0 ) {

                                hydmx = aw[hydp[hydk][hydj]][4];
                                hydmy = aw[hydp[hydk][hydj]][5];

                                if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 ); else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                hydq2 += hydm * hydfcorr; // updating PHASE 2 discharge
                            }

                            hydm0 = pow( pow(aw[hydp[hydk][hydj]][7], 2 ) + pow( aw[hydp[hydk][hydj]][8], 2), 0.5 );

                            if ( aw[hydp[hydk][hydj]][6] > sico.HFLOWMIN && hydm0 > 0 ) {

                                hydmx = aw[hydp[hydk][hydj]][7];
                                hydmy = aw[hydp[hydk][hydj]][8];

                                if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 ); else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                hydq3 += hydm * hydfcorr; // updating PHASE 3 discharge
                            }
                        }
                    }

                    fprintf(f_hydinfo[hydj], "%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                        tsum, hydh, hydv, hyde, hydq, hydh2, hydv2, hyde2, hydq2, hydh3, hydv3, hyde3, hydq3); // writing hydrograph info to file

                    if ( hydq != 0 || hydq2 != 0 || hydq3 != 0 ) {

                        if ( qtinit[hydj] == sico.UNDEF ) qtinit[hydj] = tsum;

                        fprintf(f_hydtrans[hydj], "%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                        tsum - qtinit[hydj], hydq, hydv, hydq2, hydv2, hydq3, hydv3 ); // writing hydrograph info to file in input hydrograph format
                    }
                }
            }


// -- STOP --- Writing hydrograph infos to files ----------------------------------------------------------------


// -- START -- Display and files of status of simulation --------------------------------------------------------


            if ( sico.MODEL <= 3 ) { // one-phase models


               #ifdef WITHGRASS


                    printf("   %i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.*f\n", nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max,
                        prec_vol, vol_flow/1000, prec_ekin, ekin_flow/1000000); // display


               #else


                    printf("   %i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.*f", nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max,
                        prec_vol, vol_flow/1000, prec_ekin, ekin_flow/1000000);


                #endif


                fflush(stdout);

                fprintf(f_summary, "%i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.*f\n", nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max,
                    prec_vol, vol_flow, prec_ekin, ekin_flow); // summary file
                    
            } else if ( sico.MODEL == 7 ) { // multi-phase model


               #ifdef WITHGRASS


                    printf("   %i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.*f\t%.*f\t%.*f\n",
                        nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max, prec_hflow, hflow_max2, vflow_max2,
                        prec_hflow, hflow_max3, vflow_max3, prec_vol, vol_flow/1000, prec_vol, vol_flow2/1000, prec_vol, vol_flow3/1000, prec_ekin, ekin_flow/1000000); // display


                #else


                    printf("   %i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.*f\t%.*f\t%.*f",
                        nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max, prec_hflow, hflow_max2, vflow_max2,
                        prec_hflow, hflow_max3, vflow_max3, prec_vol, vol_flow/1000, prec_vol, vol_flow2/1000, prec_vol, vol_flow3/1000, prec_ekin, ekin_flow/1000000);


                #endif


                fflush(stdout);
                
                fprintf(f_summary, "%i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.*f\t%.*f\t%.*f\n",
                    nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max, prec_hflow, hflow_max2, vflow_max2,
                    prec_hflow, hflow_max3, vflow_max3, prec_vol, vol_flow, prec_vol, vol_flow2, prec_vol, vol_flow3, prec_ekin, ekin_flow); // summary file
            }

            fprintf(f_volumes, "%i\t%.1f", nout, tsum ); // volumes file
    
            for ( z=0; z<nzones; z++ ) {
    
                fprintf(f_volumes, "\t%.3f\t%.3f", vol_zone1[z], vol_czone1[z] );
                if ( sico.MODEL == 7 ) fprintf(f_volumes, "\t%.3f\t%.3f\t%.3f\t%.3f", vol_zone2[z], vol_czone2[z], vol_zone3[z], vol_czone3[z] );
            }
    
            fprintf(f_volumes, "\n");

            sprintf(path, "%s%snout%d.txt", outfiles, prefix, xint); // file of number of time steps, control for success, and volumes
            f_nout=fopen(path, "w");
            fprintf(f_nout, "%i\n%i\n%.2f\n%.0f\n%0f\n%0f\n%.0f\n%0f\n%0f\n", nout, csuccess, vmax, vol_entr, vol_entr2, vol_entr3, vol_edge, vol_edge2, vol_edge3 );
            fclose(f_nout);


// -- STOP --- Display and files of status of simulation --------------------------------------------------------


// -- START -- Preparing and writing output raster maps and velocity fields -------------------------------------


            #ifdef WITHGRASS


                for ( i=0; i<sico.IMAX; i++ ) {

                    if ( cdomain[i] != 0 ) { // if cell is not an edge cell, converting depths to heights

                        if ( sico.MODEL <= 3 ) {

                            hflowi = fconvout( i, aw, 1, 1, betaxy[i], sico );
                            hentri = fconvout( i, aw, 1, 2, betaxy[i], sico );
                            
                        } else if ( sico.MODEL == 7 ) {

                            hflowi = fconvout( i, aw, 1, 1, betaxy[i], sico );
                            hflowi2 = fconvout( i, aw, 2, 1, betaxy[i], sico );
                            hflowi3 = fconvout( i, aw, 3, 1, betaxy[i], sico );

                            hentri = fconvout( i, aw, 1, 2, betaxy[i], sico );
                            hentri2 = fconvout( i, aw, 2, 2, betaxy[i], sico );
                            hentri3 = fconvout( i, aw, 3, 2, betaxy[i], sico );
                        }
                        
                    } else { hflowi = 0; hentri = 0; if ( sico.MODEL == 7 ) { hflowi2 = 0; hentri2 = 0; hflowi3 = 0; hentri3 = 0; }}  // for edge cells, applying depths as heights

                    v[0] = hflowi; // mixture or PHASE 1 flow height
                    if ( aw[i][0] > sico.HFLOWMIN ) v[1] = aw[i][2] / aw[i][0]; else v[1] = 0;
                    if ( aw[i][0] > sico.HFLOWMIN ) v[2] = -aw[i][1] / aw[i][0]; else v[2] = 0; // mixture or PHASE 1 x and y velocities

                    if ( sico.MODEL <= 3 ) { // one-phase models

                        v[3] = hentri; // change of basal surface
                        v[4] = aw[i][4]; // velocity
                        v[5] = aw[i][5]; // flow kinetic energy
                        v[6] = aw[i][6]; // flow pressure
                        
                    } else if ( sico.MODEL == 7 ) { // multi-phase model

                        v[3] = hflowi2; // PHASE 2 flow height
                        if ( aw[i][3] > sico.HFLOWMIN ) v[4] = aw[i][5] / aw[i][3]; else v[4] = 0;
                        if ( aw[i][3] > sico.HFLOWMIN ) v[5] = -aw[i][4] / aw[i][3]; else v[5] = 0; // PHASE 2 x and y velocities
                        v[6] = hflowi3; // PHASE 3 flow height
                        if ( aw[i][6] > sico.HFLOWMIN ) v[7] = aw[i][5] / aw[i][6]; else v[7] = 0;
                        if ( aw[i][6] > sico.HFLOWMIN ) v[8] = -aw[i][4] / aw[i][6]; else v[8] = 0; // PHASE 3 x and y velocities
                        v[9] = hentri;
                        v[10] = hentri2;
                        v[11] = hentri3; // change of basal surface
                        v[12] = aw[i][12];
                        v[13] = aw[i][13];
                        v[14] = aw[i][14]; // flow velocities
                        v[15] = hflowi + hflowi2 + hflowi3; // total flow height
                        v[16] = aw[i][16];
                        v[17] = aw[i][17];
                        v[18] = aw[i][18];
                        v[19] = aw[i][19]; // flow kinetic energies
                        v[20] = aw[i][20];
                        v[21] = aw[i][21];
                        v[22] = aw[i][22];
                        v[23] = aw[i][23]; // flow pressures
                        v[24] = hentri + hentri2 + hentri3; // total change of basal surface
                    }

                    for ( k=0; k<nvect_red; k++ ) outv[px[i]][py[i]][k] = v[k];
                }


            #endif


            if ( sico.MULT == 0 ) { // for single model run

                if ( nout < 10 ) sprintf( madd, "000"); // fill string (for maintaining correct order in list of maps)
                else if ( nout < 100 ) sprintf( madd, "00");
                else if ( nout < 1000 ) sprintf( madd, "0");

                for ( k=0; k<nvect_red; k++ ) {
                
                    sprintf( mv, "%s%s%s%i", prefix, mv0[k], madd, nout ); // names of output raster maps

                    if ( sico.AFLAG == 1 || ( sico.MODEL <= 3 && ( k==0 || k==3 )) || ( sico.MODEL == 7 && ( k==0 || k==3 || k==6 || k==9 || k==10 || k==11 || k == 15 || k == 24 ))) {


                        #ifdef WITHGRASS


                            foutrast ( mv, outv, sico, k, tsum ); // if GRASS is used, writing GRASS raster maps


                        #endif


                        if ( sico.MODEL <= 3 ) foutasc ( aw, px, py, outmaps, mv, betaxy, k, sico ); // writing ascii raster maps (one-phase models)
                        else if ( sico.MODEL == 7 ) foutasc ( aw, px, py, outmaps, mv, betaxy, k, sico ); // writing ascii raster maps (multi-phase model)
                    }
                }

                if ( sico.MODEL <= 3 ) { // ascii raster maps of maximum flow height at time step
                    sprintf( mv, "%s%s%s%i", prefix, mv0[7], madd, nout );
                    foutasc ( aw, px, py, outmaps, mv, betaxy, 7, sico );
                    
                } else if ( sico.MODEL == 7 ) {
                    sprintf( mv, "%s%s%s%i", prefix, mv0[31], madd, nout );
                    foutasc ( aw, px, py, outmaps, mv, betaxy, 31, sico );
                }

                if ( nout == 1 && sico.MODEL <= 3 ) foutdircoord ( f_directions, 1, px, py, sico ); // file for display of flow vectors as arrows
                else if ( nout == 1 && sico.MODEL == 7 ) {
                    foutdircoord ( f_directions, 4, px, py, sico );
                    foutdircoord ( f_directions2, 5, px, py, sico );
                    foutdircoord ( f_directions3, 6, px, py, sico );
                } // writing coordinates to file

                if ( sico.MODEL <= 3 ) { for ( k=0; k<3; k++ ) foutdir ( f_directions, aw, k, 1, px, py, sico );

                } else if ( sico.MODEL == 7 ) {
                    for ( k=0; k<3; k++ ) foutdir ( f_directions, aw, k, 4, px, py, sico );
                    for ( k=3; k<6; k++ ) foutdir ( f_directions2, aw, k, 5, px, py, sico );
                    for ( k=6; k<9; k++ ) foutdir ( f_directions3, aw, k, 6, px, py, sico );
                } // writing parameters to file
            }

            if ( (int)( round( 1000 * tsum )) >= (int)( round( 1000 * tmax )) || ccontinue == 0 ) { // raster maps for last time step

                for ( k=0; k<nvect_red; k++ ) {
                
                    if (( sico.MODEL <= 3 && ( k == 0 || k == 3 )) || ( sico.MODEL == 7 && ( k==0 || k==3 || k==6 || k==9 || k==10 || k==11 || k==15 || k==24 ))) {

                        if ( sico.MULT == 0 ) sprintf( mv, "%s%s_fin", prefix, mv0[k] );
                        else sprintf( mv, "%s%s_fin%d", prefix, mv0[k], xint ); // names of output raster maps


                        #ifdef WITHGRASS


                            foutrast ( mv, outv, sico, k, tsum ); // if GRASS is used, writing GRASS raster maps


                        #endif


                        if ( sico.MODEL <= 3 ) foutasc ( aw, px, py, outmaps, mv, betaxy, k, sico ); // writing ascii raster maps (one-phase models)
                        else if ( sico.MODEL == 7 ) foutasc ( aw, px, py, outmaps, mv, betaxy, k, sico ); // writing ascii raster maps (multi-phase model)
                    }
                }
            }


// -- STOP --- Preparing and writing output raster maps and velocity fields -------------------------------------


            cflmax = 0; // resetting maximum cfl value
            nout += 1; // updating number of output time steps
        }

        nsum += 1; // updating total number of time steps


// *** End of loop over time steps ------------------------------------------------------------------------------


    }


// -- START -- Preparing and writing output raster maps of cumulative and maximum values ------------------------


    for ( k=nvect_red; k<nvect_all; k++ ) {

        if ( sico.AFLAG == 1 || ( sico.MODEL <= 3 && ( k==7 || k == 11 || k==16 || k==17 )) 
            || ( sico.MODEL == 7 && ( k==25 || k==27 || k==29 || k==31 || k==40 || k==41 || k==42 || k==43 || k==44 || k==45 || k==46 ))) {

            if ( sico.MULT == 0 ) sprintf( mv, "%s%s", prefix, mv0[k] );
            else sprintf( mv, "%s%s%d", prefix, mv0[k], xint ); // names of output raster maps


            #ifdef WITHGRASS


                for ( i=0; i<sico.IMAX; i++ ) {

                    if ( cdomain[i] != 0 ) { // if cell is not an edge cell, converting depths into heights

                        if ( sico.MODEL <= 3 && k == 7 ) pout = fconvout( i, aw, 1, 3, betaxy[i], sico );
                        else if ( sico.MODEL == 7 && k == 25 ) pout = fconvout( i, aw, 1, 3, betaxy[i], sico );
                        else if ( sico.MODEL == 7 && k == 27 ) pout = fconvout( i, aw, 2, 3, betaxy[i], sico );
                        else if ( sico.MODEL == 7 && k == 29 ) pout = fconvout( i, aw, 3, 3, betaxy[i], sico );
                        else if ( sico.MODEL == 7 && k == 31 ) pout = fconvout( i, aw, 4, 3, betaxy[i], sico );
                        else pout = aw[i][k];
                    }
                    else pout = aw[i][k]; // for edge cells, applying depths as heights

                    outv[px[i]][py[i]][k] = pout;
                }

                foutrast ( mv, outv, sico, k, tsum ); // writing GRASS raster maps


            #endif


            foutasc ( aw, px, py, outmaps, mv, betaxy, k, sico ); // writing ascii raster maps
        }
    }

    
// -- STOP --- Preparing and writing output raster maps of cumulative and maximum values ------------------------


// -- START -- Cleaning system ----------------------------------------------------------------------------------


    #ifdef WITHGRASS


        Segment_release(&seg_elev); // releasing segment data (if GRASS is used)
        if ( sico.RELM == 1 ) Segment_release(&seg_hrelease);
        if ( sico.RELV == 1 ) { Segment_release(&seg_vinx); Segment_release(&seg_viny); }
        if ( sico.ENTR == 1 ) Segment_release(&seg_hentrmax);
        if ( sico.MODEL == 7 && sico.RELM2 == 1 ) Segment_release(&seg_hrelease2);
        if ( sico.MODEL == 7 && sico.RELM3 == 1 ) Segment_release(&seg_hrelease3);
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) { Segment_release(&seg_vinx2); Segment_release(&seg_viny2); }
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) { Segment_release(&seg_vinx3); Segment_release(&seg_viny3); }
        if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) Segment_release(&seg_hentrmax2);
        if ( sico.MODEL == 7 && sico.ENTR3 == 1 ) Segment_release(&seg_hentrmax3);
        if ( sico.ZONES == 1 ) Segment_release(&seg_zones);
        if ( sico.CENTR == 1 ) Segment_release(&seg_centr);
        if ( sico.CVSHEAR == 1 ) Segment_release(&seg_cvshear);
        if ( sico.PHI == 1 ) Segment_release(&seg_phi);
        if ( sico.PHI2 == 1 ) Segment_release(&seg_phi2);
        if ( sico.PHI3 == 1 ) Segment_release(&seg_phi3);
        if ( sico.DELTAB == 1 ) Segment_release(&seg_deltab);
        if ( sico.TUFRI == 1 ) Segment_release(&seg_tufri);
        if ( sico.DELTA == 1 ) Segment_release(&seg_delta);
        if ( sico.DELTA2 == 1 ) Segment_release(&seg_delta2);
        if ( sico.DELTA3 == 1 ) Segment_release(&seg_delta3);
        if ( sico.NYSS == 1 ) Segment_release(&seg_nyss);
        if ( sico.NYFS == 1 ) Segment_release(&seg_nyfs);
        if ( sico.NYFF == 1 ) Segment_release(&seg_nyff);
        if ( sico.AMBDRAG == 1 ) Segment_release(&seg_ambdrag);
        if ( sico.FLUFRI == 1 ) Segment_release(&seg_flufri);
        if ( sico.TRANSSSFS == 1 ) Segment_release(&seg_transssfs);
        if ( sico.TRANSSSFF == 1 ) Segment_release(&seg_transssff);
        if ( sico.TRANSFSFF == 1 ) Segment_release(&seg_transfsff);
        if ( sico.TRELEASE == 1 ) Segment_release(&seg_trelease);
        if ( sico.TRELSTOP == 1 ) Segment_release(&seg_trelstop);
        if ( sico.STOPTIME == 1 ) Segment_release(&seg_stoptime);
        if ( sico.TSLIDE == 1 ) Segment_release(&seg_tslide);

        free( v ); free_dmatrix3(outv, sico.M, sico.N); // freeing memory


    #endif


    fclose ( f_summary );
    fclose ( f_volumes );
    if ( hydrograph == 1 ) {

        fclose ( f_hydout );
        for ( i = 0; i < hydnin + hydnout; i++ ) {
            if ( sico.MULT == 0 ) { 

                fprintf(f_hydtrans[i], "%.1f\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\n", tsum - qtinit[i] + tout );
                fprintf(f_hydtrans[i], "100000.0\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\n"); // finalizing hydrograph output files in input format
                fclose( f_hydinfo[i] ); fclose( f_hydtrans[i] );
            }
        }
    }

    if ( sico.MULT == 0 ) fclose ( f_directions );
    if ( sico.MULT == 0 && sico.MODEL == 7 ) { fclose ( f_directions2 ); fclose ( f_directions3 ); } // closing files

    free(sico.MAINMAPSET); free(mv); free(mv0); free(madd); free(prefix); free(outmaps); free(outfiles); free(in[0]); free(in); // freeing operational arrays

    free(flowpar); free(ib); free(ibasket[0]); free(ibasket); free(icheck[0]); free(icheck); free(cflux); free(cdomain); free(cdomain2); free(cstopped);
        // freeing flow parameter and control arrays

    free(elevname); free(pelev); free(pelev0); free(px); free(py); free(dx); free(dy); free(betax); free(betay); free(betaxy); // freeing terrain arrays
    
    free(aw[0]); free(awt[0]); free(af[0]); free(ag[0]); free(as[0]); free(ad[0]); free(asigma_x[0]); free(asigma_y[0]); 
    free(asigma_f[0]); free(asigma_g[0]); free(wintd[0]); free(wintdtest[0]); free(f[0]); free(g[0]); free(s[0]); free(wintelev[0]);
    free(aw); free(awt); free(af); free(ag); free(as); free(ad); free(asigma_x); free(asigma_y); 
    free(asigma_f); free(asigma_g); free(asigma_xelev); free(asigma_yelev); free(wintd); free(wintdtest); free(f); free(g); free(s);  free(wintelev);
    free_dmatrix3(winta, sico.IMAX, 4); free_dmatrix3(wintb, sico.IMAX, 6); free_dmatrix3(wintc, sico.IMAX, 4); free_dmatrix3(d, sico.IMAX, 6);
        // freeing state variable and numerical scheme arrays

    if ( hydrograph == 1 ) {
        free_dmatrix3(hydhyd, hydtmaxx+1, 5);
        free(hydp[0]), free(hydp);
        free(hydelev); free(hydalpha); free(hydx_metric); free(hydy_metric); free(hydl); free(hydtmax); free(hydi); free(hydx); free(hydy);
        for ( hydj=0; hydj<hydnin; hydj++ ) free(hydreleasename[hydj]); // freeing hydrograph arrays
    }
    
    if ( sico.RELM == 1 ) free(hreleasename); 
    if ( sico.RELV == 1 ) { free(vinxname); free(vinyname); }
    if ( sico.TRELEASE == 1 ) { free(treleasename); free( ptrelease ); } 
    if ( sico.TRELSTOP == 1 ) { free(trelstopname); free( ptrelstop ); }
    free(phrelease); 
    if ( sico.TRELEASE == 1 ) free(qhrelease); 
    free(pvinx); 
    free(pviny);

    if ( sico.MODEL == 7 ) {
    
        if ( sico.RELM2 == 1 ) free(hreleasename2); 
        if ( sico.RELV2 == 1 ) { free(vinxname2); free(vinyname2); }
        if ( sico.RELM3 == 1 ) free(hreleasename3); 
        if ( sico.RELV3 == 1 ) { free(vinxname3); free(vinyname3); }
        free(phrelease2); 
        if ( sico.TRELEASE == 1 ) free(qhrelease2); 
        free(pvinx2); 
        free(pviny2);
        free(phrelease3); 
        if ( sico.TRELEASE == 1 ) free(qhrelease3); 
        free(pvinx3); 
        free(pviny3); // freeing release arrays
    }

    free(phentrmax);
    if ( sico.ENTR == 1 ) free(hentrmaxname);
    if ( sico.ZONES == 1 ) free(zonesname);
    free(pzones);
    if ( sico.CENTR == 1 ) { free(centrname); free(pcentr); } 
    if ( sico.CVSHEAR == 1 ) { free(cvshearname); free(pcvshear); } 
    if ( sico.DELTAB == 1 ) { free(deltabname); free(pdeltab); }
    if ( sico.MODEL == 7 ) { free(phentrmax2); free(phentrmax3); 
    if ( sico.ENTR2 == 1 ) free(hentrmaxname2); 
    if ( sico.ENTR3 == 1 ) free(hentrmaxname3); } // freeing arrays for entrainment
    
    if ( sico.STOPTIME == 1 ) free(stoptimename); 
    free( pstoptime );
    if ( sico.TSLIDE == 1 ) free(tslidename);
    free( ptslide );
    free( pxslide );
    free(anx);
    free(anu); // freeing arrays for stopping and initial sliding

    if ( sico.PHI == 1 ) { free(phiname); free(pphi); }  
    if ( sico.PHI2 == 1 ) { free(phi2name); free(pphi2); } 
    if ( sico.PHI3 == 1 ) { free(phi3name); free(pphi3); } 
    if ( sico.FLUFRI == 1 ) { free(flufriname); free(pflufri); } 
    if ( sico.DELTA == 1 ) { free(deltaname); free(pdelta); } 
    if ( sico.DELTA2 == 1 ) { free(delta2name); free(pdelta2); } 
    if ( sico.DELTA3 == 1 ) { free(delta3name); free(pdelta3); } 
    if ( sico.TUFRI == 1 ) free(tufriname); // freeing friction arrays
    
    if ( sico.NYSS == 1 ) { free(nyssname); free(pnyss); } 
    if ( sico.NYFS == 1 ) { free(nyfsname); free(pnyfs); } 
    if ( sico.NYFF == 1 ) { free(nyffname); free(pnyff); } 
    if ( sico.AMBDRAG == 1 ) { free(ambdragname); free(pambdrag); } // freeing viscosity and ambient drag arrays
    
    if ( sico.TRANSSSFS == 1 ) free(transssfsname); 
    if ( sico.TRANSSSFF == 1 ) free(transssffname); 
    if ( sico.TRANSFSFF == 1 ) free(transfsffname);
        // freeing phase transformation arrays

    if ( frictiograph == 1 ) { free(frictioname); free(frifri[0]); free(frifri); } // freeing frictiograph arrays
    if ( transformograph == 1 ) { free(transformoname); free(tratra[0]); free(tratra); } // freeing transformograph arrays
  
    if ( sico.CURVCTRL > 2 ) { 
    
        free(cedge[0]); free(cready[0]); free(cedge); free(cready); free(cedge0); free(cneighbours);
        if ( sico.MODEL == 7 ) {
    
            free(cedge2[0]); free(cready2[0]); free(cedge2); free(cready2); free(cedge02); free(cneighbours2); 
            free(cedge3[0]); free(cready3[0]); free(cedge3); free(cready3); free(cedge03); free(cneighbours3); // freeing diffusion control arrays
        }
    }
    
    free( path ); free( wkdir );


// -- STOP --- Cleaning system ----------------------------------------------------------------------------------


    time_stop = clock(); // time at end of model execution
    time_elapsed = ( ( float ) ( time_stop - time_start ) ) / CLOCKS_PER_SEC; // time needed for model execution
    printf("Model execution completed in %.2f seconds.\n", time_elapsed);
    fflush(stdout);

    return 0;
}
