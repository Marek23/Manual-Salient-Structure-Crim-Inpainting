#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"

#define printfFnc(...) { mexPrintf(__VA_ARGS__); mexEvalString("drawnow;");}

struct Point {
    int x;
    int y;
    int salNum;
};

struct SalPoint {
    int x;
    int y;
    int filled;
};

float **qM, **pM, **ssd3, **subC; //p_s

float ***qI, ***pI; //p_s 3d

float **pmask, **pC, **ppom, **pe1mask; // nx + p_r

float ***PI; // nx + p_r 3d

float ***I;
int   ***SI; // nx 3d

float **e1mask, **e2mask, **mask, **C; // nx

struct Point *points; //

struct SalPoint **sals;
int *salLen;

float *mu, *mv, *ml, *PIgu, *PIgv, *D;

void initVars(int p_r, int p_s, int nx, int ny, int nz, int salNumb) {
    int iii,jjj;

    qM   = calloc(p_s, sizeof(float *));
    pM   = calloc(p_s, sizeof(float *));
    ssd3 = calloc(p_s, sizeof(float *));
    subC = calloc(p_s, sizeof(float *));
    for(iii = 0; iii < p_s; iii++) {
        qM  [iii] = calloc(p_s, sizeof(float));
        pM  [iii] = calloc(p_s, sizeof(float));
        ssd3[iii] = calloc(p_s, sizeof(float));
        subC[iii] = calloc(p_s, sizeof(float));
    }

    qI = calloc(p_s, sizeof(float **));
    pI = calloc(p_s, sizeof(float **));
    for(iii = 0; iii < p_s; iii++) {
        qI[iii] = calloc(p_s, sizeof(float *));
        pI[iii] = calloc(p_s, sizeof(float *));
        for(jjj = 0; jjj < p_s; jjj++) {
            qI[iii][jjj] = calloc(nz, sizeof(float));
            pI[iii][jjj] = calloc(nz, sizeof(float));
        }
    }

    pmask    = calloc(2*p_r + nx, sizeof(float  *));
    pC       = calloc(2*p_r + nx, sizeof(float  *));
    ppom     = calloc(2*p_r + nx, sizeof(float  *));
    pe1mask  = calloc(2*p_r + nx, sizeof(float  *));
    for(iii = 0; iii < (2*p_r + nx); iii++) {
        pmask  [iii] = calloc((2*p_r + ny), sizeof(float  ));
        pC     [iii] = calloc((2*p_r + ny), sizeof(float  ));
        ppom   [iii] = calloc((2*p_r + ny), sizeof(float  ));
        pe1mask[iii] = calloc((2*p_r + ny), sizeof(float  ));
    }

    PI    = calloc(2*p_r + nx, sizeof(float **));
    for(iii = 0; iii < (2*p_r + nx); iii++) {
        PI   [iii] = calloc((2*p_r + ny), sizeof(float *));
        for(jjj = 0; jjj <  (2*p_r + ny); jjj++) {
            PI[iii][jjj] = calloc(nz, sizeof(float));
        }
    }

    I  = calloc(nx, sizeof(float **));
    SI = calloc(nx, sizeof(int **));
    for(iii = 0; iii < nx; iii++) {
        I [iii] = calloc(ny, sizeof(float *));
        SI[iii] = calloc(ny, sizeof(int *));
        for(jjj = 0; jjj < ny; jjj++) {
            I [iii][jjj] = calloc(nz, sizeof(float));
            SI[iii][jjj] = calloc(nz, sizeof(int));
        }
    }

    e1mask = calloc(nx, sizeof(float *));
    e2mask = calloc(nx, sizeof(float *));
    mask   = calloc(nx, sizeof(float *));
    C      = calloc(nx, sizeof(float *));
    for(iii = 0; iii < nx; iii++) {
        e1mask[iii] = calloc(ny, sizeof(float));
        e2mask[iii] = calloc(ny, sizeof(float));
        mask[iii]   = calloc(ny, sizeof(float));
        C[iii]      = calloc(ny, sizeof(float));
    }

    points = calloc(nx*ny, sizeof(struct Point));

    sals   = calloc(salNumb, sizeof(struct SalPoint *));
    salLen = calloc(salNumb, sizeof(int));
    for(iii = 0; iii < salNumb; iii++) {
        sals[iii] = calloc(4*nx, sizeof(struct SalPoint));
    }

    mu     = calloc(nx*ny, sizeof(float));
    mv     = calloc(nx*ny, sizeof(float));
    ml     = calloc(nx*ny, sizeof(float));
    PIgu   = calloc(nx*ny, sizeof(float));
    PIgv   = calloc(nx*ny, sizeof(float));
    D      = calloc(nx*ny, sizeof(float));
}

void clearVars(int p_r, int p_s, int nx, int ny, int nz, int salNumb) {
    //float **qM, **pM, **ssd3, **subC; //p_s
    int iii,jjj;
    for(iii = 0; iii < 2*p_r+1; iii++) {
        free(qM[iii]);
        free(pM[iii]);
        free(ssd3[iii]);
        free(subC[iii]);
    }
    free(qM);
    free(pM);
    free(ssd3);
    free(subC);

    //float ***qI, ***pI; //p_s 3d
    for(iii = 0; iii < 2*p_r+1; iii++) {
        for(jjj = 0; jjj < 2*p_r+1; jjj++) {
            free(qI[iii][jjj]);
            free(pI[iii][jjj]);
        }
        free(qI[iii]);
        free(pI[iii]);
    }
    free(qI);
    free(pI);

    //float **pmask, **pC, **ppom; // nx + p_r
    for(iii = 0; iii < (2*p_r + nx); iii++) {
        free(pmask[iii]);
        free(pC[iii]);
        free(ppom[iii]);
        free(pe1mask[iii]);
    }
    free(pmask);
    free(pC);
    free(ppom);
    free(pe1mask);

    // float ***PI; // nx + p_r 3d
    for(iii = 0; iii < (2*p_r + nx); iii++) {
        for(jjj = 0; jjj < (2*p_r + ny); jjj++) {
            free(PI[iii][jjj]);
        }
        free(PI[iii]);
    }
    free(PI);

    // float ***I, ***SI; // nx 3d
    for(iii = 0; iii < nx; iii++) {
        for(jjj = 0; jjj < ny; jjj++) {
            free(I[iii][jjj]);
            free(SI[iii][jjj]);
        }
        free(I[iii]);
        free(SI[iii]);
    }
    free(I);
    free(SI);

    //float **e1mask, **e2mask, **mask, **c; // nx
    for(iii = 0; iii < nx; iii++) {
        free(e1mask[iii]);
        free(e2mask[iii]);
        free(mask[iii]);
        free(C[iii]);
    }
    free(e1mask);
    free(e2mask);
    free(mask);
    free(C);

    for(iii = 0; iii < salNumb; iii++) {
        free(sals[iii]);
    }
    free(sals);
    free(salLen);

    //struct Point *points; //
    free(points);
    free(mu);
    free(mv);
    free(ml);
    free(PIgu);
    free(PIgv);
    free(D);
}

void findAndDelete(int salNumb, int _x, int _y)
{
    int i;
    for(i=0; i < salLen[salNumb]; i++)
    {
        if(sals[salNumb][i].x == _x && sals[salNumb][i].y == _y)
        {
            sals[salNumb][i].filled = 1;
            return;
        }
    }
}

int existInSals(int salAmount, int _x, int _y)
{
    int i, salNumb;
    for(salNumb=0; salNumb<salAmount; salNumb++)
    {
        for(i=0; i < salLen[salNumb]; i++)
        {
            if(sals[salNumb][i].x == _x && sals[salNumb][i].y == _y && sals[salNumb][i].filled < 1)
            {
                return salNumb;
            }
        }
    }
    return -1;
}

int allFilled(int salNumb)
{
    int i,j;
    for(i=0; i < salNumb; i++)
    {
        for(j=0; j < salLen[i]; j++)
        {
            if(sals[i][j].filled < 1)
            {
                return 0;
            }
        }

    }
    return 1;
}

void padarray2d(
    int m,
    int n,
    int t_r,
    float **in,
    float **new_in)
    {
    int new_m = m+(2*t_r);
    int new_n = n+(2*t_r);

    int i;
    int j;

    int old_i;
    int old_j;

    for(old_i=0;old_i<m;old_i++)
    {
        for(old_j=0;old_j<n;old_j++)
        {
            i=old_i+t_r;
            j=old_j+t_r;

            new_in[i][j]=in[old_i][old_j];
        }
    }

    for(i=0;i<t_r;i++)
    {
        for(j=0;j<new_n;j++)
        {
            new_in[i][j]=new_in[2*t_r-i][j];
            new_in[new_m-t_r+i][j]=new_in[new_m-t_r-2-i][j];
        }
    }

    for(i=0;i<new_m;i++)
    {
        for(j=0;j<t_r;j++)
        {
            new_in[i][j]=new_in[i][2*t_r-j];
            new_in[i][new_n-t_r+j]=new_in[i][new_n-t_r-2-j];
        }
    }
    }

void imerode(int nx, int ny, float **mask, float **out, int times)
{
    int i, j, i0, j0, t, ii, jj, pres;
    for (t = 0; t < times; t++)
    {
        padarray2d(nx,ny,1,mask,ppom);
        for (i = 0; i < nx; i++)
        {
            i0 = i+1;
            for (j = 0; j < ny; j++)
            {
                j0        = j+1;
                pres      = 0;
                out[i][j] = mask[i][j];
                for (ii = (i0-1); ii < (i0+2); ii++)
                {
                    for (jj = (j0-1); jj < (j0+2); jj++)
                    {
                        if (pres == 0 && ppom[ii][jj] < 1)
                            pres = 1;
                    }
                }
                if (pres == 1)
                    out[i][j] = 0;
            }
        }
    }
}

void imdilate(int nx, int ny, float **mask, float **out, int times)
{
    int i, j, i0, j0, t, ii, jj, pres;
    for (t = 0; t < times; t++)
    {
        padarray2d(nx,ny,1,mask,ppom);
        for (i = 0; i < nx; i++)
        {
            i0 = i+1;
            for (j = 0; j < ny; j++)
            {
                j0        = j+1;
                pres      = 0;
                out[i][j] = mask[i][j];
                for (ii = (i0-1); ii < (i0+2); ii++)
                {
                    for (jj = (j0-1); jj < (j0+2); jj++)
                    {
                        if (pres == 0 && ppom[ii][jj] > 0)
                            pres = 1;
                    }
                }
                if (pres == 1)
                    out[i][j] = 1;
            }
        }
    }
}

float rgb2gray(int i, int j, float ***in)
{
    return 0.3*in[i][j][0] + 0.59*in[i][j][1] + 0.11*in[i][j][2];
}

void showMatrix2D(
    int m,
    int n,
    float **w)
    {
    int i;
    int j;

    for(i=0; i<m; i++)
    {
        for(j=0;j<n;j++)
        {
            printfFnc("%.1f ", w[i][j]);
        }
    printf("\n");
    }
    }

void showMatrix3D(
    int m,
    int n,
    int c,
    float ***w)
    {
    int i;
    int j;
    int k;
    for(k=0; k<c; k++)
    {
        for(i=0; i<m; i++)
        {
            for(j=0;j<n;j++)
            {
                printf("%.1f ", w[i][j][k]);
            }
        printf("\n");
        }
    printf("\n \n");
    }
    }

void padarray3d(
    int m,
    int n,
    int c,
    int t_r,
    float ***in,
    float ***new_in)
    {
    int new_m = m+(2*t_r);
    int new_n = n+(2*t_r);

    int i;
    int j;
    int k;
    int old_i;
    int old_j;

    for(old_i=0;old_i<m;old_i++)
    {
        for(old_j=0;old_j<n;old_j++)
        {
            i=old_i+t_r;
            j=old_j+t_r;
            for(k=0;k<c;k++)
            {
                new_in[i][j][k]=in[old_i][old_j][k];
            }
        }
    }

    for(i=0;i<t_r;i++)
    {
        for(j=0;j<new_n;j++)
        {
            for(k=0;k<c;k++)
            {
                new_in[i][j][k]=new_in[2*t_r-i][j][k];
                new_in[new_m-t_r+i][j][k]=new_in[new_m-t_r-2-i][j][k];
            }
        }
    }

    for(i=0;i<new_m;i++)
    {
        for(j=0;j<t_r;j++)
        {
            for(k=0;k<c;k++)
            {
                new_in[i][j][k]=new_in[i][2*t_r-j][k];
                new_in[i][new_n-t_r+j][k]=new_in[i][new_n-t_r-2-j][k];
            }
        }
    }
    }

void subMatrix3D(
    int p_r,
    int m,
    int n,
    int c,
    int k,
    float ***in,
    int i0,
    int j0,
    float **out)
    {
    int i ,j;
    for(i=0;i<(2*p_r+1);i++)
    {
        for(j=0;j<(2*p_r+1);j++)
        {
            out[i][j]=in[i0-p_r+i][j0-p_r+j][k];
        }
    }
    }

void sub3Matrix3D(
    int p_r,
    int m,
    int n,
    int c,
    float ***in,
    int i0,
    int j0,
    float ***out)
    {
    int i ,j, k;
    for(i=0;i<(2*p_r+1);i++)
    {
        for(j=0;j<(2*p_r+1);j++)
        {
            for(k=0;k<c;k++)
            {
                out[i][j][k]=in[i0-p_r+i][j0-p_r+j][k];
            }
        }
    }
    }

void subMatrix(
    int p_r,
    int m,
    int n,
    float **in,
    int i0,
    int j0,
    float **out)
    {
    int i ,j;
    for(i=0;i<(2*p_r+1);i++)
    {
        for(j=0;j<(2*p_r+1);j++)
        {
            out[i][j]=in[i0-p_r+i][j0-p_r+j];
        }
    }
    }

void prodMatrix(
    int m,
    int n,
    float **in1,
    float **in2,
    float **prod)
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            prod[i][j]=in1[i][j]*in2[i][j];
        }
    }
    }

void prod3Matrix(
    int m,
    int n,
    float **in1,
    float **in2,
    float **in3,
    float **prod)
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            prod[i][j]=in1[i][j]*in2[i][j]*in3[i][j];
        }
    }
    }

void subt2Matrix(
    int m,
    int n,
    float **from,
    float **q,
    float **out)
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            out[i][j]=(from[i][j]-q[i][j])*(from[i][j]-q[i][j]);
        }
    }
    }

float sumMatrix(
    int m,
    int n,
    float **in)
    {
    int i ,j;
    float sum=0;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            sum+=in[i][j];
        }
    }
    return sum;
    }

float sumsqr(
    float a,
    float b,
    float c)
    {
    return sqrt(a*a + b*b + c*c);
    }

int any(
    int m,
    int n,
    float **in)
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(in[i][j]>0)
            {
                return 1;
            }
        }
    }
    return 0;
    }

int allOne(
    int m,
    int n,
    float **in)
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(in[i][j]<1)
            {
                return 0;
            }
        }
    }
    return 1;
    }

int allNonZero(
    int m,
    int n,
    float **in)
    {
    int i,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(in[i][j] > 0)
            {
                return 0;
            }
        }
    }
    return 1;
    }

void updatePhi(
    int M,
    int N,
    int c,
    float ***u,
    float ** phi)
    {
        int i,j;
        for(i=0;i<M;i++)
        {
            for(j=0;j<N;j++)
            {
                if(u[i][j][0]!=0   &&
                   u[i][j][1]!=255 &&
                   u[i][j][2]!=0   &&
                   phi[i][j]  ==0)
                {
                    phi[i][j]=1;
                }
            }
        }
    }

void subWeight(
    int m,
    int n,
    int c,
    int s_s,
    float *****w,
    int iw,
    int jw,
    int kw,
    float **out)
    {
        int i,j;
        for(i=0;i<s_s;i++)
        {
            for(j=0;j<s_s;j++)
            {
                out[i][j]=w[i][j][iw][jw][kw];
            }
        }
    }

void mexFunction(int numOut, mxArray *pmxOut[],
                 int numIn, const mxArray *pmxIn[])
{
    printfFnc("Wejście w algorytm. \n");
    int nx               = mxGetScalar(pmxIn[0]);
    int ny               = mxGetScalar(pmxIn[1]);
    int nz               = mxGetScalar(pmxIn[2]);
    double *I_i          = mxGetPr(pmxIn[3]);
    double *SI_i         = mxGetPr(pmxIn[4]);
    double *mask_i       = mxGetPr(pmxIn[5]);
    double *C_i          = mxGetPr(pmxIn[6]);

    int p_r              = mxGetScalar(pmxIn[7]);
    int s_r              = mxGetScalar(pmxIn[8]);
    float alfa           = mxGetScalar(pmxIn[9]);
    int salAmount        = mxGetScalar(pmxIn[10]);

    int x,y,min_x,max_x,x0,y0,min_y,max_y,i,j,ii,jj,it,in4d,in2d,k,l,c_it,p,salNumb;
    int pnx = 2*p_r + nx;
    int pny = 2*p_r + ny;
    
    printfFnc("nx %d ny %d nz %d p_r %d s_r %d alfa %f salAmount %d \n", nx, ny, nz, p_r, s_r, alfa, salAmount);

    int p_s = 2*p_r+1;
    initVars(p_r,p_s,nx,ny,nz,salAmount);

    printfFnc("Koniec inicjalizacji zmiennych. \n");
    c_it=0;
    it  =0;
    while(c_it<nz)
    {
        while(it<(c_it+1)*nx*ny)
        {
            in2d=0;
            while(in2d<nx*ny)
            {
                I [in2d%nx][in2d/nx][c_it] = (float)I_i [it];
                SI[in2d%nx][in2d/nx][c_it] = (int)SI_i[it++];
                in2d++;
            }
        }
        c_it++;
    }

    for (i = 0; i < nx*ny; i++)
    {
        mask[i%nx][i/nx] = (float)mask_i[i];
        C   [i%nx][i/nx] = (float)C_i   [i];
    }

    padarray2d(nx,ny,p_r,C,pC);

    // Salient structures points
    for (i=0; i<nx; i++)
    {
        for (j=0; j<ny; j++)
        {
            if(SI[i][j][0] == 1
            && SI[i][j][1] == 0
            && SI[i][j][2] == 0)
            {
                sals[0][salLen[0]].x      = i;
                sals[0][salLen[0]].y      = j;
                sals[0][salLen[0]].filled = 0;
                if (mask[i][j] == 1) {
                    sals[0][salLen[0]].filled = 1;
                }
                salLen[0]++;
            }

            if(SI[i][j][0] == 0
            && SI[i][j][1] == 1
            && SI[i][j][2] == 0)
            {
                sals[1][salLen[1]].x      = i;
                sals[1][salLen[1]].y      = j;
                sals[1][salLen[1]].filled = 0;
                if (mask[i][j] == 1) {
                    sals[1][salLen[1]].filled = 1;
                }
                salLen[1]++;
            }

            if(SI[i][j][0] == 0
            && SI[i][j][1] == 0
            && SI[i][j][2] == 1)
            {
                sals[2][salLen[2]].x      = i;
                sals[2][salLen[2]].y      = j;
                sals[2][salLen[2]].filled = 0;
                if (mask[i][j] == 1) {
                    sals[2][salLen[2]].filled = 1;
                }
                salLen[2]++;
            }

            if(SI[i][j][0] == 0
            && SI[i][j][1] == 1
            && SI[i][j][2] == 1)
            {
                sals[3][salLen[3]].x      = i;
                sals[3][salLen[3]].y      = j;
                sals[3][salLen[3]].filled = 0;
                if (mask[i][j] == 1) {
                    sals[3][salLen[3]].filled = 1;
                }
                salLen[3]++;
            }

            if(SI[i][j][0] == 1
            && SI[i][j][1] == 1
            && SI[i][j][2] == 0)
            {
                sals[4][salLen[4]].x      = i;
                sals[4][salLen[4]].y      = j;
                sals[4][salLen[4]].filled = 0;
                if (mask[i][j] == 1) {
                    sals[4][salLen[4]].filled = 1;
                }
                salLen[4]++;
            }

            if(SI[i][j][0] == 1
            && SI[i][j][1] == 0
            && SI[i][j][2] == 1)
            {
                sals[5][salLen[5]].x      = i;
                sals[5][salLen[5]].y      = j;
                sals[5][salLen[5]].filled = 0;
                if (mask[i][j] == 1) {
                    sals[5][salLen[5]].filled = 1;
                }
                salLen[5]++;
            }
        }
    }

    printfFnc("SalLen[0] %d SalLen[1] %d SalLen[2] %d. \n", salLen[0], salLen[1], salLen[2]);

    printfFnc("Koniec przepisania zmiennych. \n");
    int prev     = -1;
    int done     = -1;
    int salDone  = -1;
    int giveInfo = 0;

    //SAL INPAINTING
    printfFnc("Sal inpainting. \n");
    while (salDone != 1 && done != 1 && prev != sumMatrix(nx,ny,mask))
    {
        // if (giveInfo == 10)
        // {
        printfFnc("Zaawansowanie: %f \n", sumMatrix(nx,ny,mask)/nx/ny*100);
        // giveInfo = 0;
        // }
        prev = sumMatrix(nx,ny,mask);
        padarray2d(nx,ny,p_r,mask,pmask);
        padarray3d(nx,ny,nz,p_r,I,PI);
        imdilate(nx,ny,mask,e1mask,1);

        int len = 0;

        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                if(e1mask[i][j] != mask[i][j])
                {
                    if(i > 0 && i < (nx-1)
                    && j > 0 && j < (ny-1))
                    {
                        int pom = existInSals(salAmount,i,j);
                        if( pom > -1 )
                        {
                            points[len].x      = i;
                            points[len].y      = j;
                            points[len].salNum = pom;
                            len++;
                        }
                    }
                }
            }
        }

        //printfFnc("Wyznaczone pointy: %d. \n", len);
        
        for (i=0; i<len; i++)
        {
            sub3Matrix3D(p_r,pnx,pny,nz,PI,points[i].x+p_r,points[i].y+p_r,pI);
            subMatrix(p_r,pnx,pny,pmask,points[i].x+p_r,points[i].y+p_r,pM);

            //showMatrix3D(p_s,p_s,nz,pI);
            //showMatrix2D(p_s,p_s,pM);
            float min_sum = 99999999;
            int qx        = -1;
            int qy        = -1;
            for(j=0; j<salLen[points[i].salNum]; j++)
            {
                x0 = sals[points[i].salNum][j].x + p_r;
                y0 = sals[points[i].salNum][j].y + p_r;
                subMatrix(p_r,pnx,pny,pmask,x0,y0,qM);
                if(sals[points[i].salNum][j].filled==1 && allOne(p_s,p_s,qM))
                {
                    sub3Matrix3D(p_r,pnx,pny,nz,PI,x0,y0,qI);

                    //prepare ssd3
                    for(ii=0; ii < p_s; ii++)
                    {
                        for(jj=0; jj < p_s; jj++)
                        {
                            ssd3[ii][jj] = 0;
                        }
                    }
                    for(ii=0; ii < p_s; ii++)
                    {
                        for(jj=0; jj < p_s; jj++)
                        {
                            if(pM[ii][jj] == 1)
                            {
                                ssd3[ii][jj] = sqrt(
                                    ((pI[ii][jj][0] - qI[ii][jj][0])*(pI[ii][jj][0] - qI[ii][jj][0])) +
                                    ((pI[ii][jj][1] - qI[ii][jj][1])*(pI[ii][jj][1] - qI[ii][jj][1])) +
                                    ((pI[ii][jj][2] - qI[ii][jj][2])*(pI[ii][jj][2] - qI[ii][jj][2])));
                            }
                        }
                    }
                    float ssd3s = sumMatrix(p_s,p_s,ssd3);
                    if (min_sum > ssd3s)
                    {
                        min_sum = ssd3s;
                        qx      = sals[points[i].salNum][j].x;
                        qy      = sals[points[i].salNum][j].y;
                    }
                }
            }

            if (qx != -1 && qy != -1)
            {
                //printfFnc("qx %d qy %d\n",qx,qy);
                //printfFnc("qI\n");
                sub3Matrix3D(p_r,pnx,pny,nz,PI,qx+p_r,qy+p_r,qI);
                for(ii=0; ii < p_s; ii++)
                {
                    for(jj=0; jj < p_s; jj++)
                    {
                        if(points[i].x-p_r+ii > -1
                        && points[i].x-p_r+ii < nx
                        && points[i].y-p_r+jj > -1
                        && points[i].y-p_r+jj < ny
                        && mask[points[i].x+ii-p_r][points[i].y+jj-p_r] == 0)
                        {
                            I[   points[i].x-p_r+ii][points[i].y-p_r+jj][0] = qI[ii][jj][0];
                            I[   points[i].x-p_r+ii][points[i].y-p_r+jj][1] = qI[ii][jj][1];
                            I[   points[i].x-p_r+ii][points[i].y-p_r+jj][2] = qI[ii][jj][2];
                            mask[points[i].x-p_r+ii][points[i].y-p_r+jj]    = 1;
                            findAndDelete(points[i].salNum,points[i].x+ii-p_r,points[i].y+jj-p_r);
                        }
                    }
                }
            } else {
                printfFnc("Brak wzorca. \n");
                done = 1;
            }
        }
        if (allFilled(salAmount) == 1)
        {

            done = 1;
        }
    }

    done = 0;
    //printfFnc("allNonZero %d prev %d sumMatrix %f done %d", allNonZero(nx,ny,mask), prev, sumMatrix(nx,ny,mask), done);

    //CRIM INPAINTING
    while (allOne(nx,ny,mask) == 0 && prev != sumMatrix(nx,ny,mask) && done != 1)
    {
        if (giveInfo == 10)
        {
            printfFnc("Zaawansowanie: %f \n", sumMatrix(nx,ny,mask)/nx/ny*100);
            giveInfo = 0;
        }
        giveInfo++;
        prev = sumMatrix(nx,ny,mask);
        padarray2d(nx,ny,p_r,mask,pmask);
        padarray3d(nx,ny,nz,p_r,I,PI);
        //printfFnc("padArray. \n");
        imerode(nx,ny,  mask,e1mask,1);
        imerode(nx,ny,e1mask,e2mask,1);
        padarray2d(nx,ny,p_r,e1mask,pe1mask);
        //printfFnc("MASK \n");
        //showMatrix2D(nx,ny,mask);
        //printfFnc("E1MASK \n");
        //showMatrix2D(nx,ny,e1mask);
        //printfFnc("E2MASK \n");
        //showMatrix2D(nx,ny,e2mask);
        //printfFnc("Imerode. \n");
        int len = 0;

        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                if(e1mask[i][j] != e2mask[i][j])
                {
                    if(i > 0 && i < (nx-1)
                    && j > 0 && j < (ny-1))
                    {
                        points[len].x = i;
                        points[len].y = j;
                        len++;
                    }

                }
            }
        }
        //printfFnc("Pointy: %d. \n", len);
        int maxX     = -1;
        int maxY     = -1;
        int maxPrior = -1;
        int init     = 0;
        for (p = 0; p < len; p++)
        {
            x  = points[p].x;
            y  = points[p].y;
            x0 = x + p_r;
            y0 = y + p_r;
            //printfFnc("Point x0: %d y0: %d x: %d y: %d. pnx %d pny %d nx %d ny %d\n", x0, y0, x, y, pnx, pny, nx, ny);
            mv[p] = (
                (0.3*(pe1mask[x0+1][y0-1]-pe1mask[x0-1][y0-1])/2) +
                (0.4*(pe1mask[x0+1][y0  ]-pe1mask[x0-1][y0  ])/2) +
                (0.3*(pe1mask[x0+1][y0+1]-pe1mask[x0-1][y0+1])/2)
            );

            //printfFnc("mv \n");
            //printfFnc("[x0-1] %d [y0-1] %d \n", x0-1, y0-1);
            //printfFnc("[x0  ] %d [y0  ] %d \n", x0  , y0  );
            //printfFnc("[x0+1] %d [y0+1] %d \n", x0+1, y0+1);

            // printfFnc("pe1mask[x0-1][y0+1]: %f \n", pe1mask[x0-1][y0+1]);
            // printfFnc("pe1mask[x0-1][y0-1]: %f \n", pe1mask[x0-1][y0-1]);
            // printfFnc("pe1mask[x0  ][y0+1]: %f \n", pe1mask[x0  ][y0+1]);

            // printfFnc("pe1mask[x0  ][y0-1]: %f \n", pe1mask[x0  ][y0-1]);
            // printfFnc("pe1mask[x0  ][y0+1]: %f \n", pe1mask[x0  ][y0+1]);
            // printfFnc("pe1mask[x0+1][y0+1]: %f \n", pe1mask[x0+1][y0+1]);
            // printfFnc("pe1mask[x0+1][y0-1]: %f \n", pe1mask[x0+1][y0-1]);

            mu[p] = (
                (0.3*(pe1mask[x0-1][y0+1]-pe1mask[x0-1][y0-1])/2) +
                (0.4*(pe1mask[x0  ][y0+1]-pe1mask[x0  ][y0-1])/2) +
                (0.3*(pe1mask[x0+1][y0+1]-pe1mask[x0+1][y0-1])/2)
            )*(-1);

            //printfFnc("mu \n");
            if ((mu[p]*mu[p] + mv[p]*mv[p]) > 0)
            {
                ml[p] = sqrt(mu[p]*mu[p] + mv[p]*mv[p]);
            } else {
                ml[p] = 0;
            }
            
            if (ml[p] > 0)
            {
                mu[p] = mu[p]/ml[p];
                mv[p] = mv[p]/ml[p];
            } else {
                mu[p] = 0.0;
                mv[p] = 0.0;
            }

            PIgv[p] = (
                (0.3*(rgb2gray(x0+1,y0-1,PI)-rgb2gray(x0-1,y0-1,PI))/2) +
                (0.4*(rgb2gray(x0+1,y0  ,PI)-rgb2gray(x0-1,y0  ,PI))/2) +
                (0.3*(rgb2gray(x0+1,y0+1,PI)-rgb2gray(x0-1,y0+1,PI))/2)
            );
            //printfFnc("PIgv \n");
            PIgu[p] = (
                (0.3*(rgb2gray(x0-1,y0+1,PI)-rgb2gray(x0-1,y0-1,PI))/2) +
                (0.4*(rgb2gray(x0  ,y0+1,PI)-rgb2gray(x0  ,y0-1,PI))/2) +
                (0.3*(rgb2gray(x0+1,y0+1,PI)-rgb2gray(x0+1,y0-1,PI))/2)
            );
            //printfFnc("mv: %f mu: %f ml: %f PIgv: %f PIgu: %f \n", mv[p], mu[p], ml[p], PIgv[p], PIgu[p]);

            if ( ( (PIgu[p]*mu[p])*(PIgu[p]*mu[p]) + (PIgv[p]*mv[p])*(PIgv[p]*mv[p]) ) > 0)
            {
                D[p] = sqrt( (PIgu[p]*mu[p])*(PIgu[p]*mu[p]) + (PIgv[p]*mv[p])*(PIgv[p]*mv[p]) ) / alfa;
            } else {
                D[p] = 0;
            }

            //printfFnc("D[p] : %f \n",D[p]);
            subMatrix(p_r,nx,ny,pC,x0,y0,subC);
            //printfFnc("Po subC \n");
            if (maxPrior < sumMatrix(p_s,p_s,subC)*D[p])
            {
                init = 1;
                maxPrior = sumMatrix(p_s,p_s,subC)*D[p];
                maxX     = points[p].x;
                maxY     = points[p].y;
            }
        }

        if (init == 1)
        {
            // printfFnc("maxX %d maxY %d \n", maxX, maxY);
            if (s_r > 9000) {
                min_x = 0;
                min_y = 0;
                max_x = nx;
                max_y = ny;
            } else {

                if ( maxX-s_r < 0)
                {
                    min_x = 0;
                } else {
                    min_x = maxX-s_r;
                }

                if ( maxY-s_r < 0)
                {
                    min_y = 0;
                } else {
                    min_y = maxY-s_r;
                }

                if ( maxX+s_r > nx )
                {
                    max_x = nx;
                } else {
                    max_x = maxX+s_r;
                }

                if ( maxY+s_r > ny)
                {
                    max_y = ny;
                } else {
                    max_y = maxY+s_r;
                }
            }

            //printfFnc("min_x %d max_x %d min_y %d max_y %d \n",min_x, max_x, min_y, max_y);

            sub3Matrix3D(p_r,pnx,pny,nz,PI,maxX+p_r,maxY+p_r,pI);
            //printfFnc("Pobrane pI\n");
            subMatrix(p_r,pnx,pny,pmask,maxX+p_r,maxY+p_r,pM);
            //printfFnc("Pobrane pM\n");
            float min_sum = 99999999;
            int qx        = -1;
            int qy        = -1;
            for(x=min_x; x < max_x; x++)
            {
                x0 = x + p_r;
                for(y=min_y; y < max_y; y++)
                {
                    y0 = y + p_r;
                    subMatrix(p_r,pnx,pny,pmask,x0,y0,qM);
                    // printfFnc("Pobrane qM\n");
                    if (allOne(p_s,p_s,qM))
                    {
                        sub3Matrix3D(p_r,pnx,pny,nz,PI,x0,y0,qI);
                        // printfFnc("Pobrane qI\n");
                        for(i=0; i < p_s; i++)
                        {
                            for(j=0; j < p_s; j++)
                            {
                                if(pM[i][j] == 1)
                                {
                                    ssd3[i][j] = sqrt(
                                        ((pI[i][j][0] - qI[i][j][0])*(pI[i][j][0] - qI[i][j][0])) +
                                        ((pI[i][j][1] - qI[i][j][1])*(pI[i][j][1] - qI[i][j][1])) +
                                        ((pI[i][j][2] - qI[i][j][2])*(pI[i][j][2] - qI[i][j][2])));
                                }
                            }
                        }
                        float ssd3s = sumMatrix(p_s,p_s,ssd3);
                        if (min_sum > ssd3s)
                        {
                            min_sum = ssd3s;
                            qx      = x;
                            qy      = y;
                        }
                    }
                }
            }
            //printfFnc("qx %d qy %d \n",qx, qy);
            if (qx != -1 && qy != -1)
            {
                sub3Matrix3D(p_r,pnx,pny,nz,PI,qx+p_r,qy+p_r,qI);
                for(x=0; x < p_s; x++)
                {
                    for(y=0; y < p_s; y++)
                    {
                        if(maxX-p_r+x > -1
                        && maxX-p_r+x < nx
                        && maxY-p_r+y > -1
                        && maxY-p_r+y < ny
                        && mask[maxX+x-p_r][maxY+y-p_r] == 0)
                        {
                            I[maxX-p_r+x][maxY-p_r+y][0] = qI[x][y][0];
                            I[maxX-p_r+x][maxY-p_r+y][1] = qI[x][y][1];
                            I[maxX-p_r+x][maxY-p_r+y][2] = qI[x][y][2];
                            //printfFnc("Aktualizacja obrazu. \n");
                        }
                        if(maxX+x-p_r > -1
                        && maxX+x-p_r < nx
                        && maxY+y-p_r > -1
                        && maxY+y-p_r < ny)
                        {
                            mask[maxX+x-p_r][maxY+y-p_r] = 1;
                            //printfFnc("Aktualizacja maski. \n");
                        }
                    }
                }
            } else {
                printfFnc("Brak wzorca. \n");
                done = 1;
            }
        } else {
            printfFnc("Bez zmiany. \n");
            done = 1;
        }
        //printfFnc("allNonZero %d prev %d sumMatrix %f done %d", allNonZero(nx,ny,mask), prev, sumMatrix(nx,ny,mask), done);
    }

    pmxOut[0] = mxCreateDoubleMatrix(1,nx*ny*nz,mxREAL);
    double *ret1;
    ret1 = mxGetPr(pmxOut[0]);

    c_it=0;
    it  =0;
    while(c_it<nz)
    {
        while(it<(c_it+1)*nx*ny)
        {
            in2d=0;
            while(in2d<nx*ny)
            {
                ret1[it++] = I[in2d%nx][in2d/nx][c_it];
                in2d++;
            }
        }
        c_it++;
    }
    
    printfFnc("Początek czyszczenia zmiennych. \n");

    clearVars(p_r,p_s,nx,ny,nz,salAmount);
}
