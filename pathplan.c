/* Copyright (C) 2023, Georg Schildbach.
 * -------------------------------------------------------------------------------------------------
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 * NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * -------------------------------------------------------------------------------------------------
 */

/* Macros */

#define maxit 2000               /* maximum number of RRT iterations */
#define maxnod 2000              /* maximum number of constructed nodes */
#define singtol (1e-10)          /* singularity tolerance */
#define NumPi 3.141592653589793  /* numerical value of pi */

/* Includes */

#include "mex.h"
#include "math.h"
#include "stdlib.h"


/* Computational Routine ---------------------------------------------------------------------------
 * Georg Schildbach, 10/May/2021 --- CCD + RRT path planning method
 * Continuous collision checker for arcs
 * -------------------------------------------------------------------------------------------------
 * collcheck(Obs,M,phi,r,flag)
 * -------------------------------------------------------------------------------------------------
 * INPUTS:
 * - nObs: number of vertices of the polytope
 * - Obs: List of vertices of the polytope (vertices on convex hull only, no repeated vertices!),
 *        vector of size (nObs*2)
 * - M: center point (x,y) of the arc, vector of size 2
 * - phi: initial angle phi[0] and final angle phi[1] of the arc, vector of size 2
 * - r: inner radius r[0] and outer radius r[1] of the arc, vector of size 2
 * -------------------------------------------------------------------------------------------------
 * OUTPUTS:
 * - flag: indicates an overlap between the arc and the obstacle (1 = true / 0 = false)
 * -------------------------------------------------------------------------------------------------
 */

void collcheck(const double *nObs, const double *Obs, const double *M, const double *phi, const double *r, bool *flag)

{
    
    static bool ints;               /* intersection indicator */
    static bool c1;
    static bool c2;
    static size_t k;                /* counter variable */
    static short int n;
    static double IntsP1[4];        /* list of intersection points with line 1 */
    static unsigned short nIntsP1;  /* number of intersection points with line 1 */
    static double rIntsP1[2];       /* radius of intersection points with line 1 */
    static double IntsP2[4];        /* list of intersection points with line 2 */
    static unsigned short nIntsP2;  /* number of intersection points with line 2 */
    static double rIntsP2[2];       /* radius of intersection points with line 2 */
    static double alpha;
    static double beta;
    static double v1[2];
    static double v2[2];
    static double e[2];
    static double d[2];
    static double w[2];
    static double psi;

    /* 1) Initialize -----------------------------------------------------------------------------*/
 
    *flag = 0;         /* initialize collision indicator variable */
    nIntsP1 = 0;
    nIntsP2 = 0;

    for ( k=0; k<((int)*nObs); k++ ) {   /* obstacle edge number k */
        v1[0] = Obs[2*k  ];      /* vertex 1 */
        v1[1] = Obs[2*k+1];      /* vertex 1 */
        if ( k+1<((int)*nObs) ) {
            v2[0] = Obs[2*k+2];  /* vertex 2 */
            v2[1] = Obs[2*k+3];  /* vertex 2 */
        }
        else {
            v2[0] = Obs[0];      /* vertex 2 */
            v2[1] = Obs[1];      /* vertex 2 */
        }

    /* 2) Intersection points with cone lines ----------------------------------------------------*/

        ints = 0;  /* intersection indicator variable */
        e[0] = v2[0] - v1[0];
        e[1] = v2[1] - v1[1];
        d[0] = M[0] - v1[0];
        d[1] = M[1] - v1[1];
        if ( d[1]*e[0]<d[0]*e[1] ) {  /* switch vertices s.t. v1 has a smaller angle */
            w[0] = v1[0];
            w[1] = v1[1];
            v1[0] = v2[0];
            v1[1] = v2[1];
            v2[0] = w[0];
            v2[1] = w[1];
            e[0] = v2[0] - v1[0];
            e[1] = v2[1] - v1[1];
            d[0] = M[0] - v1[0];
            d[1] = M[1] - v1[1];
        }

        /* Cone line 1 */
        w[0] = cos(phi[0]);
        w[1] = sin(phi[0]);
        if ( e[0]!=0 ) {
            beta = (d[1] - (e[1]/e[0]*d[0])) / ((e[1]/e[0]*w[0]) - w[1]);
            alpha = (d[0]/e[0]) + (w[0]/e[0]*beta);
        }
        else if ( e[1]!=0 ) {
            beta = (d[0] - (e[0]/e[1]*d[1])) / ((e[0]/e[1]*w[1]) - w[0]);
            alpha = (d[1]/e[1]) + (w[1]/e[1]*beta);
        }
        if ( (alpha>=0) && (alpha<=1) ) {  /* intersection is within edge */
            IntsP1[2*nIntsP1  ] = v1[0] + (alpha*e[0]) - M[0];  /* intersection point relative to M */
            IntsP1[2*nIntsP1+1] = v1[1] + (alpha*e[1]) - M[1];
            alpha = (IntsP1[2*nIntsP1]*IntsP1[2*nIntsP1]) + (IntsP1[2*nIntsP1+1]*IntsP1[2*nIntsP1+1]);
            alpha = sqrt(alpha);
            if ( beta>=0 ) {
                rIntsP1[nIntsP1] = + alpha;  /* radius of intersection point */
            }
            else {
                rIntsP1[nIntsP1] = - alpha;  /* radius of intersection point */
            }
            if ( (rIntsP1[nIntsP1]>=r[0]) && (rIntsP1[nIntsP1]<=r[1]) ) {
                ints = 1;
                *flag = 1;
            }
            nIntsP1++;
        }

        /* Cone line 2 */
        w[0] = cos(phi[1]);
        w[1] = sin(phi[1]);
        if ( e[0]!=0 ) {
            beta = (d[1] - (e[1]/e[0]*d[0])) / ((e[1]/e[0]*w[0]) - w[1]);
            alpha = (d[0]/e[0]) + (w[0]/e[0]*beta);
        }
        else if ( e[1]!=0 ) {
            beta = (d[0] - (e[0]/e[1]*d[1])) / ((e[0]/e[1]*w[1]) - w[0]);
            alpha = (d[1]/e[1]) + (w[1]/e[1]*beta);
        }
        if ( (alpha>=0) && (alpha<=1) ) {  /* intersection is within edge */
            IntsP2[2*nIntsP2  ] = v1[0] + (alpha*e[0]) - M[0];  /* intersection point relative to M */
            IntsP2[2*nIntsP2+1] = v1[1] + (alpha*e[1]) - M[1];
            alpha = (IntsP2[2*nIntsP2]*IntsP2[2*nIntsP2]) + (IntsP2[2*nIntsP2+1]*IntsP2[2*nIntsP2+1]);
            alpha = sqrt(alpha);
            if ( beta>=0 ) {
                rIntsP2[nIntsP2] = + alpha; /* radius of intersection point */
            }
            else {
                rIntsP2[nIntsP2] = - alpha; /* radius of intersection point */
            }
            if ( (rIntsP2[nIntsP2]>=r[0]) && (rIntsP2[nIntsP2]<=r[1]) ) {
                ints = 1;
                *flag = 1;
            }
            nIntsP2++;
        }

        if ( ints == 0 ) {  /* no intersection with cone lines */
            alpha = (e[0]*e[0]) + (e[1]*e[1]);
            alpha = sqrt(alpha);
            e[0] = e[0] / alpha;
            e[1] = e[1] / alpha;
            beta = ((M[0]-v1[0])*e[0]) + ((M[1]-v1[1])*e[1]);
            d[0] = v1[0] - M[0] + beta*e[0];  /* projection of M onto the edge, relative to M */
            d[1] = v1[1] - M[1] + beta*e[1];
            w[0] = (d[0]*d[0]) + (d[1]*d[1]);
            w[0] = sqrt(w[0]);
            if ( beta > alpha ) {
                w[1] = alpha;
            }
            else if ( beta < 0 ) {
                w[1] = 0;
            }
            else {
                w[1] = beta;
            }
            e[0] = v1[0] - M[0] + w[1]*e[0];  /* projection of M onto the edge, relative to M */
            e[1] = v1[1] - M[1] + w[1]*e[1];
            w[1] = (e[0]*e[0]) + (e[1]*e[1]);
            w[1] = sqrt(w[1]);

    /* 3) Projection point in [r_in,r_out] -------------------------------------------------------*/

            if ( (w[1]>=r[0]) && (w[1]<=r[1]) ) {  /* calculate psi of projection point  */
                if      ( e[0] > 0 ) {
                    psi = atan(e[1]/e[0]);
                }
                else if ( e[0] < 0 ) {
                    psi = atan(e[1]/e[0]) + NumPi;
                }
                else {
                    if      ( e[1] > 0 ) {
                        psi = + NumPi/2;
                    }
                    else if ( e[1] < 0 ) {
                        psi = - NumPi/2;
                    }
                    else {
                        psi = 0;
                    }
                }
                n = (short int) ((psi-phi[0])/(2*NumPi));
                if ( psi < phi[0] ) {
                    n--;
                }
                psi = psi - (n*2*NumPi);  /* psi in [phi[0],phi[0]+2*pi) */
                if ( (psi>=phi[0]) && (psi<=phi[1]) ) {  /* projection point is in the cone */
                    *flag = 1;
                }
            }

     /* 4) Projection point in [0,r_in] ----------------------------------------------------------*/

            if ( (w[0]<=r[0]) && (*flag==0) ) {
                if      ( d[0] > 0 ) {
                    psi = atan(d[1]/d[0]);
                }
                else if ( d[0] < 0 ) {
                    psi = atan(d[1]/d[0]) + NumPi;
                }
                else {
                    if      ( d[1] > 0 ) {
                        psi = + NumPi/2;
                    }
                    else if ( d[1] < 0 ) {
                        psi = - NumPi/2;
                    }
                    else {
                        psi = 0;
                    }
                }
                d[0] = v1[0] - M[0];
                d[1] = v1[1] - M[1];
                e[0] = (d[0]*d[0]) + (d[1]*d[1]);
                e[0] = sqrt(e[0]);  /* radius of vertex 1 */
                d[0] = v2[0] - M[0];
                d[1] = v2[1] - M[1];
                e[1] = (d[0]*d[0]) + (d[1]*d[1]);
                e[1] = sqrt(e[1]);  /* radius of vertex 2 */

                /* Vertex 1: check (e[0] < r[0]) iff w[0]/r[0] in [d[0],d[1]] */  
                n = (short int) ((psi-phi[0])/(2*NumPi));
                if ( psi < phi[0] ) {
                    n--;
                }
                psi = psi - (n*2*NumPi);  /* psi in [phi[0],phi[0]+2*pi) */
                if ( psi <= phi[0]+(NumPi/2) ) {  /* set lower bound */
                    d[0] = cos(psi-phi[0]);
                }
                else if ( psi <= phi[1]+(NumPi/2) ) {
                    d[0] = 0;
                }
                else {
                    d[0] = +1;
                }
                if ( psi <= phi[1] ) {  /* set upper bound */
                    d[1] = +1;
                }
                else if ( psi <= phi[1]+(NumPi/2) ) {
                    d[1] = cos(phi[1]-psi);
                }
                else {
                    d[1] = -1;
                }
                if ( (w[0]/r[0]>=d[0]) && (w[0]/r[0]<=d[1]) ) {
                    if ( (beta>=0) && (e[0]>=r[0]) ) {
                        if ( (beta<=alpha) || (e[1]<=r[1]) ) {
                            ints = 1;
                            *flag = 1;
                        }
                    }
                }
                if ( psi <= phi[1]-(3/2*NumPi) ) {  /* additional upper bound */
                    if ( w[0]/r[0] <= cos(phi[1]-psi) ) {
                        if ( (beta>=0) && (e[0]>=r[0]) ) {
                            if ( (beta<=alpha) || (e[1]<=r[1]) ) {
                                ints = 1;
                                *flag = 1;
                            }
                        }
                    }
                }

                /* Vertex 2: check (e[1] < r[0]) iff w[0]/r[0] in [d[0],d[1]] */
                n = (short int) ((psi-phi[1])/(2*NumPi));
                if ( psi >= phi[1] ) {
                    n++;
                }
                psi = psi - (n*2*NumPi);  /* psi in (phi[1]-2*pi,phi[1]] */
                if ( psi >= phi[0] ) {   /* set upper bound */
                    d[1] = +1;
                }
                else if ( psi >= phi[0]-(NumPi/2) ) {
                    d[1] = cos(psi-phi[0]);
                }
                else {
                    d[1] = -1;
                }
                if ( psi >= phi[1]-(NumPi/2) ) {  /* set lower bound */
                    d[0] = cos(psi-phi[1]);
                }
                else if ( psi >= phi[0]-(NumPi/2) ) {
                    d[0] = 0;
                }
                else {
                    d[0] = +1;
                }
                if ( (w[0]/r[0]>=d[0]) && (w[0]/r[0]<=d[1]) ) {
                    if ( (beta<=alpha) && (e[1]>=r[0]) ) {
                        if ( (beta>=0) || (e[0]<=r[1]) ) {
                            ints = 1;
                            *flag = 1;
                        }
                    }
                }
                if ( psi >= phi[0]+(3/2*NumPi) ) {  /* additional upper bound */
                    if ( w[0]/r[0] <= cos(psi-phi[0]) ) {
                        if ( (beta<=alpha) && (e[1]>=r[0]) ) {
                            if ( (beta>=0) || (e[0]<=r[1]) ) {
                                ints = 1;
                                *flag = 1;
                            }
                        }
                    }
                } 
            }
        }

        if ( *flag == 1 ) {
            break;  /* collision already detected */
        }
    }

    /* 5) Full inlusion by the obstacle ----------------------------------------------------------*/

    if ( *flag == 0 ) {  /* check for full inclusion, and set flag = 1  */
        if ( (nIntsP1>1) && (nIntsP2>1) ) {
            c1 = 0;
            c2 = 0;
                if ( nIntsP1 > 2 ) {
                    nIntsP1 = 2;
                }
            for ( k=0; k<nIntsP1; k++ ) {
                if      ( rIntsP1[k] <= r[0] ) {
                    c1 = 1;
                }
                else if ( rIntsP1[k] >= r[1] ) {
                    c2 = 1;
                }
            }
            if ( (c1==1) && (c2==1) ) {
                c1 = 0;
                c2 = 0;
                if ( nIntsP2 > 2 ) {
                    nIntsP2 = 2;
                }
                for ( k=0; k<nIntsP2; k++ ) {
                    if      ( rIntsP2[k] <= r[0] ) {
                        c1 = 1;
                    }
                    else if ( rIntsP2[k] >= r[1] ) {
                        c2 = 1;
                    }
                }
                if ( (c1==1) && (c2==1) ) {
                    *flag = 1;
                }
            }
        }
    }

}

/* Gateway Function --------------------------------------------------------------------------------
 * Georg Schildbach, 18/Dec/2022 --- CCD + RRT path planning method
 * Matlab mex-function
 * -------------------------------------------------------------------------------------------------
 * [flag] = collcheck(Obs,M,phi,r)
 * -------------------------------------------------------------------------------------------------
 * INPUTS:
 * - car: parameters of the car (w,lf,lr,lff,lrr,delta_max), vector of size 6
 * - xlim: map size (x_min, x_max, y_min, y_max), vector of size 4
 * - sObs: number of vertex points for each obstacle, vector of size nObs
 * - Obst: array with obstacle vertices in counter-clockwise order, vector of size nObs*2*max(sObs)
 * - Init: initial node (x, y, psi), vector of size 3
 * - Targ: target node (x, y, psi), vector of size 3 
 * - DeltaTarg: target tolerance (Delta x, Delta y, Delta psi), vector of size 3
 * -------------------------------------------------------------------------------------------------
 * OUTPUTS:
 * - Nodes: array containint the optimal path information, vector of size 3*optnNod
 * -------------------------------------------------------------------------------------------------
 */

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *car;          /* parameters of the car */
    double *xlim;         /* map size in R^2 */
    double *sObs;         /* number of vertices for each obstacle */
    double *Obst;         /* vertices of obstacles in R^2 */
    double *Init;         /* initial node (x, y, psi) */
    double *Targ;         /* target node (x, y, psi) */
    double *DeltaTarg;    /* target tolerance (Delta x, Delta y, Delta psi) */
    double *Nodes;        /* Nodes of the planned path */

    size_t nObs;          /* number of obstacles */
    int nNod = 2;         /* number of nodes, including initial and target node */
    int aNod = 1;         /* most advanced node index, based on distance to target (0,...,nNod-1) */

    /* auxiliary variables */
    bool flag;
    bool coll;
    bool posturn;
    bool newnode;
    bool targetnode;
    int i, j, k, l;
    int N;
    double R, R_min;
    double ang, dist;
    double P[2];
    double v[2];
    double q[2];
    double M[2];
    double r[2];
    double phi[2];

    int strat[5];
    strat[0] = 0;
    strat[1] = 0;
    strat[2] = 0;
    strat[3] = 0;
    strat[4] = 0;

    /* 1) Create and check input and outputs -----------------------------------------------------*/
    
    /* 1.1) Check number of inputs and outputs */

    if ( nrhs!=7 ) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:nrhs","Seven inputs required.");
    }
    if ( nlhs!=1 ) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:nlhs","One output required.");
    }
    
    /* 1.2) Make sure the inputs are of correct type */
    
    if ( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notDouble","Input 1 must be type double.");
    }

    if ( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notDouble","Input 2 must be type double.");
    }
    
    if ( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notDouble","Input 3 must be type double.");
    }
    
    if ( !mxIsDouble(prhs[3]) || 
         mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notDouble","Input 4 must be type double.");
    }

    if ( !mxIsDouble(prhs[4]) || 
         mxIsComplex(prhs[4])) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notDouble","Input 5 must be type double.");
    }
    
    if ( !mxIsDouble(prhs[5]) || 
         mxIsComplex(prhs[5])) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notDouble","Input 6 must be type double.");
    }
    
    if ( !mxIsDouble(prhs[6]) || 
         mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notDouble","Input 7 must be type double.");
    }
    
    /* 1.3) Check that number of rows in both inputs is 1 */
    
    if (mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notRowVector","Input 1 must be a row vector.");
    }
    
    if (mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notRowVector","Input 2 must be a row vector.");
    }
    
    if (mxGetM(prhs[2])!=1) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notRowVector","Input 3 must be a row vector.");
    }
    
    if (mxGetM(prhs[3])!=1) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notRowVector","Input 4 must be a row vector.");
    }

    if (mxGetM(prhs[4])!=1) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notRowVector","Input 5 must be a row vector.");
    }
    
    if (mxGetM(prhs[5])!=1) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notRowVector","Input 6 must be a row vector.");
    }
    
    if (mxGetM(prhs[6])!=1) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notRowVector","Input 7 must be a row vector.");
    }
    
    /* 1.4) Check the number of columns in the all inputs */
    
    if (mxGetN(prhs[0])!=6) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notVehicleParameters","Input 1 (vehicle parameters) must consist of 6 elements.");
    }

    if (mxGetN(prhs[1])!=4) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notMapLimits","Input 2 (map limits) must consist of 4 elements.");
    }

    nObs = mxGetN(prhs[2]);

    k = mxGetN(prhs[3]);  
    i = (int) k/2;
    if ( k!=2*i ) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notVertexVector","Input 4 (obstacle vertices) must have an even number of elements.");
    }
    
    if ( mxGetN(prhs[4])!=3 ) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notConfigurationVector","Input 5 (inital point) must consist of three elements.");
    }

    if ( mxGetN(prhs[5])!=3 ) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notConfigurationVector","Input 5 (target point) must consist of three elements.");
    }

    if ( mxGetN(prhs[6])!=3 ) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notConfigurationVector","Input 6 (target tolerance) must consist of three elements.");
    }
    
    /* 1.5) Get the input values */

    car = mxGetPr(prhs[0]);
    xlim = mxGetPr(prhs[1]);
    sObs = mxGetPr(prhs[2]);
    Obst = mxGetPr(prhs[3]);
    Init = mxGetPr(prhs[4]);
    Targ = mxGetPr(prhs[5]);
    DeltaTarg = mxGetPr(prhs[6]);
    
    /* 1.6) Check Obst vector for consistent size */
    
    j = 0;
    for ( i=0; i<nObs; i++ ) {
        if ( j<sObs[i] ) {
            j = sObs[i];
        }
    }
    if ( k!=2*nObs*j ) {
        mexErrMsgIdAndTxt("RRTPathPlanner:pathplan:notVertexVector","Input 4 (obstacle vertices) must be of consistent dimensions with sObs.");
    }
    
    /* 1.7) Create the output */
    
    plhs[0] = mxCreateDoubleMatrix(6,maxnod,mxREAL);
    Nodes = mxGetPr(plhs[0]);

    /* 1.8) Initialize */

    /* target node */
    Nodes[ 0] = Targ[0];
    Nodes[ 1] = Targ[1];
    Nodes[ 2] = Targ[2];
    Nodes[ 3] = pow(2,50);
    Nodes[ 4] = -1;
    Nodes[ 5] = 0;

    /* initial node */
    Nodes[ 6] = Init[0];
    Nodes[ 7] = Init[1];
    Nodes[ 8] = Init[2];
    Nodes[ 9] = 0;
    Nodes[10] = -1;
    Nodes[11] = sqrt((Targ[0]-Init[0])*(Targ[0]-Init[0])+(Targ[1]-Init[1])*(Targ[1]-Init[1]));

    /* minimum turning radius */
    R_min = (car[1]+car[2])/tan(NumPi/180*car[5])*cos(atan(car[2]/(car[1]+car[2])*tan(NumPi/180*car[5])));
    
    /* 2) Path Planning --------------------------------------------------------------------------*/
 
    for ( k=1; k<maxit; k++ ) {

        /* 2.1) Generate new random node */

        r[0] = (double)rand() / (double)RAND_MAX * 100;
        targetnode = 0;

        if ( Nodes[4] == -1 ) {  /* no connection to target point has been found */
            if ( (r[0]<5) || (k==maxit) ) {     /* try target point */
                P[0] = Targ[0];
                P[1] = Targ[1];
                targetnode = 1;
                strat[0] = strat[0] + 1;
            }
            else if ( r[0] < 20 ) {             /* progress from best known node towards target */
                r[1] = (double)rand() / (double)RAND_MAX;
                P[0] = Nodes[aNod*6  ] + r[1] * (Targ[0]-Nodes[aNod*6  ]);
                P[1] = Nodes[aNod*6+1] + r[1] * (Targ[1]-Nodes[aNod*6+1]);
                strat[1] = strat[1] + 1;
            }
            else if ( r[0] < 40 ) {             /*  progress around most advanced node */
                r[1] = (double)rand() / (double)RAND_MAX;
                P[0] = Nodes[aNod*6  ] + 0.2 * (xlim[0] + (xlim[1]-xlim[0]) * r[1]);
                r[1] = (double)rand() / (double)RAND_MAX;
                P[1] = Nodes[aNod*6+1] + 0.2 * (xlim[2] + (xlim[3]-xlim[2]) * r[1]);
                strat[2] = strat[2] + 1;
            }
            else if ( r[0] < 60 ) {             /* explore around target point */
                r[1] = (double)rand() / (double)RAND_MAX;
                P[0] = Targ[0] + 0.2 * (xlim[0] + (xlim[1]-xlim[0]) * r[1]);
                r[1] = (double)rand() / (double)RAND_MAX;
                P[1] = Targ[1] + 0.2 * (xlim[2] + (xlim[3]-xlim[2]) * r[1]);
                strat[3] = strat[3] + 1;
            }
            else {                              /* explore the whole map */
                r[1] = (double)rand() / (double)RAND_MAX;
                P[0] = xlim[0] + (xlim[1]-xlim[0]) * r[1];
                r[1] = (double)rand() / (double)RAND_MAX;
                P[1] = xlim[2] + (xlim[3]-xlim[2]) * r[1];
                strat[4] = strat[4] + 1;
            }
        }
        else {                   /* a connection to target point has already been found */
            if ( (r[0]<10) || (k==maxit) ) {    /* try target point */
                P[0] = Targ[0];
                P[1] = Targ[1];
                targetnode = 1;
                strat[0] = strat[0] + 1;
            }
            else if ( r[0] < 50 ) {              /* explore around target point */
                r[1] = (double)rand() / (double)RAND_MAX;
                P[0] = Targ[0] + 0.2 * (xlim[0] + (xlim[1]-xlim[0]) * r[1]);
                r[1] = (double)rand() / (double)RAND_MAX;
                P[1] = Targ[1] + 0.2 * (xlim[2] + (xlim[3]-xlim[2]) * r[1]);
                strat[3] = strat[3] + 1;
            }
            else {                               /* explore the whole map */
                r[1] = (double)rand() / (double)RAND_MAX;
                P[0] = xlim[0] + (xlim[1]-xlim[0]) * r[1];
                r[1] = (double)rand() / (double)RAND_MAX;
                P[1] = xlim[2] + (xlim[3]-xlim[2]) * r[1];
                strat[4] = strat[4] + 1;
            }
        }

        N = nNod;
        newnode = 1;
    
        for ( i=1; i<N; i++ ) {  /* for every stored node, except the target node */
        
            /* 2.2) Calculate connecting arc */
    
            /* Arc center point M */
            
            v[0] = cos(Nodes[i*6+2]); 
            v[1] = sin(Nodes[i*6+2]);
            q[0] = P[0] - Nodes[i*6  ];
            q[1] = P[1] - Nodes[i*6+1];
            r[0] = (v[0]*q[1])-(v[1]*q[0]);  /* determinant */
            if ( (r[0] > singtol) && (r[0] < +singtol ) ) {   /* matrix is singular */
                M[0] = (P[0]+Nodes[i*6  ])/2 - q[1]/singtol;
                M[1] = (P[1]+Nodes[i*6+1])/2 + q[0]/singtol;
            }
            else {                                            /* matrix is non-singular */
                M[0] = (  q[1]*(v[0]*Nodes[i*6  ]+v[1]*Nodes[i*6+1]) - 0.5*v[1]*(q[0]*(P[0]+Nodes[i*6  ])+q[1]*(P[1]+Nodes[i*6+1]))) / r[0];
                M[1] = (- q[0]*(v[0]*Nodes[i*6  ]+v[1]*Nodes[i*6+1]) + 0.5*v[0]*(q[0]*(P[0]+Nodes[i*6  ])+q[1]*(P[1]+Nodes[i*6+1]))) / r[0];
            }

            /* Turning direction */
    
            if ( v[0] * (M[1]-Nodes[i*6+1]) + v[1] * (Nodes[i*6  ]-M[0]) >= 0  ) {
                posturn = 1;  /* vehicle turning counter-clockwise (positive direction) */
            }
            else {
                posturn = 0;  /* vehicle turning clockwise (negative direction) */
            }
            
            /* Arc radius R */
            
            R = sqrt( (P[0]-M[0])*(P[0]-M[0]) + (P[1]-M[1])*(P[1]-M[1]) );
            r[0] = R - car[0]/2;  /* r_in */
            if ( car[3] >= car[4] ) {
                r[1] = sqrt( (R+car[0]/2)*(R+car[0]/2) + car[3]*car[3] );  /* r_out */
            }
            else {
                r[1] = sqrt( (R+car[0]/2)*(R+car[0]/2) + car[4]*car[4] );  /* r_out */
            }
            
            /* Arc limiting angles phi[0], phi[1] */
            
            if ( R >= R_min ) {
                if      ( Nodes[i*6  ]-M[0] > 0 ) {
                    phi[0] = atan((Nodes[i*6+1]-M[1])/(Nodes[i*6  ]-M[0]));
                }
                else if ( Nodes[i*6  ]-M[0] < 0 ) {
                    phi[0] = atan((Nodes[i*6+1]-M[1])/(Nodes[i*6  ]-M[0])) + NumPi;
                }
                else {
                    if      ( Nodes[i*6+1]-M[1] > 0 ) {
                        phi[0] = + NumPi/2;
                    }
                    else if ( Nodes[i*6+1]-M[1] < 0 ) {
                        phi[0] = - NumPi/2;
                    }
                    else {
                        phi[0] = 0;
                    }
                }
                if      ( P[0]-M[0] > 0 ) {
                    phi[1] = atan((P[1]-M[1])/(P[0]-M[0]));
                }
                else if ( P[0]-M[0] < 0 ) {
                    phi[1] = atan((P[1]-M[1])/(P[0]-M[0])) + NumPi;
                }
                else {
                    if      ( P[1]-M[1] > 0 ) {
                        phi[1] = + NumPi/2;
                    }
                    else if ( P[1]-M[1] < 0 ) {
                        phi[1] = - NumPi/2;
                    }
                    else {
                        phi[1] = 0;
                    }
                }
                if ( posturn == 1 ) {   /* vehicle turning counter-clockwise */
                    j = (int)((phi[1]-phi[0])/2/NumPi);
                    if ( phi[1]-phi[0] == j*2*NumPi ) {  /* phi[1] in [phi[0],phi[0]+2*pi) */
                        phi[1] = phi[0];
                    }
                    else if ( phi[1]-phi[0] > 0 ) {
                        phi[1] = phi[1] - j*2*NumPi;     
                    }
                    else {
                        phi[1] = phi[1] - (j-1)*2*NumPi;
                    }
                }
                else {                  /* vehicle turning clockwise */
                    j = (int)((phi[1]-phi[0])/2/NumPi);
                    if ( phi[1]-phi[0] == j*2*NumPi ) {  /* phi[1] in (phi[0]-2*pi,phi[0]] */
                        phi[1] = phi[0];
                    }
                    else if ( phi[1]-phi[0] > 0 ) {
                        phi[1] = phi[1] - (j+1)*2*NumPi;     
                    }
                    else {
                        phi[1] = phi[1] - j*2*NumPi;
                    }
                }
                ang = Nodes[i*6+2] + (phi[1]-phi[0]);
                if ( phi[1] >= phi[0] ) {
                    dist = Nodes[i*6+3] + ( R * (phi[1]-phi[0]) );
                }
                else {
                    dist = Nodes[i*6+3] + ( R * (phi[0]-phi[1]) );
                }
            }
    
            /* Check for admissible arc */
            
            coll = 1;
            if ( R >= R_min ) {
                coll = 0;
                if ( posturn == 1 ) {   /* vehicle turning counter-clockwise */
                    phi[0] = phi[0] - asin(car[4]/(R-car[0]/2));
                    phi[1] = phi[1] + asin(car[3]/(R-car[0]/2));
                }
                else {                  /* vehicle turning clockwise */
                    q[0] = phi[0];
                    phi[0] = phi[1] - asin(car[3]/(R-car[0]/2));
                    phi[1] = q[0]   + asin(car[4]/(R-car[0]/2));
                }
                l = 0;
                for ( j=0; j<nObs; j++ ) {
                    collcheck(&sObs[j],&Obst[2*l],&M[0],&phi[0],&r[0],&flag);
                    l = l + sObs[j];
                    if ( flag == 1 ) {
                        coll = 1;
                        break;
                    }
                }
            }
            
            /* 2.3) Add node */
    
            /* Case 1: Add new node or update node */
    
            if ( coll==0 && targetnode==0 ) {
                if ( newnode==1 || dist<Nodes[(nNod-1)*6+3] ) {
                    if ( newnode == 1 ) {
                        nNod++;
                        newnode = 0;
                    }
                    Nodes[(nNod-1)*6  ] = P[0];
                    Nodes[(nNod-1)*6+1] = P[1];
                    Nodes[(nNod-1)*6+2] = ang;
                    Nodes[(nNod-1)*6+3] = dist;
                    Nodes[(nNod-1)*6+4] = i + 1;
                    Nodes[(nNod-1)*6+5] = sqrt( (Targ[0]-P[0])*(Targ[0]-P[0]) + (Targ[1]-P[1])*(Targ[1]-P[1]) );
                    if ( Nodes[(nNod-1)*6+5] < Nodes[aNod*6+5] ) {
                        aNod = nNod - 1;
                    }
                }
            }
    
            /* Case 2: P is target node */
    
            if ( coll==0 && targetnode==1 ) {
                if ( ang>Targ[2]-DeltaTarg[2] && ang<Targ[2]+DeltaTarg[2] ) {
                    if ( dist < Nodes[3] ) {  /* new optimal path has been found! */
                        Nodes[3] = dist;
                        Nodes[4] = i + 1;
                    }
                }
            }
        }
        
        /* 2.4) Check for termination */
        
        if ( nNod >= maxnod ) {  /* maximum number of nodes reached */
            break;
        }
    }

}


