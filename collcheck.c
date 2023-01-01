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

#define NumPi 3.141592653589793  /* numerical value of pi */

/* Includes */

#include "mex.h"
#include "math.h"

/* Gateway Function --------------------------------------------------------------------------------
 * Georg Schildbach, 10/May/2021 --- Continuous collision checker for arcs
 * Matlab mex-function
 * -------------------------------------------------------------------------------------------------
 * [flag] = collcheck(Obs,M,phi,r)
 * -------------------------------------------------------------------------------------------------
 * INPUTS:
 * - Obs: List of vertices of the polytope (vertices on convex hull only, no repeated vertices!),
 *        vector of size (nObs*2)
 * - M: center point (x,y) of the arc, vector of size 2
 * - phi: initial angle phi[0] and final angle phi[1] of the arc, vector of size 2
 * - r: inner radius r[0] and outer radius r[1] of the arc, vector of size 2
 * -------------------------------------------------------------------------------------------------
 * OUTPUTS:
 * - flag: indicates an overlap between the arc and the obstacle (true / false)
 * -------------------------------------------------------------------------------------------------
 */

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *Obs;           /* vertices of Obs in R^2 */
    double *M;             /* center point (x,y) of the arc */
    double *phi;           /* initial angle phi[0] and final angle phi[1] of the arc */
    double *r;             /* inner radius r[0] and outer radius r[1] of the arc */
    size_t nObs;           /* number of vertices of Obs */
    bool *flag;            /* flag */
    
    bool ints;               /* intersection indicator */
    bool c1;
    bool c2;
    size_t k;                /* counter variable */
    short int n;
    double IntsP1[4];        /* list of intersection points with line 1 */
    unsigned short nIntsP1;  /* number of intersection points with line 1 */
    double rIntsP1[2];       /* radius of intersection points with line 1 */
    double IntsP2[4];        /* list of intersection points with line 2 */
    unsigned short nIntsP2;  /* number of intersection points with line 2 */
    double rIntsP2[2];       /* radius of intersection points with line 2 */
    double alpha;
    double beta;
    double v1[2];
    double v2[2];
    double e[2];
    double d[2];
    double w[2];
    double psi;


    /* 1) Create and check input and outputs -----------------------------------------------------*/
    
    /* 1.1) Check number of inputs and outputs */

    if(nrhs!=4) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:nrhs","Four inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:nlhs","One output required.");
    }
    
    /* 1.2) Make sure the inputs are type double */
    
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notDouble","Input 1 must be type double.");
    }
    
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notDouble","Input 2 must be type double.");
    }
    
    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notDouble","Input 3 must be type double.");
    }
    
    if( !mxIsDouble(prhs[3]) || 
         mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notDouble","Input 4 must be type double.");
    }
    
    /* 1.3) Check that number of rows in both inputs is 1 */
    
    if(mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notRowVector","Input 1 must be a row vector.");
    }
    
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notRowVector","Input 2 must be a row vector.");
    }
    
    if(mxGetM(prhs[2])!=1) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notRowVector","Input 3 must be a row vector.");
    }
    
    if(mxGetM(prhs[3])!=1) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notRowVector","Input 4 must be a row vector.");
    }
    
    /* 1.4) Check the number of columns in the all inputs */
    
    k = mxGetN(prhs[0]);
    nObs = (size_t) k/2;
    if(k!=2*nObs) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notVertexVector","Input 1 must have an even number of elements.");
    }
    
    if(mxGetN(prhs[1])!=2) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notVertexVector","Input 2 (center point M) must consist of two elements.");
    }
    
    if(mxGetN(prhs[2])!=2) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notVertexVector","Input 3 (inner and outer radius r) must consist of two elements.");
    }
    
    if(mxGetN(prhs[3])!=2) {
        mexErrMsgIdAndTxt("ContCollisionCheck:collcheck:notVertexVector","Input 4 (start and final angle phi) must consist of two elements.");
    }
    
    /* 1.5) Get the input values */

    Obs = mxGetPr(prhs[0]);
    M = mxGetPr(prhs[1]);
    phi = mxGetPr(prhs[2]);
    r = mxGetPr(prhs[3]);
    
    /* 1.6) Create the output */
    
    plhs[0] = mxCreateLogicalScalar(0);
    flag = mxGetLogicals(plhs[0]);
    
    /* 2) Collision Checking ---------------------------------------------------------------------*/
 
    *flag = false;          /* initialize collision indicator variable */
    nIntsP1 = 0;
    nIntsP2 = 0;

    for ( k=0; k<nObs; k++ ) {   /* obstacle edge number k */
        v1[0] = Obs[2*k  ];      /* vertex 1 */
        v1[1] = Obs[2*k+1];      /* vertex 1 */
        if ( k+1<nObs ) {
            v2[0] = Obs[2*k+2];  /* vertex 2 */
            v2[1] = Obs[2*k+3];  /* vertex 2 */
        }
        else {
            v2[0] = Obs[0];      /* vertex 2 */
            v2[1] = Obs[1];      /* vertex 2 */
        }

    /* 2.1) Intersection points with cone lines */

        ints = false;  /* intersection indicator variable */
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
                ints = true;
                *flag = true;
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
                ints = true;
                *flag = true;
            }
            nIntsP2++;
        }

        if ( !ints ) {  /* no intersection with cone lines */
            alpha = (e[0]*e[0]) + (e[1]*e[1]);
            alpha = sqrt(alpha);
            e[0] = e[0] / alpha;
            e[1] = e[1] / alpha;
            beta = ((M[0]-v1[0])*e[0]) + ((M[1]-v1[1])*e[1]);
            d[0] = v1[0] - M[0] + beta*e[0];  /* projection of M onto the edge, relative to M */
            d[1] = v1[1] - M[1] + beta*e[1];
            w[0] = (d[0]*d[0]) + (d[1]*d[1]);
            w[0] = sqrt(w[0]);
            if ( beta>alpha ) {
                w[1] = alpha;
            }
            else if ( beta<0 ) {
                w[1] = 0;
            }
            else {
                w[1] = beta;
            }
            e[0] = v1[0] - M[0] + w[1]*e[0];  /* projection of M onto the edge, relative to M */
            e[1] = v1[1] - M[1] + w[1]*e[1];
            w[1] = (e[0]*e[0]) + (e[1]*e[1]);
            w[1] = sqrt(w[1]);

    /* 2.2) Projection point in [r_in,r_out] */

            if ( (w[1]>=r[0]) && (w[1]<=r[1]) ) {  /* calculate psi of projection point  */
                if      ( e[0]>0 ) {
                    psi = atan(e[1]/e[0]);
                }
                else if ( e[0]<0 ) {
                    psi = atan(e[1]/e[0]) + NumPi;
                }
                else {
                    if      ( e[1]>0 ) {
                        psi = + NumPi/2;
                    }
                    else if ( e[1]<0 ) {
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
                    *flag = true;
                }
            }

     /* 2.3) Projection point in [0,r_in] */

            if ( (w[0]<=r[0]) && !*flag ) {
                if      ( d[0]>0 ) {
                    psi = atan(d[1]/d[0]);
                }
                else if ( d[0]<0 ) {
                    psi = atan(d[1]/d[0]) + NumPi;
                }
                else {
                    if      ( d[1]>0 ) {
                        psi = + NumPi/2;
                    }
                    else if ( d[1]<0 ) {
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
                if ( psi<=phi[0]+(NumPi/2) ) {  /* set lower bound */
                    d[0] = cos(psi-phi[0]);
                }
                else if ( psi<=phi[1]+(NumPi/2) ) {
                    d[0] = 0;
                }
                else {
                    d[0] = +1;
                }
                if ( psi<=phi[1] ) {  /* set upper bound */
                    d[1] = +1;
                }
                else if ( psi<=phi[1]+(NumPi/2) ) {
                    d[1] = cos(phi[1]-psi);
                }
                else {
                    d[1] = -1;
                }
                if ( (w[0]/r[0]>=d[0]) && (w[0]/r[0]<=d[1]) ) {
                    if ( (beta>=0) && (e[0]>=r[0]) ) {
                        if ( (beta<=alpha) || (e[1]<=r[1]) ) {
                            ints = true;
                            *flag = true;
                        }
                    }
                }
                if ( psi <= phi[1]-(3/2*NumPi) ) {  /* additional upper bound */
                    if ( w[0]/r[0] <= cos(phi[1]-psi) ) {
                        if ( (beta>=0) && (e[0]>=r[0]) ) {
                            if ( (beta<=alpha) || (e[1]<=r[1]) ) {
                                ints = true;
                                *flag = true;
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
                            ints = true;
                            *flag = true;
                        }
                    }
                }
                if ( psi >= phi[0]+(3/2*NumPi) ) {  /* additional upper bound */
                    if ( w[0]/r[0] <= cos(psi-phi[0]) ) {
                        if ( (beta<=alpha) && (e[1]>=r[0]) ) {
                            if ( (beta>=0) || (e[0]<=r[1]) ) {
                                ints = true;
                                *flag = true;
                            }
                        }
                    }
                } 
            }
        }

        if ( *flag ) {
            break;  /* collision already detected */
        }
    }

    /* 2.4) Full inlusion by the obstacle */

    if ( !*flag ) {  /* check for full inclusion, and set flag = true  */
        if ( (nIntsP1>1) && (nIntsP2>1) ) {
            c1 = false;
            c2 = false;
                if ( nIntsP1 > 2 ) {
                    nIntsP1 = 2;
                }
            for ( k=0; k<nIntsP1; k++ ) {
                if      ( rIntsP1[k] <= r[0] ) {
                    c1 = true;
                }
                else if ( rIntsP1[k] >= r[1] ) {
                    c2 = true;
                }
            }
            if ( c1 && c2 ) {
                c1 = false;
                c2 = false;
                if ( nIntsP2 > 2 ) {
                    nIntsP2 = 2;
                }
                for ( k=0; k<nIntsP2; k++ ) {
                    if      ( rIntsP2[k] <= r[0] ) {
                        c1 = true;
                    }
                    else if ( rIntsP2[k] >= r[1] ) {
                        c2 = true;
                    }
                }
                if ( c1 && c2 ) {
                    *flag = true;
                }
            }
        }
    }

}


