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

#define dbl_eps 2.2204460492503131e-16 /* machine epsilon */
#define tol 2.2204460492503131e-14     /* 100*machine epsilon */
#define dbl_max 1.6e+308               /* highest double number */
#define maxit 100                      /* maximum number of iterations */

/* Includes */

#include "mex.h"
#include "math.h"

/* Computational Routine ---------------------------------------------------------------------------
 * Georg Schildbach, 06/May/2021 --- Gilbert-Johnson-Keerthi (GJK) Collision Detection in 2D
 * Calculates the value of the support function of a polytope X
 * -------------------------------------------------------------------------------------------------
 * suppfct(v,X,nX,k)
 * -------------------------------------------------------------------------------------------------
 * v: vector of search direction, to be maximized, vector of size 2
 * X: vector containing the vertices of X in R^2, vector of size 2*nX
 * nX: number of vertices in X, integer scalar
 * c: storage variable, vector of size 2
 * n: counter variable, integer scalar
 * k: number of a vertex realizing the value of the support function, integer scalar
 * -------------------------------------------------------------------------------------------------
 * Determines the vertex k that maximizes the support function argmax(v'*X(:,k))
 * -------------------------------------------------------------------------------------------------
 */

void suppfct(double *v, double *X, size_t *nX, double *c, size_t *n, size_t *k)
{

    c[0] = -dbl_max;

    for ( *n=0; *n<*nX; (*n)++) {
        c[1] = (v[0]*X[2*(*n)]) + (v[1]*X[2*(*n)+1]);
        if (c[1] > c[0]) {
            c[0] = c[1];
            *k = *n;
        }
    }
    
}

/* Gateway Function --------------------------------------------------------------------------------
 * Georg Schildbach, 05/May/2021 --- Gilbert-Johnson-Keerthi (GJK) Collision Detection in 2D
 * Matlab mex-function
 * -------------------------------------------------------------------------------------------------
 * [flag,dist] = gjkdist(X,Y)
 * -------------------------------------------------------------------------------------------------
 * INPUTS:
 * - V-polytope 1: vertices of X in R^2 (1 x 2*nX)
 * - V-polytope 2: vertices of Y in R^2 (1 x 2*nY)
 * -------------------------------------------------------------------------------------------------
 * OUTPUTS:
 * - flag: indicates an overlap between X and Y (true / false)
 * - dist: distance between X and Y (= 0 if there is a collision)
 * -------------------------------------------------------------------------------------------------
 */

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *X;             /* vertices of X in R^2 */
    double *Y;             /* vertices of Y in R^2 */
    size_t nX;             /* size of X */
    size_t nY;             /* size of Y */
    bool *flag;            /* flag */
    double *dist;          /* distance */
    
    size_t k;
    size_t n;
    size_t kX;           /* selected node of X */
    size_t kY;           /* selected node of Y */
    size_t nS = 0;       /* size of the simplex */
    double v[2];
    double w[2];
    double lambda[3];
    double D[3];
    double S[6];         /* simplex */

    /* 1) Create and check input and outputs -----------------------------------------------------*/
    
    /* 1.1) Check number of inputs and outputs */

    if(nrhs!=2) {
        mexErrMsgIdAndTxt("GJKCollisionCheck:gjkdist:nrhs","Two inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("GJKCollisionCheck:gjkdist:nlhs","Two outputs required.");
    }
    
    /* 1.2) Make sure the inputs are type double */
    
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("GJKCollisionCheck:gjkdist:notDouble","Input 1 must be type double.");
    }
    
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("GJKCollisionCheck:gjkdist:notDouble","Input 2 must be type double.");
    }
    
    /* 1.3) Check that number of rows in both inputs is 1 */
    
    if(mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("GJKCollisionCheck:gjkdist:notRowVector","Input 1 must be a row vector.");
    }
    
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("GJKCollisionCheck:gjkdist:notRowVector","Input 2 must be a row vector.");
    }
    
    /* 1.4) Get the number of columns in both inputs and check that they are even */
    
    k = mxGetN(prhs[0]);
    nX = (size_t) k/2;
    if(k!=2*nX) {
        mexErrMsgIdAndTxt("GJKCollisionCheck:gjkdist:notVertexVector","Input 1 must have an even number of elements.");
    }
    
    k = mxGetN(prhs[1]);
    nY = (size_t) k/2;
    if(k!=2*nY) {
        mexErrMsgIdAndTxt("GJKCollisionCheck:gjkdist:notVertexVector","Input 2 must have an even number of elements.");
    }
    
    /* 1.5) Get the vertices of X and Y in R^2 */

    X = mxGetPr(prhs[0]);
    Y = mxGetPr(prhs[1]);
    
    /* 1.6) Create the outputs */
    
    plhs[0] = mxCreateLogicalScalar(0);
    plhs[1] = mxCreateDoubleScalar(0);
    flag = mxGetLogicals(plhs[0]);
    dist = mxGetPr(plhs[1]);
    
    /* 2) GJK Algorithm --------------------------------------------------------------------------*/
    
    S[0] = X[0] - Y[0];  /* first node, x */
    S[1] = X[1] - Y[1];  /* first node, y */
    
    for (k=0; k<maxit; k++) {
        
    /* 2.1) Simplex has one node */
        
        if (nS == 0) {
            if ( (fabs(S[0])<tol) && (fabs(S[1])<=tol) ) {
                *flag = true;
                *dist = 0;
                break;  /* optimum found and equal to zero */
            }
            else {                
                suppfct(&S[0],Y,&nY,&v[0],&n,&kY);
                w[0] = -S[0];
                w[1] = -S[1];
                suppfct(&w[0],X,&nX,&v[0],&n,&kX);
                v[0] = X[kX*2  ] - Y[kY*2  ];
                v[1] = X[kX*2+1] - Y[kY*2+1];
                if ( (fabs(S[0]-v[0])<tol) && (fabs(S[1]-v[1])<=tol) ) {
                    *flag = false;
                    *dist = (v[0]*v[0]) + (v[1]*v[1]);
                    *dist = sqrt(*dist);
                    break;  /* optimum found and equal to zero */
                }
                else {
                    nS = 1;
                    S[2] = v[0];  /* add second node */
                    S[3] = v[1];
                }
            }
        }
        
    /* 2.2) Simplex consisting of two nodes */
 
        else if (nS == 1) {
            w[0] = S[2] - S[0];
            w[1] = S[3] - S[1];
            lambda[1] = (w[0]*w[0]) + (w[1]*w[1]);
            lambda[1] = sqrt(lambda[1]);
            if ( lambda[1]==0 ) {
                nS = 0;
            }
            else {
                w[0] = w[0] / lambda[1];
                w[1] = w[1] / lambda[1];
                lambda[0] = - (S[0]*w[0]) - (S[1]*w[1]);
                if ( lambda[0]<0 ) {
                    lambda[0] = 0;
                }
                else if ( lambda[0]>lambda[1] ) {
                    lambda[0] = lambda[1];
                }
                w[0] = S[0] + lambda[0]*w[0];  /* closest point on the node connection to the origin */
                w[1] = S[1] + lambda[0]*w[1];
                if ( (fabs(w[0])<=tol) && (fabs(w[1])<=tol) ) {
                    *flag = true;
                    *dist = 0;
                    break;  /* optimum found and equal to zero */
                }
                else {
                    suppfct(&w[0],Y,&nY,&v[0],&n,&kY);
                    w[0] = -w[0];
                    w[1] = -w[1];
                    suppfct(&w[0],X,&nX,&v[0],&n,&kX);
                    v[0] = X[kX*2  ] - Y[kY*2  ];
                    v[1] = X[kX*2+1] - Y[kY*2+1];
                    if ( (fabs(v[0]-S[0])<=tol) && (fabs(v[1]-S[1])<=tol) ) {
                        *flag = false;
                        *dist = (w[0]*w[0]) + (w[1]*w[1]);
                        *dist = sqrt(*dist);
                        break;  /* optimum found and not equal to zero */
                    }
                    else if ( (fabs(v[0]-S[2])<=tol) && (fabs(v[1]-S[3])<=tol) ) {
                        *flag = false;
                        *dist = (w[0]*w[0]) + (w[1]*w[1]);
                        *dist = sqrt(*dist);
                        break;  /* optimum found and not equal to zero */
                    }
                    else {
                        nS = 2;
                        S[4] = v[0]; /* add third node */
                        S[5] = v[1]; /* add third node */
                    }
                }
            }
        }
                                
                                
	/* 2.4) Simplex consisting of three nodes */

        else if ( nS==2 ) {
            D[0] = S[0]*S[3] - S[2]*S[1];
            D[1] = S[0]*S[5] - S[4]*S[1];
            D[2] = S[2]*S[5] - S[4]*S[3];
            if ( D[0]!=0 ) {
                lambda[0] = - D[2] / D[0];
                lambda[1] = D[1] / D[0];
                lambda[2] = -1;
            }
            else if ( D[1]!=0 ) {
                lambda[0] = D[2] / D[1];
                lambda[1] = -1;
                lambda[2] = D[0] / D[1];
            }
            else if ( D[2]!=0 ) {
                lambda[0] = -1;
                lambda[1] = D[1] / D[2];
                lambda[2] = - D[0] / D[2];
            }
            else {  /* points must be co-linear */
                if      ( (S[0]==S[4]) || (S[1]==S[5]) ) {
                    nS = 1;
                }
                else if ( (S[2]==S[4]) || (S[3]==S[5]) ) {
                    nS = 1;
                }
                else if ( (S[0]==S[2]) || (S[1]==S[3]) ) {
                    nS = 1;
                    S[0] = S[4];
                    S[1] = S[5];
                }
                else {
                    lambda[0] = (S[0]-S[2]) / (S[4]-S[2]);
                    lambda[1] = (S[2]-S[4]) / (S[0]-S[4]);
                    lambda[2] = (S[4]-S[0]) / (S[2]-S[0]);
                    if      ( (lambda[0]>=0) && (lambda[0]<=1) ) {
                        nS = 1;
                        S[0] = S[4];
                        S[1] = S[5];
                    }
                    else if ( (lambda[1]>=0) && (lambda[1]<=1) ) {
                        nS = 1;
                        S[2] = S[4];
                        S[3] = S[5];
                    }
                    else if ( (lambda[2]>=0) && (lambda[2]<=1) ) {
                        nS = 1;
                    }
                }
            }
            if ( nS==2 ) {  /* the regular case */
                D[0] = lambda[0] + lambda[1] + lambda[2];
                lambda[0] = lambda[0] / D[0];
                lambda[1] = lambda[1] / D[0];
                lambda[2] = lambda[2] / D[0];
                if ( (lambda[0]>-tol) && (lambda[1]>-tol) && (lambda[2]>-tol) ) {
                    *flag = true;
                    *dist = 0;
                    break;  /* optimum found and equal to zero */
                }
                else {
                    w[0] = S[2] - S[0];
                    w[1] = S[3] - S[1];
                    lambda[1] = (w[0]*w[0]) + (w[1]*w[1]);
                    lambda[1] = sqrt(lambda[1]);
                    w[0] = w[0] / lambda[1];
                    w[1] = w[1] / lambda[1];
                    lambda[0] = - (S[0]*w[0]) - (S[1]*w[1]);
                    if ( lambda[0]<0 ) {
                        lambda[0] = 0;
                    }
                    else if ( lambda[0]>lambda[1] ) {
                        lambda[0] = lambda[1];
                    }
                    w[0] = S[0] + (lambda[0]*w[0]);  /* closest point on the node connection to the origin */
                    w[1] = S[1] + (lambda[0]*w[1]);
                    lambda[2] = (w[0]*w[0]) + (w[1]*w[1]);
                    lambda[2] = sqrt(lambda[2]);
                    v[0] = S[4] - S[2];
                    v[1] = S[5] - S[3];
                    lambda[1] = (v[0]*v[0]) + (v[1]*v[1]);
                    lambda[1] = sqrt(lambda[1]);
                    v[0] = v[0] / lambda[1];
                    v[1] = v[1] / lambda[1];
                    lambda[0] = - (S[2]*v[0]) - (S[3]*v[1]);
                    if ( lambda[0]<0 ) {
                        lambda[0] = 0;
                    }
                    else if ( lambda[0]>lambda[1] ) {
                        lambda[0] = lambda[1];
                    }
                    v[0] = S[2] + (lambda[0]*v[0]);  /* closest point on the node connection to the origin */
                    v[1] = S[3] + (lambda[0]*v[1]);
                    lambda[1] = (v[0]*v[0]) + (v[1]*v[1]);
                    lambda[1] = sqrt(lambda[1]);
                    if ( lambda[1]<lambda[2] ) {
                        lambda[2] = lambda[1];
                        w[0] = v[0];
                        w[1] = v[1];
                        nS = 0;  /* farthest node of the simplex from the origin */
                    }
                    v[0] = S[0] - S[4];
                    v[1] = S[1] - S[5];
                    lambda[1] = (v[0]*v[0]) + (v[1]*v[1]);
                    lambda[1] = sqrt(lambda[1]);
                    v[0] = v[0] / lambda[1];
                    v[1] = v[1] / lambda[1];
                    lambda[0] = - (S[4]*v[0]) - (S[5]*v[1]);
                    if ( lambda[0]<0 ) {
                        lambda[0] = 0;
                    }
                    else if ( lambda[0]>lambda[1] ) {
                        lambda[0] = lambda[1];
                    }
                    v[0] = S[4] + (lambda[0]*v[0]);  /* closest point on the node connection to the origin */
                    v[1] = S[5] + (lambda[0]*v[1]);
                    lambda[1] = (v[0]*v[0]) + (v[1]*v[1]);
                    lambda[1] = sqrt(lambda[1]);
                    if ( lambda[1]<lambda[2] ) {
                        lambda[2] = lambda[1];
                        w[0] = v[0];
                        w[1] = v[1];
                        nS = 1;  /* farthest node of the simplex from the origin */
                    }
                    suppfct(&w[0],Y,&nY,&v[0],&n,&kY);
                    w[0] = -w[0];
                    w[1] = -w[1];
                    suppfct(&w[0],X,&nX,&v[0],&n,&kX);
                    v[0] = X[kX*2  ] - Y[kY*2  ];
                    v[1] = X[kX*2+1] - Y[kY*2+1];
                    if ( (fabs(v[0]-S[0])<=tol) && (fabs(v[1]-S[1])<=tol) ) {
                        *flag = false;
                        *dist = (w[0]*w[0]) + (w[1]*w[1]);
                        *dist = sqrt(*dist);
                        break;  /* optimum found and not equal to zero */
                    }
                    else if ( (fabs(v[0]-S[2])<=tol) && (fabs(v[1]-S[3])<=tol) ) {
                        *flag = false;
                        *dist = (w[0]*w[0]) + (w[1]*w[1]);
                        *dist = sqrt(*dist);
                        break;  /* optimum found and not equal to zero */
                    }
                    else if ( (fabs(v[0]-S[4])<=tol) && (fabs(v[1]-S[5])<=tol) ) {
                        *flag = false;
                        *dist = (w[0]*w[0]) + (w[1]*w[1]);
                        *dist = sqrt(*dist);
                        break;  /* optimum found and not equal to zero */
                    }
                    else {
                        S[nS*2  ] = v[0]; /* replace farthest node */
                        S[nS*2+1] = v[1]; /* replace farthest node */
                        nS = 2;
                    }
                }
            }
        }
    }

}


