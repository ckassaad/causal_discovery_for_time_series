#include <string.h>
#include "mex.h"

/* WARNING: THIS FUNCTION PERFORMS NO ERROR CHECKING WHATSOEVER! Don't even
 * think of calling genvar_mex itself! It is designed to be called via the
 * genvar.m Matlab wrapper function (utils subfolder); all error checking - not
 * to mention documentation - happens there. To compile, see mvgc_makemex.m in
 * the utils subfolder */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double* a = mxGetPr(prhs[0]);
	double* z = mxGetPr(prhs[1]);
	
	const size_t n  = mxGetM(prhs[1]); 
	const size_t m  = mxGetN(prhs[1]);
	const mwSize d  = mxGetNumberOfDimensions(prhs[0]);
	const mwSize p  = (d < 3 ? 1 : mxGetDimensions(prhs[0])[2]);
	const size_t nn = n*n; 
    
    size_t t, k, i, j;
    double *xt, *xtk, *xti, *ak, *aki, w;
	
	double* x = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL)); /* allocate output array x  */
	memcpy(x,z,n*m*sizeof(double));                                  /* copy input z to output x */

	xt = x+n*p;
	for (t=p; t<m; ++t) {
		xtk = xt-n;
		ak = a;
		for (k=0; k<p; ++k) {
			xti = xt;
			aki = ak;
			for (i=0; i<n; ++i) {
				w = 0.0;
				for (j=0; j<n; ++j) w += *(aki+n*j) * *(xtk+j);
				*xti++ += w;
				++aki;
			}
			xtk -= n;
			ak += nn;
		}
		xt += n;
	}
}
