#include <math.h>
#include "mex.h"
#define EPS  1.0e-6

/*
  
  NOTE: as of Feb 7, 2006, this is better than SurfaceGeom's 

  Inputs from matlab: faces -- a Fx3 array
		      vert -- a Vx3 array
                      vertVals -- a Vx1 array
  Outputs: 6 Sx1 arrays

*/

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  /* ************** Inputs *******************/
  double *faces, *vert, *vertVals;
  
  /* ************** Outputs *******************/
  double *Xs, *Ys, *Zs;
  double *Xe, *Ye, *Ze;

  /* ************** Functions **************** */
  void findContour(double *faces, 
		   double *vert,
		   double *vertVals,
		   int F, 
		   int V, 
		   int *p_NumSegments,
		   int maxNumSegments,
		   double *contourLarge);

  /* ************** Others ******************* */
  int counter;
  int ss;
  int V, F;
  int numSegments;
  int maxNumSegments;
  double *contourLarge;

  /* Inputs */
  faces = mxGetPr(prhs[0]);
  vert = mxGetPr(prhs[1]);
  vertVals = mxGetPr(prhs[2]);

  /* Get dimensions */
  F = mxGetM(prhs[0]);
  V = mxGetM(prhs[1]);
  maxNumSegments = F; 
  
  contourLarge = mxMalloc( 2*3*maxNumSegments*sizeof(double) );

  findContour(faces,
	      vert,
	      vertVals,
	      F, 
	      V, 
	      &numSegments,
	      maxNumSegments,
	      contourLarge);

  /* Outputs */
  plhs[0] = mxCreateDoubleMatrix(numSegments, 1, mxREAL);
  Xs = mxGetPr(plhs[0]);

  plhs[1] = mxCreateDoubleMatrix(numSegments, 1, mxREAL);
  Ys = mxGetPr(plhs[1]);

  plhs[2] = mxCreateDoubleMatrix(numSegments, 1, mxREAL);
  Zs = mxGetPr(plhs[2]);

  plhs[3] = mxCreateDoubleMatrix(numSegments, 1, mxREAL);
  Xe = mxGetPr(plhs[3]);

  plhs[4] = mxCreateDoubleMatrix(numSegments, 1, mxREAL);
  Ye = mxGetPr(plhs[4]);

  plhs[5] = mxCreateDoubleMatrix(numSegments, 1, mxREAL);
  Ze = mxGetPr(plhs[5]);

  /* Copy output into Xs, etc */
  counter = 0;
  for(ss = 0; ss < numSegments; ss++)
    {
      Xs[ss] = contourLarge[counter++];
      Ys[ss] = contourLarge[counter++];
      Zs[ss] = contourLarge[counter++];
      
      Xe[ss] = contourLarge[counter++];
      Ye[ss] = contourLarge[counter++];
      Ze[ss] = contourLarge[counter++];
    }

  /* Free up contourLarge */
  mxFree(contourLarge);

}

void findContour(double *faces, 
		 double *vert,
		 double *vertVals,
		 int F, 
		 int V, 
		 int *p_NumSegments,
		 int maxNumSegments,
		 double *contourLarge)
{

  int t;
  int nA, nB, nC;
  int counter = 0;
  double s0, diff;
  double vec[3];
  int numHits;
  int zeroVert[3];

  for(t = 0; t < F; t++)
    {
      if(counter >= 2*3*maxNumSegments)
	{
	  mexPrintf("\nIn findContourMEX.c: need to allocate more space");
	  mexPrintf(" for contourLarge[]\n");
	  mexPrintf(" -- try increasing maxNumSegments and recompiling.\n\n");
	  /* mexErrMsgTxt("Exiting ..."); */ /* uncomment if not R14 sp1 */
	}
      
      nA = faces[t]       - 1;
      nB = faces[t + F]   - 1;
      nC = faces[t + F*2] - 1;
      
      if( mxIsFinite(vertVals[nA]) &&
	  mxIsFinite(vertVals[nB]) &&
	  mxIsFinite(vertVals[nC]) )
	{
	  numHits = 0;

	  if(vertVals[nA]*vertVals[nB] < 0)
	    {
              diff = vertVals[nA] - vertVals[nB]; 
              if( fabs(diff) < EPS )
                {
                  s0 = 0.5; 
                }
              else
                {
                  s0 = vertVals[nA]/diff;
                }

	      vec[0] = vert[nA      ] + s0*(vert[nB      ] - vert[nA      ]);
	      vec[1] = vert[nA + V  ] + s0*(vert[nB + V  ] - vert[nA + V  ]);
	      vec[2] = vert[nA + V*2] + s0*(vert[nB + V*2] - vert[nA + V*2]);
	      
	      contourLarge[counter++] = vec[0]; 
	      contourLarge[counter++] = vec[1]; 
	      contourLarge[counter++] = vec[2]; 

	      numHits++;
	    }
	  
	  if(vertVals[nA]*vertVals[nC] < 0)
	    {
              diff = vertVals[nA] - vertVals[nC];
              if( fabs(diff) < EPS )
                {
                  s0 = 0.5; 
                }
              else
                {
                  s0 = vertVals[nA]/diff;
                }
	      
	      vec[0] = vert[nA      ] + s0*(vert[nC      ] - vert[nA      ]);
	      vec[1] = vert[nA + V  ] + s0*(vert[nC + V  ] - vert[nA + V  ]);
	      vec[2] = vert[nA + V*2] + s0*(vert[nC + V*2] - vert[nA + V*2]);
	      
	      contourLarge[counter++] = vec[0]; 
	      contourLarge[counter++] = vec[1]; 
	      contourLarge[counter++] = vec[2]; 
	      
	      numHits++;
	    }
	  
	  if(vertVals[nB]*vertVals[nC] < 0)
	    {
              diff = vertVals[nB] - vertVals[nC];
              if( fabs(diff) < EPS )
                {
                  s0 = 0.5;
                }
              else
                {
                  s0 = vertVals[nB]/diff;
                }
	      
	      vec[0] = vert[nB      ] + s0*(vert[nC      ] - vert[nB      ]);
	      vec[1] = vert[nB + V  ] + s0*(vert[nC + V  ] - vert[nB + V  ]);
	      vec[2] = vert[nB + V*2] + s0*(vert[nC + V*2] - vert[nB + V*2]);
	      
	      contourLarge[counter++] = vec[0]; 
	      contourLarge[counter++] = vec[1]; 
	      contourLarge[counter++] = vec[2]; 
	      
	      numHits++;
	    }

          /* Take care of the 0 + - case */
          if( numHits == 1 )
            {
              if( vertVals[nA] == 0 )
                {
		  contourLarge[counter++] = vert[nA      ];
		  contourLarge[counter++] = vert[nA + V  ];
		  contourLarge[counter++] = vert[nA + V*2];
                  numHits++;
                }
              else if( vertVals[nB] == 0 )
                {
		  contourLarge[counter++] = vert[nB      ];
		  contourLarge[counter++] = vert[nB + V  ];
		  contourLarge[counter++] = vert[nB + V*2];
                  numHits++;
                }
              else if( vertVals[nC] == 0 )
                {
		  contourLarge[counter++] = vert[nC      ];
		  contourLarge[counter++] = vert[nC + V  ];
		  contourLarge[counter++] = vert[nC + V*2];
                  numHits++;
                }
            }

          /* Take care of the 0 0 +- case (not efficient, I know) */
          if( numHits == 0 )
            {
              if( (vertVals[nA] == 0) && (vertVals[nB] == 0) )
                {
		  contourLarge[counter++] = vert[nA      ];
		  contourLarge[counter++] = vert[nA + V  ];
		  contourLarge[counter++] = vert[nA + V*2];
                  numHits++;

		  contourLarge[counter++] = vert[nB      ];
		  contourLarge[counter++] = vert[nB + V  ];
		  contourLarge[counter++] = vert[nB + V*2];
                  numHits++;
                }
              else if( (vertVals[nA] == 0) && (vertVals[nC] == 0) )
                {
		  contourLarge[counter++] = vert[nA      ];
		  contourLarge[counter++] = vert[nA + V  ];
		  contourLarge[counter++] = vert[nA + V*2];
                  numHits++;

		  contourLarge[counter++] = vert[nC      ];
		  contourLarge[counter++] = vert[nC + V  ];
		  contourLarge[counter++] = vert[nC + V*2];
                  numHits++;
                }
              else if( (vertVals[nB] == 0) && (vertVals[nC] == 0) )
                {
		  contourLarge[counter++] = vert[nB      ];
		  contourLarge[counter++] = vert[nB + V  ];
		  contourLarge[counter++] = vert[nB + V*2];
                  numHits++;

		  contourLarge[counter++] = vert[nC      ];
		  contourLarge[counter++] = vert[nC + V  ];
		  contourLarge[counter++] = vert[nC + V*2];
                  numHits++;
                }
            }

          if( (numHits == 1) || (numHits == 3) )
            {
              mexPrintf("\nIn findContourMEX.c: numHits should be 0 or 2!");
              /* mexErrMsgTxt("Exiting ..."); */
            }

	}/* if the val for each vertex in the triangle is finite */

    }/* for each triangle */

      *p_NumSegments = counter/6;

}

void codeGraveyard(double *faces, 
                   double *vert,
                   double *vertVals,
                   int F, 
                   int V, 
                   int *p_NumSegments,
                   int maxNumSegments,
                   double *contourLarge)
{

  int t;
  int nA, nB, nC;
  int counter = 0;
  double s0;
  double vec[3];
  int numHits;
  int zeroVert[3];

  /* --------------------------------------------------- */
  /* If we have a situation where the vertex value is basically zero */
  /* on two vertices (but not on the third), the idea here was to draw */
  /* the edge between the two zero-value vertices */
  /* --------------------------------------------------- */

  /* if numHits == 0, check for the case of 2 zero verts */
  if(numHits == 0)
    {
      if( fabs(vertVals[nA]) < EPS )
        zeroVert[numHits++] = nA;
          
      if( fabs(vertVals[nB]) < EPS )
        zeroVert[numHits++] = nB;
          
      if( fabs(vertVals[nC]) < EPS )
        zeroVert[numHits++] = nC;
          
      if(numHits == 2)
        {
          nA = zeroVert[0];
          nB = zeroVert[1];
          
          contourLarge[counter++] = vert[nA      ];
          contourLarge[counter++] = vert[nA + V  ];
          contourLarge[counter++] = vert[nA + V*2];
          
          contourLarge[counter++] = vert[nB      ];
          contourLarge[counter++] = vert[nB + V  ];
          contourLarge[counter++] = vert[nB + V*2];
        }
          
    } /* if numHits == 0, check for the case of 2 zero verts */
	  
}  
