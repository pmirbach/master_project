#include "fintrf.h"
C======================================================================
C     matri1.f
C     Computational function that takes a matrix and doubles the elements
C     This is a MEX file for MATLAB.
C======================================================================

C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetPr
      mwPointer mxGetM, mxGetN

C     Pointers to input/output mxArrays:
      mwPointer x_ptr, y_ptr

C     Array information:
      mwPointer mrows, ncols
      mwSize size


C     Arguments for computational routine:
      real*8  x_input, y_output


C-----------------------------------------------------------------------
C     Read Data 
C     Get the size of the input array.
      mrows = mxGetM(prhs(1))
      ncols = mxGetN(prhs(1))
      size = mrows*ncols

C     Create Fortran array from the input argument.
      x_ptr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(x_ptr,x_input,size)

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(ncols,mrows,0)
      y_ptr = mxGetPr(plhs(1))


C     Call the computational subroutine.
      call timestwo(y_output, x_input)

C     Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(y_output,y_ptr,size)     

      return
      end

C-----------------------------------------------------------------------
C     Computational routine

      subroutine matr1(mat_out,mat_in,m,n)

      mwSize m, n
      real*8  mat_in(m,n), mat_out(m,n)
      mwSize i, j

      do i=1,m

         do j=1,n

            mat_out(m,n) = 2.0 * mat_in(i,j)

         end do

      end do

      return
      end