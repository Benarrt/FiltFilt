#include <cstdint>
#include <vector>
#include <algorithm>
#include <iostream>

using vectord = std::vector<double>;

//  ****************************************************************************
//  ***                               DOUBLE                                 ***
//  ****************************************************************************

void CoreDoubleN(const double *X, uint32_t MX, uint32_t NX, const double *a, const double *b,
                 uint32_t order, double *Z, double *Y)
{
  // Direct form II transposed method for general filter length.
  // Implemented as time domain difference equations.
  // INPUT:
  //   X:  Double array. Operation happens of 1st dimension.
  //   MX: Number of elements in the 1st dimension
  //   NX: Number of columns, considers mutliple dimensions.
  //   a, b: Double vector of filter parameters. Both have nParam elements.
  //       The first element of a is 1.
  //   Z:  DOUBLE array, initial conditions.
  //   nParam: Number of filter parameters, order of filter + 1.
  // OUTPUT:
  //   Z:  DOUBLE array, final conditions.
  //   Y:  Double array, allocated by the caller.
   
  double Xi, Yi;
  uint32_t i, j, R;
  
  i = 0;
  while (NX--) {                         // Next slice
     R = i + MX;                         // End of the column
     while (i < R) {
        Xi = X[i];                       // Get signal
        Yi = b[0] * Xi + Z[0];           // Filtered value
        for (j = 1; j < order; j++) {    // Update conditions
           Z[j - 1] = b[j] * Xi + Z[j] - a[j] * Yi;
        }
        Z[order - 1] = b[order] * Xi - a[order] * Yi;
        
        Y[i++] = Yi;                      // Write to output
     }
     Z += order;                          // Next condition vector
  }
  
  return;
}

// =============================================================================
void CoreDouble2(const double *X, uint32_t MX, uint32_t NX, const double *a, const double *b,
                 double *Z, double *Y)
{
  // Filter with loop unrolled for 2 parameters (filter order 1).
  // Same input as the CoreDoubleN, but ommited [order], because it is 1.
   
  double Xi, Yi, z0, a1 = a[1];
  uint32_t i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi - a1 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
  }

  return;
}

// =============================================================================
void CoreDouble3(const double *X, uint32_t MX, uint32_t NX, const double *a, const double *b,
                 double *Z, double *Y)
{
  double Xi, Yi, z0, z1, a1 = a[1], a2 = a[2];
  uint32_t i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi      - a2 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
  }

  return;
}

// =============================================================================
void CoreDouble4(const double *X, uint32_t MX, uint32_t NX, const double *a, const double *b,
                 double *Z, double *Y)
{
  double Xi, Yi, z0, z1, z2, a1 = a[1], a2 = a[2], a3 = a[3];
  uint32_t i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi      - a3 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
  }

  return;
}

// =============================================================================
void CoreDouble5(const double *X, uint32_t MX, uint32_t NX, const double *a, const double *b,
                 double *Z, double *Y)
{
  double Xi, Yi, z0, z1, z2, z3, a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  uint32_t i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     z3 = Z[3];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi + z3 - a3 * Yi;
        z3 = b[4] * Xi      - a4 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
     *Z++ = z3;
  }

  return;
}

// =============================================================================
void CoreDouble6(const double *X, uint32_t MX, uint32_t NX, const double *a, const double *b,
                 double *Z, double *Y)
{
  double Xi, Yi, z0, z1, z2, z3, z4,
         a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  uint32_t i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     z3 = Z[3];
     z4 = Z[4];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi + z3 - a3 * Yi;
        z3 = b[4] * Xi + z4 - a4 * Yi;
        z4 = b[5] * Xi      - a5 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
     *Z++ = z3;
     *Z++ = z4;
  }

  return;
}

// =============================================================================
void CoreDouble7(const double *X, uint32_t MX, uint32_t NX, const double *a, const double *b,
                 double *Z, double *Y)
{
  // Still 33% faster than the loop method.
  double Xi, Yi, z0, z1, z2, z3, z4, z5,
         a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5], a6 = a[6];
  uint32_t i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     z3 = Z[3];
     z4 = Z[4];
     z5 = Z[5];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi + z3 - a3 * Yi;
        z3 = b[4] * Xi + z4 - a4 * Yi;
        z4 = b[5] * Xi + z5 - a5 * Yi;
        z5 = b[6] * Xi      - a6 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
     *Z++ = z3;
     *Z++ = z4;
     *Z++ = z5;
  }
  
  return;
}

// *****************************************************************************
// ***                             DOUBLE REVERSE                            ***
// *****************************************************************************

void CoreDoubleNR(const double *X, uint32_t MX, uint32_t NX, const double *a, const double *b,
                  uint32_t order, double *Z, double *Y)
{
  // Method for general filter length.
  // Signal X is process backwards, but a, b, and  Z have standard direction.
   
  double Xi, Yi;
  uint32_t j, R;
  int64_t i;
  
  R = 0;
  while (NX--) {
     i = R + MX - 1;
     while (i >= R) {
        Xi = X[i];
        Yi = b[0] * Xi + Z[0];
        for (j = 1; j < order; j++) {
           Z[j - 1] = b[j] * Xi + Z[j] - a[j] * Yi;
        }
        Z[order - 1] = b[order] * Xi - a[order] * Yi;
        
        Y[i--] = Yi;
     }
     Z += order;
     R += MX;
  }
  
  return;
}

void FilterX(const vectord& b, const vectord& a, const vectord& x, const vectord& z, bool reverse, vectord& outY, vectord& outZ)
{
   outY.resize(x.size());
   outZ = z;

   if(a.size() != b.size() || a[0] != 1)
   {
      std::cout << "Bad coefficents" << std::endl;
      return;
   }

   size_t order = a.size()-1;
   if(z.size() != order)
   {
      std::cout << "z.size() != order " << z.size() << " != " << order  << std::endl;
      return;
   }

   if (reverse)
   {
      CoreDoubleNR(x.data(), x.size(), 1, a.data(), b.data(), order, outZ.data(), outY.data());
   }
   else
   {
      switch (order) {
         case 1:   CoreDouble2(x.data(), x.size(), 1, a.data(), b.data(), outZ.data(), outY.data());  break;
         case 2:   CoreDouble3(x.data(), x.size(), 1, a.data(), b.data(), outZ.data(), outY.data());  break;
         case 3:   CoreDouble4(x.data(), x.size(), 1, a.data(), b.data(), outZ.data(), outY.data());  break;
         case 4:   CoreDouble5(x.data(), x.size(), 1, a.data(), b.data(), outZ.data(), outY.data());  break;
         case 5:   CoreDouble6(x.data(), x.size(), 1, a.data(), b.data(), outZ.data(), outY.data());  break;
         case 6:   CoreDouble7(x.data(), x.size(), 1, a.data(), b.data(), outZ.data(), outY.data());  break;
         default:  CoreDoubleN(x.data(), x.size(), 1, a.data(), b.data(), order, outZ.data(), outY.data());
      }
   }
}