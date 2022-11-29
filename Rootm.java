///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Program file name: Rootm.java                                         //
//                                                                       //
// Tao Pang 2006                                                       //
//                                                                       //
// Last modified: January 18, 2006                                       //
//                                                                       //
// (1) This Java program is part of the book, "An Introduction to        //
//     Computational Physics, 2nd Edition," written by Tao Pang and      //
//     published by Cambridge University Press on January 19, 2006.      //
//                                                                       //
// (2) No warranties, express or implied, are made for this program.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// An example of searching the root via the secant method
// for f_i[x_k] for i=1,2,...,n.

import java.lang.*;
public class Rootm {
  public static void main(String argv[]) {
    int n = 2;            // 2 variables
    int ni = 10;
    double del = 1e-6;    // tolerance
    double x[] = {1.5, 1.5}; // initial point
    secantm(ni, x, del);

    // Output the root obtained 
    System.out.println("The root is at x = " + x[0]
      + "; y = " + x[1]);
  }

  // Method to carry out the multivariable secant search.
  public static void secantm(int ni, double x[],
    double del) {
    int n = x.length;
    double h = 2e-5; // derivative step
    int index[] = new int[n];
    double a[][] = new double[n][n];

    int k = 0;
    double dx = 0.1; // initial newton step
    while ((Math.abs(dx)>del) && (k<ni)) { // ni max steps
      double b[] = f(x);
      for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
          double hx = x[j]*h;
          x[j] += hx;
          double c[] = f(x);
          a[i][j] = (c[i]-b[i])/hx; // Discrete Jacobian
        }
      }
      for (int i=0; i<n; ++i) b[i] = -b[i];
      double d[] = solve(a, b, index); // Find newton step
      dx = 0;
      for (int i=0; i<n; ++i) {
        dx += d[i]*d[i]; // update newton step
        x[i] += d[i];
      }
      dx = Math.sqrt(dx/n);
      k++;
    }

    // bug
    if (k==n) {
      System.out.println("Convergence not" +
      " found after " + ni + " iterations");
    }
  }

  // Method to solve the equation a[][] x[] = b[] with
  // the partial-pivoting Gaussian elimination scheme.
  public static double[] solve(double a[][], double b[],
    int index[]) {
    int n = b.length;
    double x[] = new double[n];

    // Invoke the partial-pivoting Gaussian elimination
    gaussian(a, index);

    // Rescale array b[i] with the ratios stored
    for(int i=0; i<n-1; ++i) {
      for(int j =i+1; j<n; ++j) {
        b[index[j]] -= a[index[j]][i]*b[index[i]];
      }
    }

    // Perform the backward substitutions
    x[n-1] = b[index[n-1]]/a[index[n-1]][n-1];
    for (int i=n-2; i>=0; --i) {
      x[i] = b[index[i]];
      for (int j=i+1; j<n; ++j) {
        x[i] -= a[index[i]][j]*x[j];
      }
      x[i] /= a[index[i]][i];
    }
    return x;
  }

  // Method to perform the partial-pivoting Gaussian
  // elimination.  index[] records the pivoting order.
  public static void gaussian(double a[][], int index[]) {
    int n = index.length;
    double c[] = new double[n];

    // Initialize the index
    for (int i=0; i<n; ++i) index[i] = i;

    // Find the rescaling factors, one from each row
    for (int i=0; i<n; ++i) {
      double c1 = 0;
      for (int j=0; j<n; ++j) {
        double c0 = Math.abs(a[i][j]);
        if (c0 > c1) c1 = c0;
      }
      c[i] = c1;
    }

    // Search the pivoting element from each column
    int k = 0;
    for (int j=0; j<n-1; ++j) {
      double pi1 = 0;
      for (int i=j; i<n; ++i) {
        double pi0 = Math.abs(a[index[i]][j]);
        pi0 /= c[index[i]];
        if (pi0 > pi1) {
          pi1 = pi0;
          k = i;
        }
      }

      // Interchange rows according to the pivoting order
      int itmp = index[j];
      index[j] = index[k];
      index[k] = itmp;
      for (int i=j+1; i<n; ++i) {
        double pj = a[index[i]][j]/a[index[j]][j];

        // Record pivoting ratios below the diagonal
        a[index[i]][j] = pj;

        // Modify other elements accordingly
        for (int l=j+1; l<n; ++l)
          a[index[i]][l] -= pj*a[index[j]][l];
      }
    }
  }

  // Method to provide function f_i[x_k].
  public static double[] f(double x[]) {
    double fx[] = new double[2];
    fx[0] = Math.exp(x[0]*x[0])*Math.log(x[1])-x[0]*x[0];
    fx[1] = Math.exp(x[1])*Math.log(x[0])-x[1]*x[1];
    return fx;
  }
}
