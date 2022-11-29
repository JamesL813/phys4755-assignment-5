///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Program file name: Root.java                                          //
//                                                                       //
// Tao Pang 2006                                                         //
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

// An example of searching a root via the secant method
// for f(x)=exp(x)*ln(x)-x*x=0.

import java.lang.*;

public class Root {
  public static void main(String argv[]) {
    double del = 1e-6, a = 1, b = 2;
    double dx = (b - a) / 10, x = (a + b) / 2;
    int n = 6;
    double[] r1 = {0.0, 0.0, 0.0};
    double[] r2 = {0.0, 0.0, 0.0};
    double[] r3 = {0.0, 0.0, 0.0};
    x = secant(n, del, x, dx);
    System.out.println("Root obtained: " + x);
  }

  // Method to carry out the secant search.

  public static double secant(int n, double del,
      double x, double dx) {
    int k = 0;
    double x1 = x + dx;
    while ((Math.abs(dx) > del) && (k < n)) {
      double d = f(x1) - f(x);
      double x2 = x1 - f(x1) * (x1 - x) / d;
      x = x1;
      x1 = x2;
      dx = x1 - x;
      k++;
    }
    if (k == n)
      System.out.println("Convergence not" +
          " found after " + n + " iterations");
    return x1;
  }

  // Method to provide function f(x)=exp(x)*log(x)-x*x.

  public static double f(int i, int j, double r[][]) {
    double e = 0.0;

    double epsilon0 = 0.0;

    double v0 = 1.09e3;
    double r0 = 0.321;

    return eta(i, j) * (e * e / (4 * Math.PI * epsilon0 * rDiff(i, j, r)))
        + delta(i, j, r) * v0 * Math.pow(-rDiff(i, j, r), r0);
  }

  private static int rDiff(double r1[], double r2[]) {
    return Math.sqrt(Math.pow(r1[0] + r2[0], 2) + Math.pow(r1[1] + r2[1], 2) + Math.pow(r1[2] + r2[2], 2));
  }

  private static int eta(int c1, int c2) {
    if (c1 * c2 < 0) {
      return -1;
    } else {
      return 1;
    }
  }

  private static int delta(double c1, double c2) {
    if (c1 * c2 < 0) {
      return 1;
    } else {
      return 0;
    }
  }
}
