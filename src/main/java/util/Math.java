package util;
import java.util.*;
import static java.lang.Math.sqrt;
import static java.lang.Math.PI;
import static java.lang.Math.hypot;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import org.ejml.data.*;
import static org.ejml.dense.fixed.CommonOps_DDF2.*;

public class Math {
    public static final String str(DMatrix2x2 m) {
        return String.format("((%g, %g), (%g, %g))", m.a11, m.a12, m.a21, m.a22);
    }
    public static final void add(double a, DMatrix2x2 A, double b, DMatrix2x2 B,
                     DMatrix2x2 dest) {
        dest.a11 = a * A.a11 + b * B.a11;
        dest.a12 = a * A.a12 + b * B.a12;
        dest.a21 = a * A.a21 + b * B.a21;
        dest.a22 = a * A.a22 + b * B.a22;
    }
    public static final boolean invertTranspose( DMatrix2x2 a , DMatrix2x2 inv ) {
        double scale = 1.0/elementMaxAbs(a);

        double a11 = a.a11*scale;
        double a12 = a.a12*scale;
        double a21 = a.a21*scale;
        double a22 = a.a22*scale;

        double m11 = a22;
        double m12 = -( a21);
        double m21 = -( a12);
        double m22 = a11;

        double det = (a11*m11 + a12*m12)/scale;

        inv.a11 = m11/det;
        inv.a21 = m21/det;
        inv.a12 = m12/det;
        inv.a22 = m22/det;

        return !Double.isNaN(det) && !Double.isInfinite(det);
    }

    public static final double[] singularValueDecomposition(DMatrix2x2 m) {
        double u1, u2, x, y, v1, v2, e, f, g, h, Q, R, a1, a2, th, phi;
        e = (m.a11 + m.a22)/2;
        f = (m.a11 - m.a22)/2;
        g = (m.a21 + m.a12)/2;
        h = (m.a21 - m.a12)/2;
        Q = hypot(e, h);
        R = hypot(f, g);
        x = Q + R;
        y = Q - R;
        a1 = atan2(g, f);
        a2 = atan2(h, e);
        th = (a2 - a1)/2;
        phi = (a2 + a1)/2;
        u1 = cos(phi);
        u2 = sin(phi);
        v1 = cos(th);
        v2 = sin(th);
        return new double[] {u1, u2, x, y, v1, v2};
    }
    public static double standardDeviation(double a[]) {
        int n = a.length;
        double sum = 0;
        double sq_sum = 0;
        for(int i = 0; i < n; ++i) {
            sum += a[i];
        }
        double mean = sum / n;
        for(int i = 0; i < n; ++i) {
            sq_sum += (a[i]-mean) * (a[i]-mean);
        }
        double variance = sq_sum / n;
        return sqrt(variance);
    }
    double mod(double x, double n) {
        return (x % n) - (x < 0 ? n : 0);
    }
    static public double circleMod(double phi) {
        if(phi < -PI) {
            return phi + 2 * PI;
        } else if(phi > PI) {
            return phi - 2 * PI;
        }
        return phi;
    }
    public static double average(double a[]) {
        int n = a.length;
        double sum = 0;
        for(int i = 0; i < n; ++i) {
            sum += a[i];
        }
        return sum / n;
    }

    public static double[] getColumn(double[][] a, int c) {
        // int m = a[0].length;
        double[] col = new double[a.length];
        for(int i = 0; i < a.length; i++) {
            col[i] = a[i][c];
        }
        return col;
    }

    public static double harmonicMean(double a, double b) {
        return 2 * a * b / (a + b + 1e-9);
    }
    public static double sqrd(double a) {
        return a*a;
    }
}
