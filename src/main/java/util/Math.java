package util;
import java.util.*;
import static java.lang.Math.sqrt;
public class Math {
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
