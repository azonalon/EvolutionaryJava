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
    public static double harmonicMean(double a, double b) {
        return 2 * a * b / (a + b + 1e-9);
    }
}
