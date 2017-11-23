package util;
import java.util.*;
import static java.lang.Math.sqrt;
public class Math {
    public static double standardDeviation(double a[]) {
        int n = a.length;
        if(n == 0)
        return 0.0;
        double sum = 0;
        double sq_sum = 0;
        for(int i = 0; i < n; ++i) {
            sum += a[i];
            sq_sum += a[i] * a[i];
        }
        double mean = sum / n;
        double variance = sq_sum / n - mean * mean;
        return sqrt(variance);
    }
}
