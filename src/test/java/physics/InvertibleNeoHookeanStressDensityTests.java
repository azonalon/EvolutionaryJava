package physics;
import org.junit.Test;
import static physics.InvertibleNeoHookeanStressDensity.*;
import static org.junit.Assert.*;
public class InvertibleNeoHookeanStressDensityTests {

    @Test
    public void computeExamples() {
        double[] r = computeStressDifferentialComponents(0.3, 0.5, 1.3, 1.2, 0.3);
        double[] r2 = computeStressGradientComponents(0.3, 0.5, 1.3, 1.2, 0.3);
        double[] t = new double[] {3.44706684599153,18.7689664089106,0.0816001282036203,24.5276960027332,-0.448611202750378,14.0356565981633};
        double[] t2 = new double[] {3.44706684599153,-9.36000316599424,-5.58988985857141};
        assertArrayEquals(t, r, 1e-10);
        assertArrayEquals(t2, r2, 1e-10);

        t = new double[] {2.9585392278565,13.7223024624783,1.14118685408538,3.83247641114176,-3.07712305191825,46.3684960624857};
        t2 = new double[] {2.9585392278565,-1.60327363841031,-13.494065091661};
        r = computeStressDifferentialComponents(1.0, 0.2, 1.3, 1.2, 0.3);
        r2 = computeStressGradientComponents(1.0, 0.2, 1.3, 1.2, 0.3);
        assertArrayEquals(t, r, 1e-10);
        assertArrayEquals(t2, r2, 1e-10);

        t = new double[] {2.9585392278565,13.7223024624782,1.14118685408533,46.3684960624857,-3.07712305191825,3.83247641114176};
        t2 = new double[] {2.9585392278565,-13.494065091661,-1.60327363841031};
        r = computeStressDifferentialComponents(0.2, 1, 1.3, 1.2, 0.3);
        r2 = computeStressGradientComponents(0.2, 1, 1.3, 1.2, 0.3);
        assertArrayEquals(t, r, 1e-10);
        assertArrayEquals(t2, r2, 1e-10);

        t = new double[] {11.9055855463672,1244.84315670161,-1210.09421816975,26.7939884373991,-7.83297749968835,27.0390187815828};
        r = computeStressDifferentialComponents(0.01, 0.002, 1.3, 1.2, 0.3);
        r2 = computeStressGradientComponents(0.01, 0.002, 1.3, 1.2, 0.3);
        t2 = new double[] {11.9055855463672,-14.5906284951007,-14.8686200033555};
        assertArrayEquals(t, r, 1e-10);
        assertArrayEquals(t2, r2, 1e-10);

        t = new double[] {12.096443131052, 1e10, -1e10,
            27.6203475133729,
            -7.93830358922197,
            26.4064283369386};
        r = computeStressDifferentialComponents(-0.02, 0.02, 1.3, 1.2, 0.3);
        r2 = computeStressGradientComponents(-0.02, 0.02, 1.3, 1.2, 0.3);
        t2 = new double[] {12.096443131052,-15.5487202819036,-14.1511871090737};
        assertArrayEquals(t, r, 1e-10);
        assertArrayEquals(t2, r2, 1e-10);
    }
}
