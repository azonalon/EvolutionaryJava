package physics;
import physics.*;
import org.junit.Test;
import org.junit.Rule;
import org.junit.rules.*;
import static org.junit.Assert.*;
import static physics.SoftBody.*;
import static java.lang.Math.PI;
import static physics.Cell.dt;
import static physics.Cell.t;
import static util.Math.*;

public class TestMathFunctions {
    double delta = 1e-4;

    @Test
    public void testInterpolate() {
        double a=-5, fa=16, dfa=-8,
               b=6, fb=49, dfb=14;
        double min = -1;

        assertEquals(min, interpolate(a, fa, dfa, b, fb, dfb), delta);
        assertNotEquals(min + 3, interpolate(a, fa, dfa, b, fb, dfb), delta);
        assertEquals(min, interpolate(b, fb, dfb, a, fa, dfa),delta);
        assertNotEquals(min + 3, interpolate(b, fb, dfb, a, fa, dfa), delta);

        a=1.64852; fa=6.03818; dfa=5.59878;
               b=0.48813; fb=2.17631; dfb=1.66492;
        min = Double.NaN;
        assertEquals(min, interpolate(b, fb, dfb, a, fa, dfa),delta);
        assertNotEquals(3, interpolate(b, fb, dfb, a, fa, dfa), delta);

        assertEquals(min, interpolate(a, fa, dfa, b, fb, dfb), delta);
        assertNotEquals(3, interpolate(a, fa, dfa, b, fb, dfb), delta);

        a=0.080544; fa=-0.94141; dfa=0.300789;
               b=1.39375; fb=-0.88757; dfb=0.943421;
        min = 0.987711;
        assertEquals(min, interpolate(b, fb, dfb, a, fa, dfa),delta);
        assertNotEquals(3, interpolate(b, fb, dfb, a, fa, dfa), delta);

        assertEquals(min, interpolate(a, fa, dfa, b, fb, dfb), delta);
        assertNotEquals(3, interpolate(a, fa, dfa, b, fb, dfb), delta);

    }
}
