package physics;
import physics.*;
import util.Math.*;
// import Math;

public class Bond {
    Cell a, b;
    double nu0;
    int rotationCounter;
    double k, E, alpha, zeta, l;
    double angle;

    public static Bond harmonicAverageBond(Cell a, Cell b, double angle) {
        Bond bond = new Bond();
        double dx  = - a.x + b.x;
        double dy  = - a.y + b.y;
        bond.nu0 = Math.atan2(dy, dx);
        a.theta = 0;
        b.theta = 0;
        bond.angle  = angle;
        bond.a=a;
        bond.b=b;
        bond.k = util.Math.harmonicMean(a.k, b.k);

        bond.E = util.Math.harmonicMean(a.E, b.E);
        bond.zeta = util.Math.harmonicMean(a.zeta, b.zeta);
        bond.alpha = util.Math.harmonicMean(a.alpha, b.alpha);
        bond.l = a.r + b.r;
        bond.rotationCounter = 0;
        return bond;
    }
    public String toString() {
        return String.format("Bond(a=%d, b=%d, phi=%f, k=%f, E=%f, l=%f)",
                            a.index, b.index, nu0, k, E, l);
    }
}
