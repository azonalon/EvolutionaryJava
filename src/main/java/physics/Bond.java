package physics;
import physics.*;
import util.Math.*;
// import Math;

public class Bond {
    Cell a, b;
    double phi0;
    int rotationCounter;
    double k, E, D, c, l;
    double angle;

    public static Bond harmonicAverageBond(Cell a, Cell b, double angle) {
        Bond bond = new Bond();
        double dx  = - a.x + b.x;
        double dy  = - a.y + b.y;
        bond.phi0 = Math.atan2(dy, dx);
        a.theta = 0;
        b.theta = 0;
        bond.angle  = angle;
        bond.a=a;
        bond.b=b;
        bond.k = util.Math.harmonicMean(a.k, b.k);

        bond.E = util.Math.harmonicMean(a.E, b.E);
        bond.D = util.Math.harmonicMean(a.D, b.D);
        bond.c = util.Math.harmonicMean(a.c, b.c);
        bond.l = a.r + b.r;
        bond.rotationCounter = 0;
        return bond;
    }
    public String toString() {
        return "Bond(\na="+a +",\nb="+b+ ")";
    }
}
