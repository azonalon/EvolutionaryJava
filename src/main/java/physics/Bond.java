package physics;
import physics.*;
import util.Math.*;

public class Bond {
    Cell a, b;
    double phi0;
    int rotationCounter;
    double k, E, D, c, l;
    double angleA, angleB;

    public static Bond harmonicAverageBond(Cell a, Cell b, double angleA, double angleB) {
        Bond bond = new Bond();
        double dx  = - a.x + b.x;
        double dy  = - a.y + b.y;
        a.theta = Math.atan2(dy, dx) + angleA;
        b.theta = Math.atan2(dy, dx) + angleB;
        bond.a=a;
        bond.b=b;
        bond.angleA = angleA;
        bond.angleB = angleB;
        bond.k = util.Math.harmonicMean(a.k, b.k);
        bond.E = util.Math.harmonicMean(a.E, b.E);
        bond.D = util.Math.harmonicMean(a.D, b.D);
        bond.c = util.Math.harmonicMean(a.c, b.c);
        bond.l = a.r + b.r;
        bond.rotationCounter = 0;
        bond.phi0 = 0;
        return bond;
    }
}
