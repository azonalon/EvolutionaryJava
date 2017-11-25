package physics;
import physics.*;
import util.Math.*;

public class Bond {
    Cell a, b;
    double phi0;
    int rotationCounter=0;
    double k, E, D, c, l;
    double unstressedAngle=0;

    public static Bond harmonicAverageBond(Cell a, Cell b, double unstressedAngle) {
        Bond bond = new Bond();
        double dx  = - a.x + b.x;
        double dy  = - a.y + b.y;
        a.theta = b.theta = Math.atan2(dy, dx) + unstressedAngle;
        bond.a=a;
        bond.b=b;
        bond.unstressedAngle = unstressedAngle;
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
