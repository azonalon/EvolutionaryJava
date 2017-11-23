package physics;
import physics.*;
import java.util.*;
// import util.Math;
import java.util.function.*;
import static java.lang.Math.PI;
import static java.lang.Math.sqrt;
// import static System.out.format;

public class SoftBody {

    SoftBody(Cell[] cells, Bond[] bonds) {
        this.cells = cells;
        this.bonds = bonds;
    }

    static Bond harmonicAverageBond(Cell a, Cell b, double unstressedAngle) {
        Bond bond = new Bond();
        double dx  = - a.x + b.x;
        double dy  = - a.y + b.y;
        a.theta = b.theta = Math.atan2(dy, dx) + unstressedAngle;
        bond.a=a;
        bond.b=b;
        bond.unstressedAngle = unstressedAngle;
        bond.k = harmonicMean(a.k, b.k);
        bond.E = harmonicMean(a.E, b.E);
        bond.D = harmonicMean(a.D, b.D);
        bond.c = harmonicMean(a.c, b.c);
        bond.l = a.r + b.r;
        bond.rotationCounter = 0;
        bond.phi0 = 0;
        return bond;
    }

    Cell[] cells;
    Bond[] bonds;
    double energy = 0;
    static double ka=0, kb=0, kc=0.2, kd=0.2;
    static double[] energies;
    static int stepCounter=0;

    /*
     * STATIC methods
     */
    static double sqrd(double a) {
        return a*a;
    }

    static class Bond {
        Cell a, b;
        double phi0;
        int rotationCounter=0;
        double k, E, D, c, l;
        double unstressedAngle=0;
    }

    static double phi0;
    static int phiCounter = 0;
    static int LEFT = 1 << 0;
    static int RIGHT= 1 << 1;
    static int UP   = 1 << 2;
    static int DOWN = 1 << 3;
    static int counter=0;
    static double dt = 0.1;
    static boolean startProgram = true;
    static Function<Cell, String> cellStatusReading =
        (Cell c)->{return String.format("%f ", c.x);};

    static Function<SoftBody, String> bodyStatusReading =
        (body)->{return body.totalEnergy() + "";};


    static double circleMod(double angle) {
        return angle;
    }
    static double harmonicMean(double a, double b) {
        return 2 * a * b / (a + b + 1e-9);
    }

    /**
     * checks if the bond forces between the cells have been added to each cell.
     * this is a symmetric operation a <--> b
     * @param  Cell a             [description]
     * @param  Cell b             [description]
     * @return      [description]
     */

    void propagateCells() {
        for(Cell c: cells) {
            System.out.print(cellStatusReading.apply(c) + " ");
            energy += c.L * c.L * c.I /2.0 + c.vx * c.vx * c.m /2.0 + c.vy * c.vy * c.m /2.0;
            c.propagate();
        }
        System.out.print(bodyStatusReading.apply(this));
        System.out.println();
        stepCounter++;
    }

    void updateForces() {
        energy=0;
        for(Cell c: cells) {
            c.fX=0;
            c.fY=0;
            c.T=0;
        }
        for(Bond bond: bonds) {
            addBondForces(bond);
        }
    }

    void addBondForces(Bond b) {
        Cell first = b.a;
        Cell second = b.b;
        double dx  = - first.x + second.x;
        double dy  = - first.y + second.y;
        double d = Math.sqrt(dx*dx + dy*dy);

        double l=b.l, E=b.E, k=b.k, c=b.c, D=b.D;

        dx = dx / d;
        dy = dy / d;
        double phi = Math.atan2(dy, dx);
        if(b.phi0 * phi < -PI) {
            if(b.phi0 > phi) {
                b.rotationCounter++;
            } else {
                b.rotationCounter--;
            }
        }
        b.phi0 = phi;
        phi += - b.unstressedAngle + b.rotationCounter * 2 * PI;
        double fShear = - 6 * l * l * E * (first.theta + second.theta - 2 * phi);
        assert(d>0.01);


        first.fX  += (d - l) * k * dx + fShear * (-1*dy);
        first.fY  += (d - l) * k * dy + fShear * (   dx);
        second.fX -= (d - l) * k * dx + fShear * (-1*dy);
        second.fY -= (d - l) * k * dy + fShear * (   dx);

        first.T  -= 2 * E * l * l * l * (2 * first.theta + second.theta - 3 * phi);
        second.T -= 2 * E * l * l * l * (2 * second.theta + first.theta - 3 * phi);


        energy += (d-l)*(d-l)* k/2.0;
        energy += ka * l * l * l * E * sqrd(first.theta + second.theta - 2 * phi);
        energy += kb * l * l * E * sqrd(first.theta + second.theta - 2 * phi);
        energy += kc * E * l * l * l * sqrd(2 * first.theta + second.theta - 3 * phi);
        energy += kd * E * l * l * l * sqrd(2 * second.theta + first.theta - 3 * phi);
        // kc = 0.2
        // kd = 0.2

        if(startProgram == true) {
            System.err.format("Parameters: k=%f, E=%f, D=%f, c=%f, l=%f, m1=%f, m2=%f\n",
                                k, E, D, c, l, first.m, second.m);
            startProgram = false;
        }


        // damping forces
        double vXRelative = first.vx - second.vx;
        double vYRelative = first.vy - second.vy;
        // subtract parts orthogonal to direction vector
        double orthogonalFraction = vXRelative * -1 * dy + vYRelative * dx;
        double dphicalc = -Math.asin(orthogonalFraction * dt/d);
        double shearDamping = 0.05 * D * (first.L - dphicalc + second.L - dphicalc);
        orthogonalFraction -= shearDamping;
        vXRelative -=  orthogonalFraction * -1 * dy;
        vYRelative -=  orthogonalFraction * dx;
        // System.err.format("orth vel.: %f, orth Damping: %f, Torque: %f, phi %f\n",orthogonalFraction, shearDamping, second.T, phi);
        // System.err.format("fShear: %f, d: %f, dphi: %f, vxR: %f, vyR: %f\n", fShear, d, dphicalc, vXRelative, vYRelative);
        // System.err.format("dx: %f, dy: %f, phi: %f\n", dx, dy, phi);
        System.err.format("th1: %f, th2: %f, phi: %f\n", first.theta, second.theta, phi);

        first.fX -= c * vXRelative;
        first.fY -= c * vYRelative;
        second.fX += c * vXRelative;
        second.fY += c * vYRelative;
        first.T -= D * (first.L - second.L);
        second.T -= D * (second.L - first.L);
    }

    public double totalEnergy() {
        return energy;
    }

    public static void main(String[] args) {
        Cell[] cells = null;
        Bond[] bonds = null;

        if("symmetric rotation oscillation".equals(args[0])) {
            Cell a, b;
            bodyStatusReading = (body) -> ("" + body.totalEnergy());
            cellStatusReading = (c) -> {
                return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r   , E,
                           1, 1, 1, 1, 0.5 , 1,
                           //x, y,vx,vy, L
                           -0.5, 0, 0, -1, -0.0 * 2 * PI);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r, E,
                             1, 1, 1, 1, 0.5, 0.1,
                           //x, y,vx,vy, L
                             0.5, 0, 0.0, 1, 0.0 * 2*PI);
            cells = new Cell[] {a, b};
            bonds = new Bond[] {harmonicAverageBond(a, b, 0.0)};
            dt = Double.parseDouble(args[1]);
            int nSteps = Integer.parseInt(args[2]);
            testSimulation(cells, bonds, nSteps);
        } else if("beam oscillation rotation".equals(args[0])) {
            Cell a, b;
            bodyStatusReading = (body) -> ("" + body.totalEnergy());
            cellStatusReading = (c) -> {
                return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r   , E, index
                           10000, 10000, 1, 1, 0.5 , 1,
                           //x, y,vx,vy, L
                           0, 0, 0, 0, -0.0 * 2 * PI);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r, E,
                             1, 1, 1, 1, 0.5, 1,
                           //x, y,vx,vy, L
                             0, -1, 0.3, 0.0, 0.0 * 2*PI);
            cells = new Cell[] {a, b};
            bonds = new Bond[] {harmonicAverageBond(a, b, 0.0)};
            dt = Double.parseDouble(args[1]);
            int nSteps = Integer.parseInt(args[2]);
            testSimulation(cells, bonds, nSteps);
        } else if("beam oscillation".equals(args[0])) {
            Cell a, b;
            bodyStatusReading = (body) -> ("" + body.totalEnergy());
            cellStatusReading = (c) -> {
                return String.format("%f %f %f ", c.y, c.theta, c.L);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r   , E, index
                         10000, 10000, 0,  2*PI/10000, 0.5 , 1,
                           //x, y,vx,vy, L
                             0, 0, 0, 0, 0.00);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r, E, index
                             1, 1, 0, 2 * PI, 0.5, 1,
                           //x, y,vx,vy, L
                             1.0, 0, 0, 0.02, 0);
            cells = new Cell[] {a, b};
            bonds = new Bond[] {harmonicAverageBond(a, b, 0.0)};
            dt = Double.parseDouble(args[1]);
            int nSteps = Integer.parseInt(args[2]);
            testSimulation(cells, bonds, nSteps);
        } else if("orbit".equals(args[0])) {
            Cell a, b;
            cellStatusReading = (c) -> {
                return String.format(" %f %f ", c.x, c.y);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r   , E
                         10000, 1, 0,  2*PI, 0.5 , 0,
                           //x, y,th,vx,vy, L
                             0, 0, 0, 0, 0);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r, E
                             1, 1, 0, 2 * PI, 0.5, 0,
                           //x, y,th,vx,vy, L
                             1.0, 0, 0, .1, 0);
            cells = new Cell[] {a, b};
            bonds = new Bond[] {harmonicAverageBond(a, b, 0.0)};
            a.right = b; b.left = a;
            dt = Double.parseDouble(args[1]);
            int nSteps = Integer.parseInt(args[2]);
            testSimulation(cells, bonds, nSteps);
        } else if("relative motion".equals(args[0])) {
            Cell a, b;
            a = new Cell(null, null, null, null,
                             1, 1, 0, 2*PI, 0.5, 1,
                             -0.7, 0, 1, 0, 0);
            b = new Cell(null, null, null, null,
                         1, 1, 0, 2*PI, 0.5, 1,
                         +0.7, 0, 1, 0, 0);
            cells = new Cell[] {a, b};
            bonds = new Bond[] {harmonicAverageBond(a, b, 0.0)};
            a.right = b; b.left = a;
            dt = Double.parseDouble(args[1]);
            int nSteps = Integer.parseInt(args[2]);
            testSimulation(cells, bonds, nSteps);
        } else if("basic rotation".equals(args[0])) {
            Cell a, b;
            cellStatusReading = (c) -> {
                return String.format(" %f %f ", c.x, c.y);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r  , E, index
                             1, 1, 1, 2*PI, 1, 0,
                           //x, y,th,vx,vy, L
                            -2, 0, 0, -2, 0);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r  , E, index
                             1, 1, 1, 2*PI, 1, 0,
                           //x, y,th,vx,vy, L
                             2, 0, 0, 2, 0);
            cells = new Cell[] {a, b};
            bonds = new Bond[] {harmonicAverageBond(a, b, 0.0)};
            a.right = b; b.left = a;
            dt = Double.parseDouble(args[1]);
            int nSteps = Integer.parseInt(args[2]);
            testSimulation(cells, bonds, nSteps);
        } else if("relative motion damping".equals(args[0])) {
            Cell a, b;
            a = new Cell(null, null, null, null,
                             1, 1, 0.1, 1, 1, 1,
                             -2, 0, 1, 0, 0);
            b = new Cell(null, null, null, null,
                         1, 1, 0.1, 1, 1, 1,
                         2, 0, 1, 0, 0);
            cells = new Cell[] {a, b};
            bonds = new Bond[] {harmonicAverageBond(a, b, 0.0)};
            a.right = b; b.left = a;
            dt = Double.parseDouble(args[1]);
            int nSteps = Integer.parseInt(args[2]);
            testSimulation(cells, bonds, nSteps);
        } else if("minimize".equals(args[0])) {
            Cell a, b;
            a = new Cell(null, null, null, null, 10000, 10000, 0, 1, 0.5 , 1, 0, 0, 0, 0, -0.0 * 2 * PI);
            b = new Cell(null, null, null, null, 1, 1, 0, 1, 0.5, 1, 0, -1, 0.3, 0.0, 0.0 * 2*PI);
            cells = new Cell[] {a, b};
            bonds = new Bond[] {harmonicAverageBond(a, b, 0.0)};
            dt = 0.001;
            int nSteps = 100;
            double[] energies = new double[nSteps+1];
            ka = Double.parseDouble(args[1]);
            kb = Double.parseDouble(args[2]);
            kc = Double.parseDouble(args[3]);
            kd = Double.parseDouble(args[4]);
            bodyStatusReading = (body) -> {
                energies[stepCounter] = body.totalEnergy();
                return "" + energies[stepCounter];
            };
            cellStatusReading = (c) -> { return "";};
            testSimulation(cells, bonds, nSteps);
            System.err.println("Run finished with " + Arrays.toString(energies));
            System.out.println(util.Math.standardDeviation(energies));
        }
    }
    static void testSimulation(Cell[] cells, Bond[] bonds, int nSteps) {
        SoftBody bod = new SoftBody(cells, bonds);
        double t = 0;
        for(int i = 0; i < nSteps; i++) {
            System.out.format("%f ", t);
            bod.updateForces();
            bod.propagateCells();
            t += dt;
        }
        System.err.println("Number of Calls " + counter);
    }
}
