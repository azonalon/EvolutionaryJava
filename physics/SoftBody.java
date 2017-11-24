package physics;
import physics.*;
import java.util.*;
// import util.Math;
import java.util.function.*;
import static java.lang.Math.PI;
import static java.lang.Math.sqrt;
// import static System.out.format;

public class SoftBody {

    public SoftBody(Cell[] cells, Bond[] bonds) {
        this.cells = cells;
        this.bonds = bonds;
    }


    Cell[] cells;
    Bond[] bonds;
    double energy = 0;
    public static double ka=0, kb=-0.0, kc=0, kd=0;
    /*
     * STATIC methods
     */
    public static double sqrd(double a) {
        return a*a;
    }


    public static double phi0;
    static int phiCounter = 0;
    static int LEFT = 1 << 0;
    static int RIGHT= 1 << 1;
    static int UP   = 1 << 2;
    static int DOWN = 1 << 3;
    public static double dt = 0.1;
    static boolean startProgram = true;
    public static Function<Cell, String> cellStatusReading =
        (Cell c)->{return String.format("%f ", c.x);};

    public static Function<SoftBody, String> bodyStatusReading =
        (body)->{return body.totalEnergy() + "";};


    /**
     * checks if the bond forces between the cells have been added to each cell.
     * this is a symmetric operation a <--> b
     * @param  Cell a             [description]
     * @param  Cell b             [description]
     * @return      [description]
     */

    public void propagateCells() {
        for(Cell c: cells) {
            System.out.print(cellStatusReading.apply(c) + " ");
            energy += c.L * c.L * c.I/2.0 + c.vx * c.vx * c.m /2.0 + c.vy * c.vy * c.m /2.0;
            c.propagate();
        }
        System.out.print(bodyStatusReading.apply(this));
        System.out.println();
    }

    public void updateForces() {
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
        phi0 = phi;
        double fShear = - 6 * l * l * E * (first.theta + second.theta - 2 * phi);
        assert(d>0.01);


        first.fX  += (d - l) * k * dx + fShear * (-1*dy);
        first.fY  += (d - l) * k * dy + fShear * (   dx);
        second.fX -= (d - l) * k * dx + fShear * (-1*dy);
        second.fY -= (d - l) * k * dy + fShear * (   dx);

        first.T  -= 2 * E * l * l * l * (2 * first.theta + second.theta - 3 * phi);
        second.T -= 2 * E * l * l * l * (2 * second.theta + first.theta - 3 * phi);


        energy += (d-l)*(d-l)* k/2.0;
        energy -= 2 * l * l * l * E * (
        sqrd(first.theta + second.theta - 2 * phi) -
        (first.theta - phi) * (second.theta - phi)
        );

        energy += ka * sqrd(first.theta + - phi) + kb;
        // energy += kb * E * l * l * l * sqrd(2 * second.theta + first.theta - 3 * phi);
        // energy += kc * E * l * l * l * sqrd(2 * second.theta + first.theta - 3 * phi);

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
        double shearDamping = 0.00 * D * (first.L - dphicalc + second.L - dphicalc);
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
}
