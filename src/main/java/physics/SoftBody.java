package physics;
import physics.*;
import java.util.*;
import static util.Math.*;
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
    static boolean startProgram = true;

    public static Consumer<Cell> cellStatusCallback = (Cell c)->{};
    public static Consumer<SoftBody> bodyStatusCallback = (SoftBody b)->{};
    public static Consumer<Bond> bondStatusCallback = (Bond b)->{};

    /**
     * checks if the bond forces between the cells have been added to each cell.
     * this is a symmetric operation a <--> b
     * @param  Cell a             [description]
     * @param  Cell b             [description]
     * @return      [description]
     */
    public void propagateCells() {
        for(Cell c: cells) {
            cellStatusCallback.accept(c);
            energy += c.L * c.L * c.I/2.0 + c.vx * c.vx * c.m /2.0 + c.vy * c.vy * c.m /2.0;
            c.propagate();
        }
        bodyStatusCallback.accept(this);
    }

    /**
     * Initializes forces to zero before each step.
     * Goes through all the bonds and adds interaction forces to both cells
     * in the bond
     */
    public void updateForces() {
        energy=0;
        for(Cell c: cells) {
            c.fX=0;
            c.fY=0;
            c.T=0;
        }
        for(Bond bond: bonds) {
            addBondForces(bond);
            bondStatusCallback.accept(bond);
        }
    }

    /**
     * Calculates and add forces for the two cells in  a bond.
     * @param Bond b
     */
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
        double fShear = + 6 * l * l * E * (first.theta + second.theta - 2 * phi);
        assert(d>0.01);


        first.fX  += (d - l) * k * dx - fShear * (-1*dy);
        first.fY  += (d - l) * k * dy - fShear * (   dx);
        second.fX += -(d - l) * k * dx + fShear * (-1*dy);
        second.fY += -(d - l) * k * dy + fShear * (   dx);

        // first.T  -= 2 * E * l * l * l * (2 * first.theta + second.theta - 3 * phi);
        // second.T -= 2 * E * l * l * l * (2 * second.theta + first.theta - 3 * phi);
        first.T  += 6 * E * l * l * l * (-second.theta - first.theta  + 2*phi - 3*(first.theta - second.theta));
        second.T += 6 * E * l * l * l * (+first.theta  + second.theta - 2*phi - 3*(second.theta - first.theta));


        energy += (d-l)*(d-l)* k/2.0;
        energy -= 2 * l * l * l * E * (
        sqrd(first.theta + second.theta - 2 * phi) -
        (first.theta - phi) * (second.theta - phi) * 1
        );
        kb = 0;
        kc = 0.0;
        // energy += ka * sqrd(first.theta + - phi) + kb;
        // energy += kb * E * l * l * l * (sqrd(2 * second.theta + first.theta - 3 * phi) -
        //                 2*(2*second.theta - 2*phi) * (first.theta - phi));
        // energy += kc * E * l * l * l * sqrd(2 * first.theta + second.theta - 3 * phi);

        // if(startProgram == true) {
            // System.err.format("Parameters: k=%f, E=%f, D=%f, c=%f, l=%f, m1=%f, m2=%f\n",
            //                     k, E, D, c, l, first.m, second.m);
            // startProgram = false;
        // }


        // damping forces
        double vXRelative = first.vx - second.vx;
        double vYRelative = first.vy - second.vy;
        // subtract parts orthogonal to direction vector
        double vOrthogonal = vXRelative * -1 * dy + vYRelative * dx;
        // double dphi = Math.asin(vOrthogonal * dt/d);
        // double shearDamping = - D * (first.L - dphicalc + second.L - dphicalc);
        // orthogonalFraction += shearDamping;
        double odf1 = 1;// - 0.1 * Math.abs(first.L - dphi);
        double odf2 = 1;// - 0.1 * Math.abs(second.L - dphi);
        vXRelative -=  vOrthogonal * -1 * dy;
        vYRelative -=  vOrthogonal * dx;
        // System.err.format("orth vel.: %f, orth Damping: %f, Torque: %f, phi %f\n",orthogonalFraction, shearDamping, second.T, phi);
        // System.err.format("fShear: %f, d: %f, dphi: %f, vxR: %f, vyR: %f\n", fShear, d, dphicalc, vXRelative, vYRelative);
        // System.err.format("dx: %f, dy: %f, phi: %f\n", dx, dy, phi);
        // System.err.format("L1: % 04.8f, L2: % 04.8f\n", first.L, second.L);
        // System.err.format("th1: %f, th2: %f, phi: %f, odf1: %f, odf2: %f\n", first.theta, second.theta, phi, odf1, odf2);

        first.fX  -= c * vXRelative * odf2;
        first.fY  -= c * vYRelative * odf2;
        second.fX += c * vXRelative * odf1;
        second.fY += c * vYRelative * odf1;
        first.T  -= D * (first.L - second.L);
        second.T += D * (first.L - second.L);
    }

    public double totalEnergy() {
        return energy;
    }
}
