package physics;
import physics.*;
import java.util.*;
import static util.Math.*;
import java.util.function.*;
import static java.lang.Math.PI;
// import static java.lang.Math.sqrt;
// import static System.out.format;

public class SoftBody {
    public SoftBody(Cell[] cells, Bond[] bonds) {
        this.cells = cells;
        this.bonds = bonds;
    }
    /**
     * Creates a SoftBody from an array of predefined cell
     * types. Goes through the whole array and makes bonds
     * where non null cells are adjacent. You have to make
     * the cells in the array connected, otherwise isolated
     * indices will be loose.
     * @param  CellCulture[][] types  grid of celltype objects
     * @return SoftBody according to celltype grid
     */
    public SoftBody(Cell[][] cells, double cellWidth, double cellHeight) {
        Vector<Bond> bondVector = new Vector<Bond>();
        Vector<Cell> cellVector = new Vector<Cell>();
        Cell upper, lower, right;
        double x=0, y=0;
        y = cellHeight * (cells.length - 1);
        for(int row = 0; row<cells.length; row++) {
            x = 0;
            // System.out.println("y" + y);
            for(int col = 0; col<cells[row].length; col++) {
                // System.out.println("x" + x);
                lower = null;
            right = null;
                if(row+1 < cells.length)
                    lower = cells[row+1][col  ];
                if(col+1 < cells[row].length)
                    right = cells[row  ][col+1];
                upper = cells[row][col];

                if(upper != null) {
                    upper.setPosition(x, y);
                    cellVector.add(upper);
                } else {
                    x += cellWidth;
                    continue;
                }
                if(right != null) {
                    bondVector.add(Bond.harmonicAverageBond(upper, right, 0));
                }
                if(lower != null) {
                    bondVector.add(Bond.harmonicAverageBond(upper, lower, -PI/2));
                }
                x += cellWidth;
            }
            y -= cellHeight;
        }
        this.bonds = bondVector.toArray(new Bond[bondVector.size()]);
        this.cells = cellVector.toArray(new Cell[cellVector.size()]);
    }

    Cell[] cells;
    Bond[] bonds;
    double energy = 0;
    public static double ka=0, kb=-0.0, kc=0, kd=0;
    static boolean startProgram = true;

    public static Consumer<Cell> cellStatusCallback = (Cell c)->{};
    public Consumer<Cell> cellForceCallback  = (Cell c)->{};
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
            energy += c.L * c.L * (c.I)/2.0 + c.vx * c.vx * c.m /2.0 + c.vy * c.vy * c.m /2.0;
            cellForceCallback.accept(c);
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
        double th1 = first.theta;
        double th2 = second.theta;
        double d = Math.sqrt(dx*dx + dy*dy);

        double l=b.l, E=b.E, k=b.k, c=b.c, D=b.D;

        dx = dx / d;
        dy = dy / d;
        double phi = circleMod(Math.atan2(dy, dx) - b.angle);
        // double phi = Math.atan2(dy, dx) - b.angle;
        // double phi = Math.atan2(dy, dx);
        // System.err.format("phi: %f, b.angle %f\n", phi, b.angle);
        if(b.phi0 * phi < -PI) {
            if(b.phi0 > phi) {
                b.rotationCounter++;
            } else {
                b.rotationCounter--;
            }
        }
        b.phi0 = phi;
        phi += b.rotationCounter * 2 * PI;
        double fShear = 6* l * l * E * (th1 + th2 - 2 * phi);
        assert d>0.01: "Cell distance too small: d="  + d;

        first.fX  += (d - l) * k * dx - fShear * (-1*dy);
        first.fY  += (d - l) * k * dy - fShear * (   dx);
        second.fX += -(d - l) * k * dx + fShear * (-1*dy);
        second.fY += -(d - l) * k * dy + fShear * (   dx);

        first.T  -= 2 * E * l * l * l * (2 * th1 + th2 - 3 * phi);
        second.T -= 2 * E * l * l * l * (2 * th2 + th1 - 3 * phi);
        // double k1=1, k2=1;
        // first.T  += E * l * l * l * ((-k1*second.theta - k2 *first.theta  + (k1+k2)*phi) - 1*(first.theta - second.theta));
        // second.T += E * l * l * l * ((+k1*first.theta  + k2 * second.theta - (k1+k2)*phi) - 1*(second.theta - first.theta));


        energy += (d-l)*(d-l)* k/2.0;
        energy += 2*l * l * l * E * (
        sqrd(th1 + th2 - 2 * phi) -
        (th1 - phi) * (th2 - phi) * 1
        );

        // damping forces
        double vXRelative = first.vx - second.vx;
        double vYRelative = first.vy - second.vy;
        // subtract parts orthogonal to direction vector
        double vOrthogonal = vXRelative * -1 * dy + vYRelative * dx;
        double odf1 = 1;// - 0.1 * Math.abs(first.L - dphi);
        double odf2 = 1;// - 0.1 * Math.abs(second.L - dphi);
        vXRelative -=  vOrthogonal * -1 * dy;
        vYRelative -=  vOrthogonal * dx;
        // System.err.format("orth vel.: %f, orth Damping: %f, Torque: %f, phi %f\n",orthogonalFraction, shearDamping, second.T, phi);
        // System.err.format("fShear: %f, d: %f, dphi: %f, vxR: %f, vyR: %f\n", fShear, d, dphicalc, vXRelative, vYRelative);
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
