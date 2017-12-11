package physics;
import physics.*;
import java.util.ArrayList;
import static util.Math.*;
import java.util.function.*;
import static java.lang.Math.PI;
import org.la4j.*;
import org.la4j.vector.*;
import org.la4j.matrix.*;
// import static java.lang.Math.sqrt;
// import static System.out.format;
import static physics.Cell.dt;

public class SoftBody {
    public Vector X0; // newly calculated state
    public Vector X1; // last state
    public Vector X2; // oldest state
    public Vector f, v; //for debug
    Vector M, MI;
    int m; // number of cells
    int s; // number of bonds


    public SoftBody(Cell[] cells, Bond[] bonds) {
        initializeState(cells, bonds);
    }

    public void initializeState(Cell[] cells, Bond[] bonds) {
        this.cells = cells;
        this.bonds = bonds;
        this.m = cells.length;
        this.s = bonds.length;
        X0 = DenseVector.zero(m * 3);
        X1 = DenseVector.zero(m * 3);
        X2 = DenseVector.zero(m * 3);
        M = DenseVector.zero(m * 3);
        MI = DenseVector.zero(m * 3);
        for(int i = 0; i < m; i++) {
            Cell c = cells[i];

            cells[i].linkToBody(this);
            cells[i].index = 3*i;
            X0.set(3 * i, c.x + c.vx*dt);
            X0.set(3 * i + 1, c.y + c.vy*dt);
            X0.set(3 * i + 2, c.theta + c.L*dt);

            X1.set(3 * i, c.x);
            X1.set(3 * i + 1, c.y);
            X1.set(3 * i + 2, c.theta);

            M.set(3*i, c.m);
            M.set(3*i + 1, c.m);
            M.set(3*i + 2, c.I);
            MI.set(3*i, 1/c.m);
            MI.set(3*i + 1, 1/c.m);
            MI.set(3*i + 2, 1/c.I);
        }
    }
    Vector innerForce(Vector X1, Vector X2) {
        Vector f = DenseVector.zero(m * 3);
        for(Bond bond: bonds) {
            int a = bond.a.index;
            int b = bond.b.index;
            double dx  = X1.get(b + 0) - X1.get(a + 0);
            double dy  = X1.get(b + 1) - X1.get(a + 1);
            double th1 = X1.get(a + 2);
            double th2 = X1.get(b + 2);
            double d = Math.sqrt(dx*dx + dy*dy);

            double l=bond.l, E=bond.E, k=bond.k;//, c=bond.c, D=bond.D;

            dx = dx / d;
            dy = dy / d;
            double phi = circleMod(Math.atan2(dx, dy) + bond.angle - PI/2);
            if(bond.phi0 * phi < -PI) {
                if(bond.phi0 > phi) {
                    bond.rotationCounter++;
                } else {
                    bond.rotationCounter--;
                }
            }
            bond.phi0 = phi;
            phi += bond.rotationCounter * 2 * PI;
            double fShear = 3.0*E/d * (th1 + th2 + 2 * phi);
            assert d>0.01: "Cell distance too small: d="  + d;

            f.set(a + 0, (d - l) * k * dx - fShear * (-1*dy));
            f.set(a + 1, (d - l) * k * dy - fShear * (   dx));
            f.set(b + 0, -(d - l) * k * dx + fShear * (-1*dy));
            f.set(b + 1, -(d - l) * k * dy + fShear * (   dx));

            f.set(a + 2, -E * (2 * th1 + th2 + 3 * phi));
            f.set(b + 2, -E * (2 * th2 + th1 + 3 * phi));

            energy += 0.5 * (l-d)*(l-d) * k + E * (sqrd(th1 + th2 + 2 * phi) -
            (th1 + phi) * (th2 + phi)  );
        }
        return f;
    }

    Vector externalForce() {
        Vector f = DenseVector.zero(m * 3);
        for(Cell c: cells) {
            f.set(c.index, c.fX);
            f.set(c.index + 1, c.fY);
            f.set(c.index + 2, c.T);
            c.fX = 0;
            c.fY = 0;
            c.T = 0;
        }
        return f;
    }

    void explicitEulerStep() {
        energy = 0;
        X2 = X1.copy();
        X1 = X0.copy();
        f = innerForce(X1, X2).add(externalForce());
        X0 = X1.multiply(2.0).add(
            MI.hadamardProduct(f).multiply(dt*dt).
            subtract(X2)
        );
        v = (X1.subtract(X2)).multiply(1/dt);
        energy += 0.5 * v.hadamardProduct(M).innerProduct(v);
        for(Cell c: cells) {
            cellStatusCallback.accept(c);
            cellForceCallback.accept(c);
        }
        bodyStatusCallback.accept(this);
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
        ArrayList<Bond> bondList = new ArrayList<Bond>();
        ArrayList<Cell> cellList = new ArrayList<Cell>();
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
                    upper.x = x;
                    upper.y = y;
                    cellList.add(upper);
                } else {
                    x += cellWidth;
                    continue;
                }
                if(right != null) {
                    bondList.add(Bond.harmonicAverageBond(upper, right, 0));
                }
                if(lower != null) {
                    bondList.add(Bond.harmonicAverageBond(upper, lower, -PI/2));
                }
                x += cellWidth;
            }
            y -= cellHeight;
        }
        initializeState(
            cellList.toArray(new Cell[cellList.size()]),
            bondList.toArray(new Bond[bondList.size()])
        );
    }

    Cell[] cells;
    Bond[] bonds;
    double energy = 0;
    static boolean startProgram = true;

    public static Consumer<Cell> cellStatusCallback = (Cell c)->{};
    public Consumer<Cell> cellForceCallback  = (Cell c)->{};
    public static Consumer<SoftBody> bodyStatusCallback = (SoftBody b)->{};
    public static Consumer<Bond> bondStatusCallback = (Bond b)->{};



    public double totalEnergy() {
        return energy;
    }
}
