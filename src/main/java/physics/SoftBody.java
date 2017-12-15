package physics;
import physics.*;
import java.util.ArrayList;
import static util.Math.*;
import java.util.function.*;
import static java.lang.Math.PI;
import org.ejml.data.*;
// import static org.ejml.sparse.csc.CommonOps_DSCC.*;
// import static org.ejml.simple.ops.SimpleOperations_DDRM.*;
import static org.ejml.dense.row.CommonOps_DDRM.*;
// import static java.lang.Math.sqrt;
// import static System.out.format;
import static physics.Cell.dt;

public class SoftBody {
    public DMatrixRMaj X0; // newly calculated state
    public DMatrixRMaj X1; // last state
    public DMatrixRMaj X2; // oldest state
    public DMatrixRMaj fInt, v, vSquared, fExt, f; //for debug
    public DMatrixRMaj M, MI;
    public DMatrixRMaj temp1, temp2;
    public DMatrixSparceCSC D;
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
        X0 = new DMatrixRMaj(m * 3, 1);
        X1 = new DMatrixRMaj(m * 3, 1);
        X2 = new DMatrixRMaj(m * 3, 1);
        fExt = new DMatrixRMaj(m * 3, 1);
        fInt = new DMatrixRMaj(m * 3, 1);
        f = new DMatrixRMaj(m * 3, 1);
        v = new DMatrixRMaj(m * 3, 1);
        vSquared = new DMatrixRMaj(m * 3, 1);
        M = new DMatrixRMaj(m * 3, 1);
        MI = new DMatrixRMaj(m * 3, 1);
        temp1 = new DMatrixRMaj(m * 3, 1);
        temp2 = new DMatrixRMaj(m * 3, 1);
        for(int i = 0; i < m; i++) {
            Cell c = cells[i];

            cells[i].linkToBody(this);
            cells[i].index = 3*i;
            X0.set(3 * i, 0, c.x + c.vx*dt);
            X0.set(3 * i + 1, 0, c.y + c.vy*dt);
            X0.set(3 * i + 2, 0, c.theta + c.L*dt);

            X1.set(3 * i, 0, c.x);
            X1.set(3 * i + 1, 0, c.y);
            X1.set(3 * i + 2, 0, c.theta);

            M.set(3*i, 0, c.m);
            M.set(3*i + 1, 0, c.m);
            M.set(3*i + 2, 0, c.I);
            MI.set(3*i, 0, 1/c.m);
            MI.set(3*i + 1, 0, 1/c.m);
            MI.set(3*i + 2, 0, 1/c.I);
        }
    }
    DMatrixSparseCSC Hessian(DMatrixRMaj X) {
        return D;
    }
    DMatrixRMaj innerForce(DMatrixRMaj X1, DMatrixRMaj X2) {
        fInt.zero();
        for(Bond bond: bonds) {
            int a = bond.a.index;
            int b = bond.b.index;
            double dx  = X1.get(b + 0, 0) - X1.get(a + 0, 0);
            double dy  = X1.get(b + 1, 0) - X1.get(a + 1, 0);
            double th1 = X1.get(a + 2, 0);
            double th2 = X1.get(b + 2, 0);
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

            fInt.set(a+0, 0, fInt.get(a+0, 0)  +  (d - l) * k * dx - fShear * (-1*dy));
            fInt.set(a+1, 0, fInt.get(a+1, 0)  +  (d - l) * k * dy - fShear * (   dx));
            fInt.set(b+0, 0, fInt.get(b+0, 0)  + -(d - l) * k * dx + fShear * (-1*dy));
            fInt.set(b+1, 0, fInt.get(b+1, 0)  + -(d - l) * k * dy + fShear * (   dx));

            fInt.set(a+2, 0, fInt.get(a+2, 0) +  -E * (2 * th1 + th2 + 3 * phi));
            fInt.set(b+2, 0, fInt.get(b+2, 0) +  -E * (2 * th2 + th1 + 3 * phi));

            energy += 0.5 * (l-d)*(l-d) * k + E * (sqrd(th1 + th2 + 2 * phi) -
            (th1 + phi) * (th2 + phi)  );
        }
        return fInt;
    }

    DMatrixRMaj externalForce() {
        fExt.zero();
        for(Cell c: cells) {
           fExt.set(c.index, 0, fExt.get(c.index, 0) +  c.fX);
           fExt.set(c.index + 1, 0, fExt.get(c.index + 1, 0) + c.fY);
           fExt.set(c.index + 2, 0, fExt.get(c.index + 2, 0) + c.T);
            c.fX = 0;
            c.fY = 0;
            c.T = 0;
        }
        return fExt;
    }

    void explicitEulerStep() {
        energy = 0;
        X2.set(X1);
        X1.set(X0);
        // zero(f, 0, m*3, 0, 1);
        add(innerForce(X1, X2), externalForce(), f);
        elementMult(MI, f, f);
        add(2, X1, -1, X2, X0);
        add(dt*dt, f, 1, X0, X0);
        // X0 = X1.multiply(2.0).add(
        //     MI.hadamardProduct(f).multiply(dt*dt).
        //     subtract(X2)
        // );
        add(-1/dt, X2, 1/dt, X1, v);
        elementMult(v, v, vSquared);
        energy += 0.5 * dot(M, vSquared);
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
