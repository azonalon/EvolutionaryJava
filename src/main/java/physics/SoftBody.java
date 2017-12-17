package physics;
import physics.*;
import java.util.ArrayList;
import static util.Math.*;
import java.util.function.*;
import static java.lang.Math.PI;
// TODO: switch to jalgo for better performance
// import org.ojalgo.matrix.store.*;

import static org.ejml.sparse.csc.CommonOps_DSCC.*;
import static org.ejml.dense.row.CommonOps_DDRM.*;
import org.ejml.data.*;

import static physics.Cell.dt;
import static physics.Cell.t;

public class SoftBody {
    public DMatrixRMaj X0; // newly calculated state
    public DMatrixRMaj X1; // last state
    public DMatrixRMaj X2; // oldest state
    public DMatrixRMaj fInt, v, vSquared, fExt, f, fDamp; //for debug
    public DMatrixRMaj M, MI;
    public DMatrixRMaj temp1, temp2;
    public DMatrixRMaj D1, D2;
    int m; // number of cells
    int s; // number of bonds
    final double t0 = 0.2;


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
        fDamp = new DMatrixRMaj(m * 3, 1);
        v = new DMatrixRMaj(m * 3, 1);
        vSquared = new DMatrixRMaj(m * 3, 1);
        M = new DMatrixRMaj(m * 3, 1);
        MI = new DMatrixRMaj(m * 3, 1);
        temp1 = new DMatrixRMaj(m * 3, 1);
        temp2 = new DMatrixRMaj(m * 3, 1);
        D1 = new DMatrixRMaj(m*3, m*3);
        D2 = new DMatrixRMaj(m*3, m*3);
        // D = new DMatrixSparseCSC(m *3, m*3, 114 * s);
        for(int i = 0; i < m; i++) {
            Cell c = cells[i];

            cells[i].linkToBody(this);
            cells[i].index = 3*i;

            X1.set(3 * i, 0, c.x);
            X1.set(3 * i + 1, 0, c.y);
            X1.set(3 * i + 2, 0, c.theta);

            X0.set(3 * i, 0, c.x + c.vx*dt);
            X0.set(3 * i + 1, 0, c.y + c.vy*dt);
            X0.set(3 * i + 2, 0, c.theta + c.L*dt);

            M.set(3*i, 0, c.m);
            M.set(3*i + 1, 0, c.m);
            M.set(3*i + 2, 0, c.I);
            MI.set(3*i, 0, 1/c.m);
            MI.set(3*i + 1, 0, 1/c.m);
            MI.set(3*i + 2, 0, 1/c.I);
        }
    }
    // DMatrixSparseCSC Hessian(DMatrixRMaj X) {
    //     return D;
    // }
    void updateHessian(int i, int j,
                        double sx, double sy, double k, double E,
                        double fS, double kdl) {
        int x1 = i + 0, y1 = i + 1, t1 = i + 2;
        int x2 = j + 0, y2 = j + 1, t2 = j + 2;
        double Dx1x1 = k + sy*sy*(6*E - kdl) - 2*fS*sx*sy;
        double Dy1y1 = k + sx*sx*(6*E - kdl) + 2*fS*sx*sy;
        double Dx1y1 =     sx*sy*(6*E + kdl) + 1*fS*(sy*sy-sx*sx);
        double Dx1t1 = -3*E*sy;
        double Dy1t1 = -3*E*sx;
        double Dt1t1 = 2*E;
        double Dt1t2 = E;
        if(t >= t0 && t < t0 + dt && s < 5) {
            System.out.format("UHessian: sx=%f, sy=%f, k=%f, E=%f, kdl=%f, fS=%f\n"
                            , sx, sy, k, E, kdl, fS);
        }

        addSymmetric(x1, x1,  Dx1x1, D1);
        addSymmetric(x1, y1,  Dx1y1, D1);
        addSymmetric(x1, t1,  Dx1t1, D1);
        addSymmetric(y1, y1,  Dy1y1, D1);
        addSymmetric(y1, t1,  Dy1t1, D1);
        addSymmetric(t1, t1,  Dt1t1, D1);
        addSymmetric(t1, t2,  Dt1t2, D1);

        addSymmetric(x2, x2, +Dx1x1, D1);
        addSymmetric(y2, y2, +Dy1y1, D1);
        addSymmetric(x2, y2, +Dx1y1, D1);
        addSymmetric(x2, t2, -Dx1t1, D1);
        addSymmetric(y2, t2, -Dy1t1, D1);
        addSymmetric(t2, t2, +Dt1t1, D1);

        addSymmetric(x1, x2, -Dx1x1, D1);
        addSymmetric(x1, y2, -Dx1y1, D1);
        addSymmetric(x1, t2, +Dx1t1, D1);
        addSymmetric(y1, x2, -Dx1y1, D1);
        addSymmetric(y1, y2, -Dy1y1, D1);
        addSymmetric(y1, t2, +Dy1t1, D1);
        addSymmetric(t1, x2, -Dx1t1, D1);
        addSymmetric(t1, y2, -Dy1t1, D1);

    }
    /**
     * M_ij -> M_ij + X
     * M_ji -> M_ij + X
     */
    void addSymmetric(int i, int j, double x, DMatrixRMaj M) {
        double x0 = M.get(i, j);
        M.set(i, j, x0 + x);
        if(i != j)
            M.set(j, i, x0 + x);
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

            double ex = dx / d;
            double ey = dy / d;
            double sx = ex / d;
            double sy = ey / d;
            double phi = circleMod(Math.atan2(ex, ey) + bond.angle - 1 * PI/2);
            if(bond.phi0 * phi < -PI) {
                if(bond.phi0 > phi) {
                    bond.rotationCounter++;
                } else {
                    bond.rotationCounter--;
                }
            }
            bond.phi0 = phi;
            phi += bond.rotationCounter * 2 * PI;

            double fS = 3.0*E * (th1 + th2 + 2 * phi);


            assert d>0.01: "Cell distance too small: d="  + d;

            double fx1 = (1 - l/d) * k * dx + fS * sy;
            double fy1 = (1 - l/d) * k * dy - fS * sx;
            fInt.set(a+0, 0, fInt.get(a+0, 0) + fx1);
            fInt.set(a+1, 0, fInt.get(a+1, 0) + fy1);
            fInt.set(b+0, 0, fInt.get(b+0, 0) - fx1);
            fInt.set(b+1, 0, fInt.get(b+1, 0) - fy1);

            fInt.set(a+2, 0, fInt.get(a+2, 0) -E * (2 * th1 + th2 + 3 * phi));
            fInt.set(b+2, 0, fInt.get(b+2, 0) -E * (2 * th2 + th1 + 3 * phi));

            double kdl = k * d * l;
            updateHessian(a, b, sx, sy, k, E, fS, kdl);
            energy += 0.5 * (l-d)*(l-d) * k + E * (sqrd(th1 + th2 + 2 * phi) -
            (th1 + phi) * (th2 + phi));
        }
        if(t >= t0 && t <= t0 + dt && s < 5) {
            // System.out.println("Bond");
            System.out.println(bonds[0].toString());
            System.out.println("X1");
            System.out.println(X1.toString());
            // System.out.println("v");
            // System.out.println(v.toString());
            System.out.println("Hessian");
            System.out.println(D1.toString());
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
        D2.set(D1);
        D1.zero();
        fDamp.zero();
        // zero(f, 0, m*3, 0, 1);
        add(+1.0/dt, X1,-1.0/dt, X2, v);
        add(innerForce(X1, X2), externalForce(), f);
        mult(D1, v, fDamp);
        // add(-1.0, fDamp, f, f);
        elementMult(MI, f, f);

        add(2, X1, -1, X2, X0);
        add(dt*dt, f, 1, X0, X0);

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
