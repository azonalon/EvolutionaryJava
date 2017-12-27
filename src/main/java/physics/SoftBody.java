package physics;
import java.nio.file.*;
import java.io.*;
import physics.*;
import java.util.ArrayList;
import static util.Math.*;
import java.util.function.*;
import static java.lang.Math.PI;
// TODO: switch to jalgo for better performance
// import org.ojalgo.matrix.store.*;

import static org.ejml.sparse.csc.CommonOps_DSCC.*;
import static org.ejml.dense.row.CommonOps_DDRM.*;
import static org.ejml.dense.row.SpecializedOps_DDRM.*;
import org.ejml.data.*;

import static physics.Cell.dt;
import static physics.Cell.t;

public class SoftBody {
    public DMatrixRMaj X0; // newly calculated state
    public DMatrixRMaj X1; // last state
    public DMatrixRMaj X2; // oldest state
    public DMatrixRMaj XHat, H;
    public DMatrixRMaj fInt, oldF, V0, vSquared, fExt, f, fDamp, dX, G;
    public DMatrixRMaj M, MI, DC; // mass, inverse mass, damping coefficient
    public DMatrixRMaj temp1, temp2;
    public DMatrixRMaj D1, D2;
    public double phi0, dphi0;
    public int i;
    int m; // number of cells
    int s; // number of bonds
    final double t0 = 2.0;
    Cell[] cells;
    Bond[] bonds;
    double energy;
    double kdamp = 1;
    double eInt; // inner energy
    static boolean startProgram = true;

    public static Consumer<Cell> cellStatusCallback = (Cell c)->{};
    public Consumer<Cell> cellForceCallback  = (Cell c)->{};
    public static Consumer<SoftBody> bodyStatusCallback = (SoftBody b)->{};
    public static Consumer<Bond> bondStatusCallback = (Bond b)->{};



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
        XHat = new DMatrixRMaj(m * 3, 1);
        G = new DMatrixRMaj(m * 3, 1);
        fExt = new DMatrixRMaj(m * 3, 1);
        fInt = new DMatrixRMaj(m * 3, 1);
        oldF = new DMatrixRMaj(m * 3, 1);
        dX = new DMatrixRMaj(m * 3, 1);
        f = new DMatrixRMaj(m * 3, 1);
        fDamp = new DMatrixRMaj(m * 3, 1);
        V0 = new DMatrixRMaj(m * 3, 1);
        vSquared = new DMatrixRMaj(m * 3, 1);
        M = new DMatrixRMaj(m * 3, m * 3);
        DC = new DMatrixRMaj(m * 3, m * 3);
        MI = new DMatrixRMaj(m * 3, m*3);
        H = new DMatrixRMaj(m * 3, m*3);
        temp1 = new DMatrixRMaj(m * 3, 1);
        temp2 = new DMatrixRMaj(m * 3, m * 3);
        D1 = new DMatrixRMaj(m*3, m*3);
        D2 = new DMatrixRMaj(m*3, m*3);
        // D = new DMatrixSparseCSC(m *3, m*3, 114 * s);
        for(int i = 0; i < m; i++) {
            Cell c = cells[i];

            cells[i].linkToBody(this);
            cells[i].index = 3*i;

            X0.set(3 * i, 0, c.x);
            X0.set(3 * i + 1, 0, c.y);
            X0.set(3 * i + 2, 0, c.theta);

            X1.set(3 * i, 0, c.x - c.vx*dt);
            X1.set(3 * i + 1, 0, c.y - c.vy*dt);
            X1.set(3 * i + 2, 0, c.theta - c.L*dt);

            // diagonal matrices
            M.set(3*i,     3*i + 0, c.m);
            M.set(3*i + 1, 3*i + 1, c.m);
            M.set(3*i + 2, 3*i + 2, c.I);
            DC.set(3*i,     3*i + 0, kdamp * c.zeta);
            DC.set(3*i + 1, 3*i + 1, kdamp * c.zeta);
            DC.set(3*i + 2, 3*i + 2, kdamp * c.alpha);
            MI.set(3*i,     3*i+0, 1/c.m);
            MI.set(3*i + 1, 3*i+1, 1/c.m);
            MI.set(3*i + 2, 3*i+2, 1/c.I);
        }
    }
    // DMatrixSparseCSC Hessian(DMatrixRMaj X) {
    //     return D;
    // }
    void addEnergy(double k, double E, double d, double l, double th1, double th2, double phi) {
        eInt += 0.5 * (l-d)*(l-d) * k + E * (sqrd(th1 + th2 - 2 * phi) -
        (th1 - phi) * (th2 - phi));
    }

    void addForce(double k, double E, double d, double l,
                  double th1, double th2, double phi,
                  double sx, double sy, double fS, int a, int b) {
        double fx1 = d*(d - l) * k * sx + fS * sy;
        double fy1 = d*(d - l) * k * sy - fS * sx;
        fInt.set(a+0, 0, fInt.get(a+0, 0) + fx1);
        fInt.set(a+1, 0, fInt.get(a+1, 0) + fy1);
        fInt.set(b+0, 0, fInt.get(b+0, 0) - fx1);
        fInt.set(b+1, 0, fInt.get(b+1, 0) - fy1);

        fInt.set(a+2, 0, fInt.get(a+2, 0) -E * (2 * th1 + th2 - 3 * phi));
        fInt.set(b+2, 0, fInt.get(b+2, 0) -E * (2 * th2 + th1 - 3 * phi));
    }

    void addStress(double k, double E, double d, double l,double sx, double sy,
                    double fS, int i, int j) {
        int x1 = i + 0, y1 = i + 1, t1 = i + 2;
        int x2 = j + 0, y2 = j + 1, t2 = j + 2;
        double kdl = k*d*l;
        double Dx1x1 = k + sy*sy*(6*E - kdl) - 2*fS*sx*sy;
        double Dy1y1 = k + sx*sx*(6*E - kdl) + 2*fS*sx*sy;
        double Dx1y1 =   - sx*sy*(6*E - kdl) + 1*fS*(sx*sx-sy*sy);
        double Dx1t1 = -3*E*sy;
        double Dy1t1 =  3*E*sx;
        double Dt1t1 =  2*E;
        double Dt1t2 =  E;
        // if(t >= t0 && t < t0 + dt && s < 5) {
        //     System.out.format("UHessian: sx=%f, sy=%f, k=%f, E=%f, kdl=%f, fS=%f\n"
        //                     , sx, sy, k, E, kdl, fS);
        // }

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


    void updateEnergyForce() {
        fInt.zero();
        eInt = 0;
        for(Bond bond: bonds) {
            int a = bond.a.index;
            int b = bond.b.index;
            double dx  = X1.get(b + 0, 0) - X1.get(a + 0, 0);
            double dy  = X1.get(b + 1, 0) - X1.get(a + 1, 0);
            double th1 = X1.get(a + 2, 0);
            double th2 = X1.get(b + 2, 0);
            double d = Math.sqrt(dx*dx + dy*dy);
            double l=bond.l, E=bond.E, k=bond.k;//, c=bond.c, D=bond.D;

            double sx = dx /d/d;
            double sy = dy /d/d;

            double nu = circleMod(Math.atan2(dy, dx) - bond.angle) - 0 * PI/2;
            if(bond.nu0 * nu < -PI) {
                if(bond.nu0 > nu) {
                    bond.rotationCounter++;
                } else {
                    bond.rotationCounter--;
                }
            }
            bond.nu0 = nu;
            nu += bond.rotationCounter * 2 * PI;
            double fS = 3.0*E * (th1 + th2 - 2 * nu);

            addEnergy(k, E, d, l, th1, th2, nu);
            addForce(k, E, d, l, th1, th2, nu, sx, sy, fS, a, b);
        }
    }

    void updateStressForceEnergy() {
        D1.zero();
        fInt.zero();
        eInt = 0;
        for(Bond bond: bonds) {
            int a = bond.a.index;
            int b = bond.b.index;
            double dx  = X0.get(b + 0, 0) - X0.get(a + 0, 0);
            double dy  = X0.get(b + 1, 0) - X0.get(a + 1, 0);
            double th1 = X0.get(a + 2, 0);
            double th2 = X0.get(b + 2, 0);
            double d = Math.sqrt(dx*dx + dy*dy);

            double l=bond.l, E=bond.E, k=bond.k;//, c=bond.c, D=bond.D;

            double ex = dx / d;
            double ey = dy / d;
            double sx = ex / d;
            double sy = ey / d;
            double nu = circleMod(Math.atan2(dy, dx) - bond.angle) - 0 * PI/2;
            if(bond.nu0 * nu < -PI) {
                if(bond.nu0 > nu) {
                    bond.rotationCounter++;
                } else {
                    bond.rotationCounter--;
                }
            }
            bond.nu0 = nu;
            nu += bond.rotationCounter * 2 * PI;
            double fS = 3.0*E * (th1 + th2 - 2 * nu);
            // assert d>0.01: "Cell distance too small: d="  + d;
            addEnergy(k, E, d, l, th1, th2, nu);
            addForce(k, E, d, l, th1, th2, nu, sx, sy, fS, a, b);
            addStress(k, E, d, l,sx, sy, fS, a, b);
        }
        // if(t >= t0 && t <= t0 + dt && s < 5) {
            // System.out.println("Bond");
            // System.out.println(bonds[0].toString());
            // System.out.println("X1");
            // System.out.println(X1.toString());
            // System.out.println("v");
            // System.out.println(v.toString());
            // System.out.println("Hessian");
            // System.out.println(D1.toString());
        // }
    }

    DMatrixRMaj updateExternalForce() {
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

    double newtonStep(DMatrixRMaj X0, DMatrixRMaj XHat, DMatrixRMaj v,
                    DMatrixRMaj G, DMatrixRMaj H, DMatrixRMaj dX) {
        double dE;
        // mult(M, G, temp1);
        // negative sign because Nabla E = -f
        // caculate damping force
        // signs are positive because Nabla E = -f = -(-DC D2 v)
        dE = dot(G, G);
        add(1.0, X0, -1, G, X0);
        return dE;
    }

    double stableNewtonStep(DMatrixRMaj X0, DMatrixRMaj XHat, DMatrixRMaj V0,
                    DMatrixRMaj G, DMatrixRMaj H) {


        double k = 1e-2;
        double alpha = 1;
        double alphaMax  = 1e3;

        // compute the step
        updateStressForceEnergy();
        calculatePhiGradient();
        calculateHessian();

        invert(H, H);
        mult(-1, H, G, dX);
        // make sure dX is suitable
        double dE = dot(G, G);
        System.out.format("dE = %8.12f , phi0=%8.12f, eInt=%8.12f\n", dE, phi0, eInt);
        double dx = dot(dX, dX);
        dphi0 = dot(dX, G);
        double dXdE = dphi0 * dphi0 * Math.signum(dphi0);
        if(dXdE < -k*dx*dE) {
            // dx is suitable
            // System.out.format("dX is suitable! dXdE=%g, dX dE=%g\n", dXdE, dx*dE);
        }
        else if(dXdE > k*dx*dE) {
            // -dx is suitable
            System.out.format("-dX is suitable! dXdE=%g, dX dE=%g\n", dXdE, dx*dE);
            alpha = -alpha;
        }
        else {
            System.out.format("Gradient is suitable! dXdE=%g, dX dE=%g\n", dXdE, dx*dE);
            dX.set(G);
        }
        strongWolfeLineSearch(alpha, alphaMax);
        dE = dot(G, G);
        return dE;
    }

    void scanLineToFile(double alphaMax, int n, String filename) {
        try {
            Path path =  Paths.get("build/test-results/physics/SoftBody/ScanLine");
            if(Files.notExists(path)){
                    Files.createDirectories(path);
            }
            path = path.resolve(filename);
            BufferedWriter writer = Files.newBufferedWriter(path);
            writer.write(String.format("#alpha phi dphi\n"));
            lineSearchStep(-alphaMax);
            double step = alphaMax/(double)n;
            double phi1 = phi0;
            for(int j=-n; j<n; j++) {
                writer.write(String.format("%g %12.12f %12.12f %12.12f\n",
                                        j*step, phi0, dphi0, (phi0 - phi1)/step));
                phi1 = phi0;
                lineSearchStep(alphaMax/(double)n);
            }
            lineSearchStep(-alphaMax);
            writer.close();
        } catch(IOException e) {
            throw new RuntimeException("Could not write data");
        }
    }

    double strongWolfeLineSearch(double alpha1, double alphaMax) {
        // double c = 1e-3;
        double alpha0 = 0;
        // int i = 0;
        scanLineToFile(4*alpha1, 100, String.format("scanline_%04.4f_%02d.dat", t, i));
        lineSearchStep(alpha1);
        return 0.0;
    }
    void lineSearchStep(double alpha)  {
        add(1.0, X0, alpha, dX, X0);
        add(1/dt, X0, -1/dt, X1, V0);
        updateStressForceEnergy();
        calculatePhiGradient();
        dphi0 = dot(dX, G);
    }
    void calculateHessian() {
        add(1, D1, 1/dt, D2, temp2);
        add(1/dt/dt, M, temp2, H);
        // addIdentity(H, H, 1);
    }
    void explicitEulerStep() { }

    void implicitEulerStep() {
        // iterator move values, use D2 for lagged damping
        X2.set(X1);
        X1.set(X0);
        fDamp.zero();
        updateExternalForce();
        // initial status update
        mult(DC, D1, D2);
        add(1/dt, X1, -1/dt, X2, V0);
        mult(M, V0, temp1);
        energy = eInt + 0.5 * dot(temp1, V0);


        for(Cell c: cells) {
            cellStatusCallback.accept(c);
            cellForceCallback.accept(c);
        }
        bodyStatusCallback.accept(this);


        // constants
        double tau = 1e-2;
        double dE = 0;
        double dEOld = 1e20;

        // move all the values which are fixed during iteration to XHat
        add(2, X1, -1, X2, XHat);
        // shift XHat by constant force.
        mult(MI, fExt, temp1);
        add(1, XHat, dt*dt, temp1, XHat);
        // set initial value for X
        X0.set(XHat);
        System.out.format("Iteration start: t=%f\n", t);
        i = 0;
        while(i <= 10) {

            dE = dot(G, G);
            dE = stableNewtonStep(X0, XHat, V0, G, H);

            if(dE <= tau) {
                System.out.format("dE = %8.12f , phi0=%8.12f.\nStop iteration\n", dE, phi0);
                break;
            } else if (dE >= dEOld) {
                System.err.format("Energy minimization did not converge! dE=%g, dE'=%g\n",
                                  dE, dEOld);
                // assert(false);
                break;
            }

            dEOld = dE;
            // add(1/dt, X0, -1/dt, X1, V0);
            i++;
        }
    }
    void calculatePhiGradient() {
        phi0 = +eInt;
        add(1, X0, -1, XHat, G);
        mult(1/dt/dt, M, G, temp1);
        // temp1 <-- M(x-xHat)/dt^2
        phi0 += dot(G, temp1)/2.0; // phi += (x-xHat)M(x-xHat)/(2dt^2)
        mult(D2, V0, fDamp);
        phi0 += dot(V0, fDamp)*dt/2.0; // phi +=  v kD v dt/2
        // phi0 *= -1;
        add(-1, fInt, 1, fDamp, f);
        // mult(-dt*dt, MI, f, temp1);
        add(temp1, f, G);
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



    public double totalEnergy() {
        return energy;
    }
}
