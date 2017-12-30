package physics;
import java.nio.file.*;
import java.io.*;
import physics.*;
import java.util.ArrayList;
import static util.Math.*;
import static org.junit.Assert.*;
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
import org.apache.logging.log4j.*;
import static org.apache.logging.log4j.Level.*;

public class SoftBody {
    public static final boolean DEVEL = true;
    private static final Logger log = LogManager.getLogger();
    private static double testNu=0;
    private static double testD =0;
    // private static double averageNewtonIteration =0;
    // private static double averageLineSearchIteration =0;
    public DMatrixRMaj X0; // newly calculated state
    public DMatrixRMaj X1; // last state
    public DMatrixRMaj X2; // oldest state
    public DMatrixRMaj XHat, H;
    public DMatrixRMaj fInt, oldF, V0, vSquared, fExt, f, fDamp, dX, G;
    public DMatrixRMaj M, MI, DC; // mass, inverse mass, damping coefficient
    public DMatrixRMaj temp1, temp2;
    public DMatrixRMaj D1, D2;
    public double phi, dphi;
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
        //     if(DEVEL) log.printf(DEBUG, "UHessian: sx=%f, sy=%f, k=%f, E=%f, kdl=%f, fS=%f\n"
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
    final void addSymmetric(int i, int j, double x, DMatrixRMaj M) {
        double x0 = M.get(i, j);
        M.set(i, j, x0 + x);
        if(i != j)
            M.set(j, i, x0 + x);
    }


    // final void updateEnergyForce() {
    //     fInt.zero();
    //     eInt = 0;
    //     for(Bond bond: bonds) {
    //         int a = bond.a.index;
    //         int b = bond.b.index;
    //         double dx  = X1.get(b + 0, 0) - X1.get(a + 0, 0);
    //         double dy  = X1.get(b + 1, 0) - X1.get(a + 1, 0);
    //         double th1 = X1.get(a + 2, 0);
    //         double th2 = X1.get(b + 2, 0);
    //         double d = Math.sqrt(dx*dx + dy*dy);
    //         double l=bond.l, E=bond.E, k=bond.k;//, c=bond.c, D=bond.D;
    //
    //         double sx = dx /d/d;
    //         double sy = dy /d/d;
    //
    //         double nu = circleMod(Math.atan2(dy, dx) - bond.angle) - 0 * PI/2;
    //         if(bond.nu0 * nu < -PI) {
    //             if(bond.nu0 > nu) {
    //                 bond.rotationCounter++;
    //             } else {
    //                 bond.rotationCounter--;
    //             }
    //         }
    //         bond.nu0 = nu;
    //         nu += bond.rotationCounter * 2 * PI;
    //         testNu = nu;
    //         double fS = 3.0*E * (th1 + th2 - 2 * nu);
    //
    //         addEnergy(k, E, d, l, th1, th2, nu);
    //         addForce(k, E, d, l, th1, th2, nu, sx, sy, fS, a, b);
    //     }
    // }

    final void updateStressForceEnergy() {
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
            if(DEVEL) testNu = nu;
            if(DEVEL) testD =   d;
            double fS = 3.0*E * (th1 + th2 - 2 * nu);
            // System.out.println(d);
            assert d>0.05: "Cell distance too small: d="  + d;
            addEnergy(k, E, d, l, th1, th2, nu);
            addForce(k, E, d, l, th1, th2, nu, sx, sy, fS, a, b);
            addStress(k, E, d, l,sx, sy, fS, a, b);
        }
        // if(t >= t0 && t <= t0 + dt && s < 5) {
            // if(DEVEL) log.debug("Bond");
            // if(DEVEL) log.debug(bonds[0].toString());
            // if(DEVEL) log.debug("X1");
            // if(DEVEL) log.debug(X1.toString());
            // if(DEVEL) log.debug("v");
            // if(DEVEL) log.debug(v.toString());
            // if(DEVEL) log.debug("Hessian");
            // if(DEVEL) log.debug(D1.toString());
        // }
    }

    final DMatrixRMaj updateExternalForce() {
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

    final double newtonStep(){
        double dE;
        // mult(M, G, temp1);
        // negative sign because Nabla E = -f
        // caculate damping force
        // signs are positive because Nabla E = -f = -(-DC D2 v)
        dE = dot(G, G);
        add(1.0, X0, -1, G, X0);
        return dE;
    }

    final double stableNewtonStep() {
        double k = 1e-1;
        double alpha = 1;
        double alphaMax  = 1e3;

        // compute the step
        add(1/dt, X0, -1/dt, X1, V0);
        updateStressForceEnergy();
        calculatePhiGradient();
        double dE = dot(G, G);
        calculateHessian();
        invert(H, H);
        mult(-1, H, G, dX);
        // make sure dX is suitable
        if(DEVEL) log.printf(DEBUG, "dE = %8.12f , phi0=%8.12f, eInt=%8.12f\n", dE, phi, eInt);
        double dx = dot(dX, dX);
        dphi = dot(dX, G);
        double dXdE = dphi * dphi * Math.signum(dphi);
        if(dXdE < -k*dx*dE) {
            // dx is suitable
            // if(DEVEL) log.printf(DEBUG, "dX is suitable! dXdE=%g, dX dE=%g\n", dXdE, dx*dE);
        }
        else if(dXdE > k*dx*dE) {
            // -dx is suitable
            if(DEVEL) log.printf(DEBUG, "-dX is suitable! dXdE=%g, dX dE=%g\n", dXdE, dx*dE);
            dphi = - dphi;
            scale(-1, dX, dX);
        }
        else {
            if(DEVEL) log.printf(DEBUG, "Gradient is suitable! dXdE=%g, dX dE=%g\n", dXdE, dx*dE);
            scale(-1, G, dX);
            dphi = -dE;
            // dX.set(G);
        }
        strongWolfeLineSearch(alpha, alphaMax);
        dE = dot(G, G);
        return dE;
    }

    final double strongWolfeLineSearch(double alpha, double alphaMax) {
        double phiS  = phi;
        double dphiS = dphi;
        double alpha0 = 0;
        double alpha1 = alpha;
        double phi0 = phi;
        double dphi0 = dphi;
        double phi1 = phi;
        double dphi1 = dphi;
        double c1  = 1e-2;
        double c2  = 1e-2;
        scanLineToFile(10, 100, String.format("scanline_%04.4f_%02d.dat", t, i));
        int j = 1;
        // lineSearchStep(alpha);
        if(DEVEL) log.debug("Line search start.");
        while(j<5) {
            lineSearchStep(alpha1 - alpha0);
            phi1 = phi;
            dphi1 = dphi;
            if(Math.abs(phi) < 1e-20 && Math.abs(dphi) < 1e-20) return alpha1;
            if(phi1 > phiS + c1*alpha1*dphiS || (j>1 && phi1 >= phi0)) {
                if(DEVEL) log.printf(DEBUG, "Case phi is too big: phiS=%g, dphi=%g, dphiWolfe=%g\n",
                                  phiS, dphiS, phiS + c1*alpha*dphiS);
                alpha1 = zoom(alpha0, phi0, dphi0, alpha1, phi1, dphi1,
                              phiS, dphiS, c1, c2,
                              alpha0, phi0, dphi0, alpha1, phi1, dphi1);
                if(DEVEL) log.printf(DEBUG, "Zoom returned with alpha=%g\n", alpha1);
                break;
            }
            if(Math.abs(dphi1) <= -c2*dphiS) {
                if(DEVEL) log.printf(DEBUG, "Case phi is small enough and dphi is small enough:\n    phiS=%g, phi=%g, dphiS=%g, dphi=%g\n",
                                  phiS, phi1, dphiS, dphi1);
                if(DEVEL) log.printf(DEBUG, "LineSearch returned with alpha=%g\n", alpha1);
                break;
            }
            if(dphi1 >= 0) {
                if(DEVEL) log.printf(DEBUG, "Case phi is small and dphi is big and positive:\n    phiS=%g, dphiS=%g, phi=%g, dphi=%g\n", phiS, dphiS, phi, dphi);
                alpha1 = zoom(alpha1, phi1, dphi1, alpha0, phi0, dphi0,
                              phiS, dphiS, c1, c2,
                              alpha0, phi0, dphi0, alpha1, phi1, dphi1);
                if(DEVEL) log.printf(DEBUG, "Zoom returned with alpha=%g\n", alpha1);
                break;
            }
            if(DEVEL) log.printf(DEBUG, "Case phi is small and dphi is big and negative:\n    dphiS=%g, dphi=%g\n", dphiS, dphi);
            alpha0 = alpha1;
            phi0 = phi1;
            dphi0 = dphi1;
            alpha1 = choose(alpha1, phi, dphi, alphaMax);
            phi1 = phi;
            dphi1 = dphi;
            j++;
        }
        if(DEVEL) log.debug("Line search end.");
        return 0.0;
    }

    final double zoom(double lo, double philo, double dphilo,
                      double hi, double phihi, double dphihi,
                      double phiS, double dphiS, double c1, double c2,
                      double alpha0, double phi0, double dphi0,
                      double alpha1, double phi1, double dphi1) {
        int j = 0;
        double alpha;
        double e = 1e-6;
        if(DEVEL) log.debug("Zoom start.");
        while(j < 30) {
            alpha = interpolate(lo, philo, dphilo, hi, phihi, dphihi);
            if(DEVEL) log.printf(DEBUG, "Interp: lo=%g, philo=%g, dlo=%g,\n"+
                              "        hi=%g, phihi=%g, dhi=%g, alpha=%g\n",
                               lo, philo, dphilo, hi, phihi, dphihi, alpha);
            assert(alpha != Double.NaN);
            lineSearchStep(alpha - alpha1);
            if(phi > phiS + c1*alpha*dphiS || (phi >= philo)) {
                hi = alpha;
                phihi = phi;
                dphihi = dphi;
            } else {
                if(Math.abs(dphi) <= -c2*dphiS) {
                    return alpha;
                }
                if(dphi*(hi - lo) >= 0) {
                    hi = lo;
                    phihi = philo;
                    dphihi = dphilo;
                }
                lo = alpha;
                philo = phi;
                dphilo = dphi;
            }
            if(Math.abs(alpha - alpha1) < e) {
                return alpha;
            }
            alpha0 = alpha1;
            phi0 = phi1;
            dphi0 = dphi1;
            phi1 = phi;
            dphi1 = dphi;
            alpha1 = alpha;
            j++;
        }
        if(DEVEL) throw new RuntimeException("Zoom did not converge");
        return -1;
    }

    final static double interpolate(double a, double fa, double dfa,
                             double b, double fb, double dfb) {
        double c = 1e-1;
        double d = b - a;
        double s = Math.signum(d);
        double d1 = dfa + dfb - 3 * (fa  - fb)/(a-b);
        double k = d1*d1 - dfa*dfb;
        double d2 = s * Math.sqrt(k);
        double x= b - d * ((dfb + d2 - d1)/(dfb - dfa + 2*d2));
        // if(Math.abs(a-x) < c*d*s) {
        //     if(DEVEL) log.printf(DEBUG, "Halving Interval: a=%g x=%g b=%g", a, x, b);
        //     x = a + 2*c*d;
        // }
        // if(Math.abs(b-x) < c*d*s) {
        //     if(DEVEL) log.printf(DEBUG, "Halving Interval: a=%g x=%g b=%g", a, x, b);
        //     x = b + 2*c*d;
        // }
        // if(DEVEL) assertTrue("Interpolated alpha is not in interval.",
        //                      x <= Math.max(a, b) && x >= Math.min(a, b));
        return x;
    }

    final static double choose(double alpha1, double phi, double dphi,
                               double alphaMax){
            // TODO: make it smart
        return Math.min(2*alpha1, alphaMax);
    }

    final void lineSearchStep(double alpha)  {
        add(1.0, X0, alpha, dX, X0);
        add(1/dt, X0, -1/dt, X1, V0);
        updateStressForceEnergy();
        calculatePhiGradient();
        dphi = dot(dX, G);
    }

    void scanLineToFile(double alphaMax, int n, String filename) {
        try {
            Path path =  Paths.get("build/test-results/physics/SoftBody/ScanLine");
            if(Files.notExists(path)){
                    Files.createDirectories(path);
            }
            path = path.resolve(filename);
            BufferedWriter writer = Files.newBufferedWriter(path);
            writer.write(String.format("#alpha phi dphi1 dphi2 bondnu\n"));
            double phiInitial = phi;
            double phi1 = phi;
            double dphi1 = dphi;
            double phiMin = phi;
            double dphiMin = phi;
            double alphaMin = 0;
            double a = 0;
            double s = 1;
            double k = 1;
            double maxStep = 0.1;
            double minStep = 1e-5;
            phi1 = phi;
            double step = 5*alphaMax/n;
            lineSearchStep(s*step);
            a += s*step;
            double maxDPhiDiff = Math.abs(dphi1 - dphi);
            for(int l=0; l<3; l++) {
                while(Math.abs(a) < alphaMax * (l == 2 ? 0:1)) {
                    if(phi < phiMin) {
                        phiMin = phi;
                        dphiMin = dphi;
                        alphaMin = a;
                    }
                    phi1 = phi;
                    dphi1 = dphi;
                    lineSearchStep(s*step);
                    a += s*step;
                    writer.write(String.format("%g %12.12f %12.12f %12.12f %12.12f %12.12f\n",
                                 a, phi, dphi, s*(phi - phi1)/step, testNu, testD));
                    if(Math.abs(dphi1 - dphi) < maxDPhiDiff) {
                        if(step < maxStep) {
                            step /= k;
                        }
                    } else {
                        if(step > minStep) {
                            step *= k;
                        }
                    }
                }
                step /= k;
                lineSearchStep(-s*step);
                a += -s*step;
                phi1 = phi;
                dphi1 = dphi;
                s=-s;
            }
            lineSearchStep(-a);
            if(DEVEL) log.printf(DEBUG, "Line Scan min value: alphaMin=%4.4f, phiMin=%g, dphiMin=%g\n", alphaMin, phiMin, dphiMin);
            writer.close();
            assertEquals(String.format("Line scan not reversed: AlphaMax=%g, n=%d, t=%g, i=%d\n", alphaMax, n, t, i),
                         phiInitial, phi, 1e-6);
        } catch(IOException e) {
            throw new RuntimeException("Could not write data");
        }
    }

    void calculateHessian() {
        add(1, D1, 1/dt, D2, temp2);
        add(1/dt/dt, M, temp2, H);
        // addIdentity(H, H, 1);
    }

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



        // constants
        double tau = 1e-2;
        double dE = 0;
        double phiOld = Double.MAX_VALUE;

        // move all the values which are fixed during iteration to XHat
        add(2, X1, -1, X2, XHat);
        // shift XHat by constant force.
        mult(MI, fExt, temp1);
        add(1, XHat, dt*dt, temp1, XHat);
        // set initial value for X
        X0.set(XHat);
        if(DEVEL) log.printf(DEBUG, "Iteration start: t=%f\n\n", t);
        i = 0;
        while(i <= 40) {
            dE = stableNewtonStep();
            if(dE <= tau) {
                if(DEVEL) log.printf(DEBUG, "dE = %8.12f , phi=%8.12f.\nStop iteration\n", dE, phi);
                energy = eInt + 0.5 * dot(temp1, V0);
                for(Cell c: cells) {
                    cellStatusCallback.accept(c);
                    cellForceCallback.accept(c);
                }
                bodyStatusCallback.accept(this);
                return;
            }
            if(DEVEL) {
                if (phi >= phiOld) {
                    assertTrue(
                        String.format("Energy minimization did not converge! phi=%g, phiOld=%g\n",
                                      phi, phiOld), false);

                    // assert(false);
                    break;
                }
            }

            phiOld = phi;
            // add(1/dt, X0, -1/dt, X1, V0);
            i++;
        }
        if(DEVEL) assertTrue(
            String.format("Energy minimization did not converge after 40 iterations! dE=%g, dE'=%g\n",
                          dE, dE),
            false
        );
    }
    void calculatePhiGradient() {
        phi = +eInt;
        add(1, X0, -1, XHat, G);
        mult(1/dt/dt, M, G, temp1);
        // temp1 <-- M(x-xHat)/dt^2
        phi += dot(G, temp1)/2.0; // phi += (x-xHat)M(x-xHat)/(2dt^2)
        mult(D2, V0, fDamp);
        phi += dot(V0, fDamp)*dt/2.0; // phi +=  v kD v dt/2
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
            // if(DEVEL) log.debug("y" + y);
            for(int col = 0; col<cells[row].length; col++) {
                // if(DEVEL) log.debug("x" + x);
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
