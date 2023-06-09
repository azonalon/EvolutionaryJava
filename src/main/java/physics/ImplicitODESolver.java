package physics;
import physics.*;

import java.nio.file.*;
import java.io.*;
import java.util.function.*;

import static util.Math.*;
import static org.junit.Assert.*;
import static org.ejml.dense.row.CommonOps_DDRM.*;
import static org.ejml.dense.row.NormOps_DDRM.*;
import org.ejml.data.*;

import org.apache.logging.log4j.*;
import static org.apache.logging.log4j.Level.*;

import static main.Main.DEVEL;

abstract class ImplicitODESolver {
    abstract double computeForce(DMatrixRMaj x, DMatrixRMaj dest);
    abstract void computeForceDifferential(DMatrixRMaj x, DMatrixRMaj dx, DMatrixRMaj dest);
    abstract void timeStepFinished();
    public double dampPotential;

    DMatrixRMaj dn, g, x0, x1, x2, v, xHat,
                temp1, temp2, M, MI, fExt,
                xAlpha, gConst, r, p;
    double dG, dPhi, dX, phi, dN;
    double kDamp;
    static double dt, t;
    int iNewton=0;
    double newtonAccuracy = 1e-5;

    private static final Logger log = LogManager.getLogger();

    final void conjugateGradientSolve(DMatrixRMaj rhs, DMatrixRMaj initialGuess,
                                      BiConsumer<DMatrixRMaj,DMatrixRMaj>
                                        computeLhs,
                                      DMatrixRMaj result) {
        computeLhs.accept(initialGuess, temp1); // temp1=Ax0
        add(1, rhs, -1, temp1, r); // temp2 = p0 = r0 = b - Ax0
        p.set(r);
        double rr1, rr2, alpha, beta;
        computeLhs.accept(p, temp1); // temp1 = Ap
        rr1 = dot(r,r);
        alpha=rr1/dot(p, temp1);
        add(1, initialGuess, alpha, p, result); // result = x1 = x0 + alpha p0

        if(DEVEL) log.debug("Conjugate Gradient Start.");
        for(int k=0; k<rhs.getNumRows(); k++) {
            // if(alpha < 0) {
            //     if(DEVEL) log.debug("Indefiniteness detected.");
            //     if(DEVEL) log.debug("Conjugate gradient stop.");
            //     break;
            // }
            add(1, r, -alpha, temp1, r);
            rr2 = dot(r, r);
            if(DEVEL) log.printf(DEBUG, "Conjugate Gradient: k=%d, ||r||=%g", k, rr2);
            if(rr2 < 1e-5 ) {
                if(DEVEL) log.debug("Conjugate gradient stop.");
                break;
            }
            beta = rr2/rr1;
            rr1 = rr2;
            add(1,r,beta, p, p);
            computeLhs.accept(p, temp1);
            alpha = rr1/dot(p, temp1);
            add(1, result, alpha, p, result);
        }
    }

    final void computeNewtonDirection(DMatrixRMaj g, DMatrixRMaj dn) {
        conjugateGradientSolve(g, g, (dx, lhs)->{
            assert(normP2(x0) > 0);
            assert(normP2(x1) > 0);
            computeForceDifferential(x0, dx, lhs);
            computeForceDifferential(x1, dx, temp2);
            add(-kDamp/dt, temp2, 1, lhs, lhs);
            elementMult(M, dx, temp2);
            add(1/dt/dt, temp2, 1, lhs, lhs);
        }, dn);
        scale(-1, dn, dn);
    };

    final void computeDPhi(DMatrixRMaj g, DMatrixRMaj dn) {
        dPhi = dot(g, dn);
    }

    final double computeOptimizeGradient(DMatrixRMaj x, DMatrixRMaj dest) {
        double energy = computeForce(x, dest);

        add(1, x, -1, xHat, temp1);
        elementMult(M, temp1, temp2);
        energy += dot(temp1, temp2)/2/dt/dt;
        add(1/dt/dt, temp2, -1, dest, dest);

        add(1/dt, x, -1/dt, x1, v);
        computeForceDifferential(x1, v, temp2);
        dampPotential = -dot(v, temp2)*kDamp/2*dt;
        energy += dampPotential;
        // if(DEVEL) log.printf(DEBUG, "Compute Gradient: damping dampPotential=%g", dampPotential);
        add(-kDamp, temp2, 1, dest, dest);
        return energy;
    }


    ImplicitODESolver(int n) {
        dn = new DMatrixRMaj(n, 1);
        g = new DMatrixRMaj(n, 1);
        r = new DMatrixRMaj(n, 1);
        p = new DMatrixRMaj(n, 1);
        v = new DMatrixRMaj(n, 1);
        x0 = new DMatrixRMaj(n, 1);
        x1 = new DMatrixRMaj(n, 1);
        x2 = new DMatrixRMaj(n, 1);
        xHat = new DMatrixRMaj(n, 1);
        xAlpha = new DMatrixRMaj(n, 1);
        gConst = new DMatrixRMaj(n, 1);
        fExt = new DMatrixRMaj(n, 1);
        fExt.zero();
        temp1 = new DMatrixRMaj(n, 1);
        temp2 = new DMatrixRMaj(n, 1);
    }

    final double newtonStep() {
        double k = 1e-2;
        double alpha = 1;
        double alphaMax  = 1e3;
        double l = 1e3;

        phi = computeOptimizeGradient(x0, g);
        dG = dot(g, g);
        if(dG < newtonAccuracy) {
            return dG;
        }
        computeNewtonDirection(g, dn);

        dN = dot(dn, dn);
        computeDPhi(g, dn);

        double dPhiSq = dPhi*Math.abs(dPhi);
        if(dPhiSq < -k*dN*dG) {
            // dN is suitable
            // if(DEVEL) log.printf(DEBUg, "dn is suitable! dNdG=%g, dn dg=%g\n", dNdG, dN*dg);
        }
        else if(dPhiSq > k*dN*dG) {
            // -dN is suitable
            if(DEVEL) log.printf(DEBUG, "-dn is suitable! dPhi*abs(dPhi)=%g, dN*dG=%g\n", dPhiSq, dN*dG);
            dPhi = -dPhi;
            scale(-1, dn, dn);
        }
        else {
            if(DEVEL) log.printf(DEBUG, "gradient is suitable! dNdG=%g, dN*dG=%g\n", dPhi*dPhi, dN*dG);
            scale(-1, g, dn);
            dPhi = -dG;
            dN = dG;
            // dn.set(g);
        }
        if(dN > l) {
            scale(l/dN, dn, dn);
            dN = l;
        }
        alpha = strongWolfeLineSearch(alpha, alphaMax);
        add(1, x0, alpha, dn, x0);
        dG = dot(g, g);
        return dG;
    }

    final double strongWolfeLineSearch(double alpha, double alphaMax) {
        double phiS  = phi;
        double dPhiS = dPhi;
        double alpha0 = 0;
        double alpha1 = alpha;
        double phi0 = phi;
        double dPhi0 = dPhi;
        double phi1 = phi;
        double dPhi1 = dPhi;
        double c1  = 1e-2;
        double c2  = 1e-2;
        // scanLineToFile(5*alpha, 100, String.format("scanline_%04.4f_%02d.dat", t, iNewton));
        int j = 1;
        // lineSearchStep(alpha);
        if(DEVEL) log.debug("Line search start.");
        while(j<5) {
            computePhiDPhi(alpha1);
            phi1 = phi;
            dPhi1 = dPhi;
            if(Math.abs(phi) < 1e-20 && Math.abs(dPhi) < 1e-20) return alpha1;
            if(phi1 > phiS + c1*alpha1*dPhiS || (j>1 && phi1 >= phi0)) {
                if(DEVEL) log.printf(DEBUG, "Case phi is too big: phiS=%g, dPhi=%g, phiWolfe=%g\n",
                                  phiS, dPhiS, phiS + c1*alpha*dPhiS);
                alpha1 = zoom(alpha0, phi0, dPhi0, alpha1, phi1, dPhi1,
                              phiS, dPhiS, c1, c2,
                              alpha0, phi0, dPhi0, alpha1, phi1, dPhi1);
                if(DEVEL) log.printf(DEBUG, "Zoom returned with alpha=%g, phi=%g, dphi=%g\n", alpha1, phi, dPhi);
                return alpha1;
            }
            if(Math.abs(dPhi1) <= -c2*dPhiS) {
                if(DEVEL) log.printf(DEBUG, "Case phi is small enough and dPhi is small enough:\n    phiS=%g, phi=%g, dPhiS=%g, dPhi=%g\n",
                                  phiS, phi1, dPhiS, dPhi1);
                if(DEVEL) log.printf(DEBUG, "LineSearch returned with alpha=%g\n", alpha1);
                return alpha1;
            }
            if(dPhi1 >= 0) {
                if(DEVEL) log.printf(DEBUG, "Case phi is small and dPhi is big and positive:\n    phiS=%g, dPhiS=%g, phi=%g, dPhi=%g\n", phiS, dPhiS, phi, dPhi);
                alpha1 = zoom(alpha1, phi1, dPhi1, alpha0, phi0, dPhi0,
                              phiS, dPhiS, c1, c2,
                              alpha0, phi0, dPhi0, alpha1, phi1, dPhi1);
                if(DEVEL) log.printf(DEBUG, "Zoom returned with alpha=%g, phi=%g, dphi=%g\n", alpha1, phi, dPhi);
                return alpha1;
            }
            if(DEVEL) log.printf(DEBUG, "Case phi is small and dPhi is big and negative:\n    dPhiS=%g, dPhi=%g\n", dPhiS, dPhi);
            alpha0 = alpha1;
            phi0 = phi1;
            dPhi0 = dPhi1;
            alpha1 = choose(alpha1, phi, dPhi, alphaMax);
            phi1 = phi;
            dPhi1 = dPhi;
            j++;
        }
        if(DEVEL) log.debug("Line search end.");
        if(DEVEL) throw new RuntimeException("Line Search did not end");
        return 0;
    }

    final double zoom(double lo, double philo, double dPhilo,
                      double hi, double phihi, double dPhihi,
                      double phiS, double dPhiS, double c1, double c2,
                      double alpha0, double phi0, double dPhi0,
                      double alpha1, double phi1, double dPhi1) {
        int j = 0;
        double alpha;
        double e = 1e-6;
        if(DEVEL) log.debug("Zoom start.");
        while(j < 30) {
            alpha = interpolate(lo, philo, dPhilo, hi, phihi, dPhihi);
            if(DEVEL) log.printf(DEBUG, "Interp: lo=%g, philo=%g, dlo=%g,\n"+
                              "        hi=%g, phihi=%g, dhi=%g, alpha=%g\n",
                               lo, philo, dPhilo, hi, phihi, dPhihi, alpha);
            assert(alpha != Double.NaN);
            computePhiDPhi(alpha);
            if(phi > phiS + c1*alpha*dPhiS || (phi >= philo)) {
                hi = alpha;
                phihi = phi;
                dPhihi = dPhi;
            } else {
                if(Math.abs(dPhi) <= -c2*dPhiS) {
                    return alpha;
                }
                if(dPhi*(hi - lo) >= 0) {
                    hi = lo;
                    phihi = philo;
                    dPhihi = dPhilo;
                }
                lo = alpha;
                philo = phi;
                dPhilo = dPhi;
            }
            if(Math.abs(alpha - alpha1) < e) {
                return alpha;
            }
            alpha0 = alpha1;
            phi0 = phi1;
            dPhi0 = dPhi1;
            phi1 = phi;
            dPhi1 = dPhi;
            alpha1 = alpha;
            j++;
        }
        if(DEVEL) throw new RuntimeException("Zoom did not converge");
        return -1;
    }

    final static double interpolate(double a, double fa, double dfa,
                             double b, double fb, double dfb) {
        double c = 5e-2;
        double d = b - a;
        double s = Math.signum(d);
        double d1 = dfa + dfb - 3 * (fa  - fb)/(a-b);
        double k = d1*d1 - dfa*dfb;
        double d2 = s * Math.sqrt(k);
        double x= b - d * ((dfb + d2 - d1)/(dfb - dfa + 2*d2));
        if(Math.abs(a-x) < c*d*s) {
            if(DEVEL) log.printf(DEBUG, "Halving Interval: a=%g x=%g b=%g", a, x, b);
            x = a + 2*c*d;
        }
        if(Math.abs(b-x) < c*d*s) {
            if(DEVEL) log.printf(DEBUG, "Halving Interval: a=%g x=%g b=%g", a, x, b);
            x = b + 2*c*d;
        }
        if(DEVEL) assertTrue("Interpolated alpha is not in interval.",
                             x <= Math.max(a, b) && x >= Math.min(a, b));
        return x;
    }

    final static double choose(double alpha1, double phi, double dphi,
                               double alphaMax){
            // TODO: make it smart
        return Math.min(2*alpha1, alphaMax);
    }

    final void computePhiDPhi(double alpha)  {
        add(alpha, dn, x0, xAlpha);
        phi = computeOptimizeGradient(xAlpha, g);
        // dPhi = dot(dn, g)/normP2(dn);
        computeDPhi(g, dn);
    }

    void scanLineToFile(double alphaMax, int n, String filename) {
        assert alphaMax >0;
        try {
            Path path =  Paths.get("build/test-results/physics/ScanLine");
            if(Files.notExists(path)){
                    Files.createDirectories(path);
            }
            path = path.resolve(filename);
            BufferedWriter writer = Files.newBufferedWriter(path);
            writer.write(String.format("#alpha phi dPhi1 dPhi2\n"));
            double phi1 = phi;
            double phiMin = phi;
            double dPhiMin = dPhi;
            double alphaMin = 0;
            double step = alphaMax/n;
            computePhiDPhi(-alphaMax - step);
            for(double alpha=-alphaMax; alpha<alphaMax; alpha+=step) {
                phi1 = phi;
                computePhiDPhi(alpha);
                if(phiMin > phi) {
                    dPhiMin = dPhi;
                    phiMin = phi;
                    alphaMin = alpha;
                }
                writer.write(String.format("% 8g % 12.8f % 12.8f % 12.8f\n",
                             alpha, phi, dPhi, (phi - phi1)/step));
            }
            computePhiDPhi(0);
            if(DEVEL) log.printf(DEBUG, "Line Scan min value: alphaMin=%4.4f,"+
                                 " phiMin=%g, dPhiMin=%g\n", alphaMin, phiMin, dPhiMin);
            // log.debug("Origin and direction" + x0 + dn);
            writer.close();
        } catch(IOException e) {
            throw new RuntimeException("Could not write data");
        }
    }

    void computeForwardEulerStep(DMatrixRMaj x2,DMatrixRMaj x1,DMatrixRMaj x0,
                                 DMatrixRMaj f) {
        x2.set(x1);
        x1.set(x0);
        add(2, x1, -1, x2, x0);
        elementMult(fExt, MI, temp1);
        add(dt*dt, temp1, x0, x0);
        add(1/dt, x1, -1/dt, x2, temp1);
        if(DEVEL) log.printf(DEBUG, "\n||v||=%g", normP2(temp1));

    }

    void implicitEulerStep() {
        double dG = 1e33;
        double phiOld = Double.MAX_VALUE;

        if(DEVEL) log.printf(DEBUG, "\nNewton iteration start. t=%g", t, dG);
        // step forward and set the initial guess
        computeForwardEulerStep(x2, x1, x0, fExt);
        xHat.set(x0);

        iNewton=0;
        while(iNewton <= 20) {
            dG = newtonStep();
            if(dG < newtonAccuracy) {
                if(DEVEL) log.printf(DEBUG, "Newton iteration stopped: i=%d, t=%g, dG=%g",
                 iNewton ,t, dG);
                return;
            }
            if(DEVEL) {
                if (phi >= phiOld*1.2) {
                    assertTrue(
                        String.format("Energy minimization did not converge!"+
                                      "phi=%g, phiOld=%g\n",
                                      phi, phiOld),
                        false
                    );
                    break;
                }
            }
            phiOld = phi;
            iNewton++;
            if(DEVEL) log.printf(DEBUG, "\nNewton iteration i=%d, t=%g, dG=%g", iNewton ,t, dG);
        }
        if(DEVEL) assertTrue(
            String.format("Energy minimization did not stop after 40 iterations!" +
                          " dE=%g, dE'=%g\n", dG, dG),
            false
        );
    }

}
