package physics;
import physics.*;

import java.nio.file.*;
import java.io.*;

import static util.Math.*;
import static org.junit.Assert.*;
import static org.ejml.dense.row.CommonOps_DDRM.*;
import static org.ejml.dense.row.NormOps_DDRM.*;
import org.ejml.data.*;

import org.apache.logging.log4j.*;
import static org.apache.logging.log4j.Level.*;

import static main.Main.DEVEL;

abstract class ImplicitODESolver {

    final void computeNewtonDirection(DMatrixRMaj x, DMatrixRMaj dn) {
        phi = computeOptimizeGradient(x, g);
        scale(-1, g, dn);
        computeDPhi(g, dn);
    };
    final void computeDPhi(DMatrixRMaj g, DMatrixRMaj dn) {
        dPhi = dot(g, dn)/normP2(dn);
    }

    final double computeOptimizeGradient(DMatrixRMaj x, DMatrixRMaj dest) {
        double energy = computeForce(x, dest);
        add(1, x, -1, xHat, temp1);
        elementMult(M, temp1, temp2);
        // energy += dot(temp1, temp2)/2/dt/dt;
        // add(1/dt/dt, temp2, -1, dest, dest);
        scale(-1, dest, g);
        return energy;
    }


    abstract double computeForce(DMatrixRMaj x, DMatrixRMaj dest);
    abstract void computeForceDifferential(DMatrixRMaj x, DMatrixRMaj dx, DMatrixRMaj dest);
    abstract void timeStepFinished();

    DMatrixRMaj dn, g, x0, x1, x2, xHat, temp1, temp2, M, MI, fExt, xAlpha, gConst;
    double dG, dPhi, dX, phi, dN;
    static double dt, t;
    int iNewton=0;

    private static final Logger log = LogManager.getLogger();
    ImplicitODESolver(int n) {
        dn = new DMatrixRMaj(n, 1);
        g = new DMatrixRMaj(n, 1);
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

        computeNewtonDirection(x0, dn);

        dN = Math.sqrt(dot(dn, dn));
        dG = Math.sqrt(dot(g, g));
        if(dPhi < -k*dN*dG) {
            // dN is suitable
            // if(DEVEL) log.printf(DEBUg, "dn is suitable! dNdG=%g, dn dg=%g\n", dNdG, dN*dg);
        }
        else if(dPhi > k*dN*dG) {
            // -dN is suitable
            if(DEVEL) log.printf(DEBUG, "-dn is suitable! dNdG=%g, dn dg=%g\n", dPhi, dN*dG);
            dPhi = - dPhi;
            scale(-1, dn, dn);
        }
        else {
            if(DEVEL) log.printf(DEBUG, "gradient is suitable! dNdG=%g, dn dg=%g\n", dPhi, dN*dG);
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
        scanLineToFile(8*alpha, 100, String.format("scanline_%04.4f_%02d.dat", t, iNewton));
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
                if(DEVEL) log.printf(DEBUG, "Zoom returned with alpha=%g\n", alpha1);
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
                if(DEVEL) log.printf(DEBUG, "Zoom returned with alpha=%g\n", alpha1);
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
            log.debug("Origin and direction" + x0 + dn);
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
        // elementMult(fExt, MI, temp1);
        // add(dt*dt, temp1, x0, x0);
    }

    void implicitEulerStep() {
        double tau = 1e-2;
        double dG = 1e33;
        double phiOld = Double.MAX_VALUE;

        // step forward and set the initial guess
        computeForwardEulerStep(x2, x1, x0, fExt);
        xHat.set(x0);

        if(DEVEL) log.printf(DEBUG, "Newton iteration start. t=%g", t);
        iNewton=0;
        while(iNewton <= 40) {
            if(dG <= tau) {
                return;
            }
            dG = newtonStep();
            if(DEVEL) {
                if (phi >= phiOld) {
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
            if(DEVEL) log.printf(DEBUG, "Newton iteration i=%d, t=%g", iNewton ,t);
        }
        if(DEVEL) assertTrue(
            String.format("Energy minimization did not stop after 40 iterations!" +
                          " dE=%g, dE'=%g\n", dG, dG),
            false
        );
    }

}
