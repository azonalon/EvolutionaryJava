package physics;

public class Leapfrog {
    abstract public static class DampedLeapfrog {
        double x1, x2, v1, v2, dt, t;
        /* position, velocity, mass, time step, damping constant */

        abstract double F(double x, double t);
        abstract double mass();
        abstract double damping();

        /**
        *  @param dt time step
        *  dt should be the same over the whole integration,
        *  otherwise it will not give stable results
        */
        public void step() {
            double c = damping();
            double m = mass();
            x1 = x2;
            v1 = v2;
            v2 = ((1.0 - c * dt/(2.0*m)) * v1 +  dt * F(x1, t)/m)/(1.0 + c * dt/(2.0*m));
            x2 = x1 + dt * v2;
        }
        public void setVelocity(double v) {
            this.v2 = v;
        }
        public double v() {
            return this.v2;
        }
        public void stepBack() {
            v2 = v1;
            x2 = x1;
        }
        public double x() {
            return this.x2;
        }

        public void init (double x0, double v0, double dt, double t){
            this.dt = dt;
            v2 = v0 - dt/2 * F(x0, t)/mass();
            x2 = x0;
        }

        public void integrate(double v0, double x0, double t0, double dt, int nSteps) {
            t = t0;
            init(x0, v0, dt, t);
            for(int i = 0; i<nSteps; i++) {
                step();
                t += dt;
                System.out.format("%4.4f %4.4f %4.4f\n", t, x2, v2);
                // System.out.format("%4.4f %4.4f %4.4f\n", t, exactSolution(t, x0), v2);
            }
        }
    }

    public static class DampedHarmonicOscillator extends DampedLeapfrog {
        double omega1, omega0, k, c, m;

        /**
        * @param m particle mass
        * @param zeta damping factor with gamma = 2 * omega0 * zeta in harmonic oscillator
        */
        DampedHarmonicOscillator(double m, double zeta, double omega0) {
            this.m = m;
            k = omega0 * omega0 * m;
            c = 2  * zeta * omega0;
            omega1 = Math.sqrt(1-zeta*zeta) * omega0;
        }

        double exactSolution(double t, double x0) {
            return x0 * Math.exp(-c/2.0*t)
            * (1 + t * omega0);
        }

        @Override
        double F(double x, double t)  {
            return -k * x;
        }
        @Override double damping() {
            return this.c;
        }
        @Override double mass() {
            return this.m;
        }
    }
    public static void main(String[] args) {
        DampedHarmonicOscillator l = new Leapfrog.DampedHarmonicOscillator(1, 0, 2 * Math.PI);
        l.integrate(0, 1, 0, 0.1, 20);
    }
}
