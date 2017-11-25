package physics.test;
import java.util.*;
import physics.*;

class Collision {
    List<TreeSet<Particle>> collisionGrid;
    Particle[] particles;
    double t;
    static int collisionBinCount = 100; // should be divisible by two
    static double gridrange = 10;
    static double binsPerGridrange = (double)collisionBinCount/gridrange;
    static double binDiameter = gridrange/collisionBinCount;

    public Collision() {
        collisionGrid = new ArrayList<TreeSet<Particle>>(collisionBinCount);
        for(int i = 0; i < collisionBinCount * 2; i++) {
            collisionGrid.add(new TreeSet<Particle>());
        }
        // binsPerGridrange = (double)collisionGridSize / gridrange;
        // binDiameter = 1/binsPerGridrange;
    }

    /**
     * Returns the index of a given point to the grid
     * This will be very inefficient if there are many particles * outside of gridrange, because particles being very far away from
     * each other will be checked against each other too.
     * @param  double x  point in space
     * @return index into grid array
     */
    public int collisionGridIndex(double x) {
        int index =  (int)((x % gridrange) * binsPerGridrange) + collisionBinCount;
        // System.out.println("Index " + index);
        return index;
    }

    static class Particle extends Leapfrog.DampedLeapfrog
                          implements Comparable<Particle>{
        double m;
        public Particle(double mass) {
            m = mass;
        }
        double mass() {return m;}
        double damping() {return 0.0;}
        double F(double x, double t) { return 0; }
        double radius() {return binDiameter/4.0;}
        public int compareTo(Particle other) {
            return Double.compare(this.x(), other.x());
        }
    }

    /**
     * output for gnuplot
     */
    public void printParticles() {
        System.out.format("%4.4f ", t);
        for(Particle p: particles) {
            System.out.format("%4.4f ", p.x());
        }
        System.out.println();
    }

    public void collisionCheck() {
        for(Particle p: particles) {
            collisionGrid.get(collisionGridIndex(p.x() + p.radius())).add(p);
            collisionGrid.get(collisionGridIndex(p.x() - p.radius())).add(p);
        }
        for(TreeSet<Particle> collided: collisionGrid) {
            if(collided.size() > 1) {
                Iterator<Particle> it = collided.iterator();
                Particle p1, p2;
                p1 = it.next();
                while(it.hasNext()) {
                    p2 = it.next();
                    if(Math.abs(p1.x() - p2.x()) - (p1.radius() + p2.radius()) < 0) {
                        System.out.println("Collision at t=" + t);
                        resolveCollision(p1, p2);
                    }
                    p1 = p2;
                }
            }
            collided.clear();
        }
    }
    void resolveCollision(Particle p1, Particle p2) {
        p1.stepBack();
        p2.stepBack();
        double newV1 = (p1.v() * (p1.mass() - p2.mass()) + 2 * p2.mass() * p2.v()) /
                       (p1.mass() + p2.mass());
        double newV2 = (p2.v() * (p2.mass() - p1.mass()) + 2 * p1.mass() * p1.v()) /
                       (p1.mass() + p2.mass());
        p1.setVelocity(newV1);
        p2.setVelocity(newV2);
        p1.step();
        p2.step();
    }

    public static void main(String[] args) {
        int nParticles = 5;
        Random r = new Random(17);
        Collision c = new Collision();
        c.t = 0;
        double dt = 0.01;
        int nSteps = 900;
        c.particles = new Particle[nParticles + 2];
        // double[] masses = new double[] {1, 1};
        Particle pStart = new Particle(100000.0);
        Particle pEnd =   new Particle(100000.0);
        pStart.init(-2, 0, dt, c.t);
        pEnd.init(+2, 0, dt, c.t);
        c.particles[nParticles] = pStart;
        c.particles[nParticles + 1] = pEnd;
        for(int i = 0; i<nParticles; i++) {
            Particle p = new Particle(1.0);
            p.init(r.nextDouble()*2-1, r.nextDouble()*2-1, dt, c.t);
            p.init(r.nextDouble()*2-1, r.nextDouble()*2-1, dt, c.t);
            c.particles[i] = p;
        }


        System.out.println("Collision grid d=" + Collision.binDiameter);
        // assert(c.particles[0].radius() * 2 > binDiameter);
        for(int i = 0; i<nSteps; i++) {
            for(Particle p: c.particles) {
                    p.step();
                    // this line makes sure the particles are not passing
                    // each other during a single time step
                    assert(p.v() * dt < 2 * p.radius());
                }
            c.collisionCheck();
            c.t += dt;
            c.printParticles();
        }
    }
}
