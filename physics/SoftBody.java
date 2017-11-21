package physics;
import physics.*;
import java.util.*;
import java.util.function.*;
import static java.lang.Math.PI;
import static java.lang.Math.sqrt;
// import static System.out.format;

public class SoftBody {
    static double sqrd(double a) {
        return a*a;
    }
    Cell[] cells;
    int [] cellStates; // these are used in BFS 'update forces' to check
                       // if bond force was calculated
    static int LEFT = 1 << 0;
    static double phi0;
    static int RIGHT= 1 << 1;
    static int UP   = 1 << 2;
    static int DOWN = 1 << 3;
    static int counter=0;
    static double dt = 0.1;
    double energy = 0;
    static boolean startProgram = true;
    static Function<Cell, String> cellStatusReading =
        (Cell c)->{return String.format("%f ", c.x);};
    static Function<SoftBody, String> bodyStatusReading =
        (body)->{return body.totalEnergy() + "";};

    SoftBody(Cell[] cells) {
        this.cells = cells;
        cellStates = new int[cells.length];
    }
    static double circleMod(double angle) {
        return angle % PI;
    }
    static double harmonicMean(double a, double b) {
        return 2 * a * b / (a + b + 1e-9);
    }

    /**
     * checks if the bond forces between the cells have been added to each cell.
     * this is a symmetric operation a <--> b
     * @param  Cell a             [description]
     * @param  Cell b             [description]
     * @return      [description]
     */

    void propagateCells() {
        for(Cell c: cells) {
            System.out.print(cellStatusReading.apply(c) + " ");
            energy += c.L * c.L * c.I /2.0 + c.vx * c.vx * c.m /2.0 + c.vy * c.vy * c.m /2.0;
            c.propagate();
        }
        System.out.print(bodyStatusReading.apply(this));
        System.out.println();
    }

    void visit(Cell a, int side) {
        assert(cellStates[a.index] & side) == 0; // dont visit twice
        cellStates[a.index] |= side;
    }
    boolean unvisited(Cell a, int side) {
        return (cellStates[a.index] & side) == 0;
    }
    boolean unvisited(Cell a) {
        return (cellStates[a.index]) == 0;
    }
    void updateForces() {
        Cell a = cells[0];
        a.fX = 0;
        a.fY = 0;
        a.T  = 0;
        energy = 0;
        Queue<Cell> q = new ArrayDeque<Cell>();
        Cell[] nextCells = new Cell[] {a.left, a.right, a.up, a.down};
        int[] directions = new int[] {RIGHT, LEFT, DOWN, UP};
        while(a != null) {
            for(int i =0; i<4; i++ ) {
                Cell next = nextCells[i];
                int direction = directions[i];
                if(next != null) {
                    if(unvisited(next)) {
                        next.fX = 0;
                        next.fY = 0;
                        next.T  = 0;
                        q.offer(next);
                    }
                    if(unvisited(next, direction)) {
                        addBondForces(a, next);
                        counter++;
                        visit(a, direction);
                    }
                }
            }
            a = q.poll();
        }
        // clean up
        Arrays.fill(cellStates, 0);
    };

    void addBondForces(Cell first, Cell second) {

        double dx  = - first.x + second.x;
        double dy  = - first.y + second.y;
        double d = Math.sqrt(dx*dx + dy*dy);
        double phi = Math.asin(dy/d);
        double dphi = phi - phi0;
        phi0 = phi;

        // double tTheta =  -first.theta + second.theta + phi;
        // rest length of the spring
        double k = harmonicMean(first.k, second.k);
        double E = harmonicMean(first.E, second.E);
        double D = harmonicMean(first.D, second.D);
        double c = harmonicMean(first.c, second.c);
        double l = first.r + second.r;
        double fShear = - 6 * l * l * E * circleMod(first.theta + second.theta - 2 * phi);
        dx = dx / d;
        dy = dy / d;


        // System.err.format("%f %f %f\n", phi, first.theta, second.theta);
        first.fX  += (d - l) * k * dx - fShear * (-1*dy);
        first.fY  += (d - l) * k * dy - fShear * (   dx);
        second.fX -= (d - l) * k * dx + 1 * fShear * (-1*dy);
        second.fY -= (d - l) * k * dy + 1 * fShear * (   dx);

        // TODO:  damp relative rotation
        first.T  += 2 * E * l * l * l * circleMod(2 * first.theta + second.theta - 3 * phi);
        second.T -= 2 * E * l * l * l * circleMod(2 * second.theta + first.theta - 3 * phi);


        energy += (d-l)*(d-l)* k/2.0;
        energy += 2 * l * l * l * E * sqrd(first.theta + second.theta + 2 *phi);

        if(startProgram == true) {
            System.err.format("Parameters: k=%f, E=%f, D=%f, c=%f, l=%f\n",
                                k, E, D, c, l);
            startProgram = false;
        }


        // damping forces
        double vXRelative = first.vx - second.vx;
        double vYRelative = first.vy - second.vy;
        // subtract parts orthogonal to direction vector
        double orthogonalFraction = vXRelative * -1 * dy + vYRelative * dx;
        double dphicalc = -Math.asin(orthogonalFraction * dt/d)/dt;
        double shearDamping = 0.05 * D * (first.L - dphicalc + second.L - dphicalc);
        orthogonalFraction -= shearDamping;
        vXRelative -=  orthogonalFraction * -1 * dy;
        vYRelative -=  orthogonalFraction * dx;
        System.err.format("orth vel.: %f, Phicalc: %f, orth Damping: %f\n", orthogonalFraction, dphicalc, shearDamping);

        first.fX -= c * vXRelative;
        first.fY -= c * vYRelative;
        second.fX += c * vXRelative;
        second.fY += c * vYRelative;
        first.T -= D * (first.L - second.L);
        second.T -= D * (second.L - first.L);
    }

    public double totalEnergy() {
        return energy;
    }
    public static void main(String[] args) {
        Cell[] cells = null;
        // Cell(Cell left, Cell right, Cell up, Cell down,
        //      double m, double I, double zeta, double omega0, double r, double E,
        //      int index, double x, double y, double theta,
        //      double vx, double vy, double L)
        if("beam oscillation rotation".equals(args[0])) {
            Cell a, b;
            bodyStatusReading = (body) -> ("" + body.totalEnergy());
            cellStatusReading = (c) -> {
                return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r   , E, index
                         10000, 10000, 1,  1, 0.5 , 1, 0,
                           //x, y,th,vx,vy, L
                             0, 0, PI/2, 0, 0, -0.0 * 2 * PI);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r, E, index
                             1, 1, 1, 1, 0.5, 1, 0,
                           //x, y,th,vx,vy, L
                             0.0, 1.0, PI/2, 0.1 * 2*PI, 0, -0.0 * 2*PI);
            cells = new Cell[] {a, b};
            a.right = b; b.left = a;
        } else if("beam oscillation".equals(args[0])) {
            Cell a, b;
            bodyStatusReading = (body) -> ("" + body.totalEnergy());
            cellStatusReading = (c) -> {
                return String.format("%f %f %f ", c.y, c.theta, c.L);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r   , E, index
                         10000, 10000, 0,  2*PI/10000, 0.5 , 1, 0,
                           //x, y,th,vx,vy, L
                             0, 0, 0, 0, 0, 0.00);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r, E, index
                             1, 1, 0, 2 * PI, 0.5, 1, 0,
                           //x, y,th,vx,vy, L
                             1.0, 0, 0, 0, 0.02, 0);
            cells = new Cell[] {a, b};
            a.right = b; b.left = a;
        } else if("orbit".equals(args[0])) {
            Cell a, b;
            cellStatusReading = (c) -> {
                return String.format(" %f %f ", c.x, c.y);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r   , E, index
                         10000, 1, 0,  2*PI/10000, 0.5 , 0, 0,
                           //x, y,th,vx,vy, L
                             0, 0, 0, 0, 0, 0);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r, E, index
                             1, 1, 0, 2 * PI, 0.5, 0, 0,
                           //x, y,th,vx,vy, L
                             1.0, 0, 0, 0, 1, 0);
            cells = new Cell[] {a, b};
            a.right = b; b.left = a;
        } else if("relative motion".equals(args[0])) {
            Cell a, b;
            a = new Cell(null, null, null, null,
                             1, 1, 0, 2*PI, 0.5, 1, 0,
                             -0.7, 0, 0, 1, 0, 0);
            b = new Cell(null, null, null, null,
                         1, 1, 0, 2*PI, 0.5, 1, 0,
                         +0.7, 0, 0, 1, 0, 0);
            cells = new Cell[] {a, b};
            a.right = b; b.left = a;
        } else if("basic rotation".equals(args[0])) {
            Cell a, b;
            cellStatusReading = (c) -> {
                return String.format(" %f %f ", c.x, c.y);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r  , E, index
                             1, 1, 1, 2*PI, 1, 0, 0,
                           //x, y,th,vx,vy, L
                            -2, 0, 0, 0, -2, 0);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r  , E, index
                             1, 1, 1, 2*PI, 1, 0, 0,
                           //x, y,th,vx,vy, L
                             2, 0, 0, 0, 2, 0);
            cells = new Cell[] {a, b};
            a.right = b; b.left = a;
        } else if("relative motion damping".equals(args[0])) {
            Cell a, b;
            a = new Cell(null, null, null, null,
                             1, 1, 0.1, 1, 1, 1, 0,
                             -2, 0, 0, 1, 0, 0);
            b = new Cell(null, null, null, null,
                         1, 1, 0.1, 1, 1, 1, 0,
                         2, 0, 0, 1, 0, 0);
            cells = new Cell[] {a, b};
            a.right = b; b.left = a;
        }

        SoftBody bod = new SoftBody(cells);
        dt = Double.parseDouble(args[1]);
        int nSteps = Integer.parseInt(args[2]);
        double t = 0;
        for(int i = 0; i < nSteps; i++) {
            System.out.format("%f ", t);
            bod.updateForces();
            bod.propagateCells();
            t += dt;
        }
        System.err.println("Number of Calls " + counter);
    }
}
