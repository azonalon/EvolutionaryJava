package physics;
import physics.*;
import java.util.*;
import java.util.function.*;
import static java.lang.Math.PI;
import static java.lang.Math.sqrt;
// import static System.out.format;

public class SoftBody {
    Cell[] cells;
    int [] cellStates; // these are used in BFS 'update forces' to check
                       // if bond force was calculated
    static int LEFT = 1 << 0;
    static int RIGHT= 1 << 1;
    static int UP   = 1 << 2;
    static int DOWN = 1 << 3;
    static int counter=0;
    static double dt = 0.1;
    static Consumer<Cell> printCellCoordinates = (Cell c)->System.out.format("%f ", c.x);

    SoftBody(Cell[] cells) {
        this.cells = cells;
        cellStates = new int[cells.length];
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
            printCellCoordinates.accept(c);
            c.propagate();
        }
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
        double averageK = (first.k + second.k)/2;

        double dx  = - first.x + second.x;
        double dy  = - first.y + second.y;
        double dLength = Math.sqrt(dx*dx + dy*dy);
        double phi = Math.asin(dy/dLength);
        // double theta1 = Math.tan((first.theta - phi)%PI/2.0) ;
        // double theta2 = Math.tan((second.theta - phi)%PI/2.0);
        double theta1 = (first.theta - phi);
        double theta2 = (second.theta - phi);

        double tTheta = theta1 + theta2;
        // rest length of the spring
        double l = first.r + second.r;
        double E = (first.E + second.E)/2;
        double D = (first.D + second.D)/2;
        double c = (first.c + second.c)/2;
        double fShear = - 6 * E * l * l * (tTheta);
        dx = dx / dLength;
        dy = dy / dLength;

        // System.out.format("\n averageK %f \n", averageK);

        System.err.format("%f %f %f\n", phi, first.theta, second.theta);
        first.fX += (dLength - l) * averageK * dx + fShear * (-1*dy);
        first.fY += (dLength - l) * averageK * dy + fShear * (   dx);
        second.fX -= (dLength - l) * averageK * dx + fShear * (-1*dy);
        second.fY -= (dLength - l) * averageK * dy + fShear * (   dx);

        // TODO:  damp relative rotation
        first.T  += 2 * E * l*l*l * (tTheta + theta1);
        second.T -= 2 * E * l*l*l * (tTheta + theta2);

        // damping forces
        double vXRelative = first.vx - second.vx;
        double vYRelative = first.vy - second.vy;
        // subtract parts orthogonal to direction vector
        double orthogonalFraction = vXRelative * -1 * dy + vYRelative * dx;
        vXRelative -= orthogonalFraction * -1 * dy;
        vYRelative -= orthogonalFraction * dx;

        first.fX -= c * vXRelative;
        first.fY -= c * vYRelative;
        second.fX += c * vXRelative;
        second.fY += c * vYRelative;
    }

    public static void main(String[] args) {
        Cell[] cells = null;
        // Cell(Cell left, Cell right, Cell up, Cell down,
        //      double m, double I, double zeta, double omega0, double r, double E,
        //      int index, double x, double y, double theta,
        //      double vx, double vy, double L)
        if("beam oscillation".equals(args[0])) {
            Cell a, b;
            printCellCoordinates = (c) -> {
                System.out.format(" %f %f ", c.x, c.y);};
            a = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r   , E, index
                         10000, 10000, 0,  2*PI/10000, 0.5 , 1, 0,
                           //x, y,th,vx,vy, L
                             0, 0, 0, 0, 0, 0);
            b = new Cell(null, null, null, null,
                           //m, I, Z, om0 , r, E, index
                             1, 1, 0, 2 * PI, 0.5, 1, 0,
                           //x, y,th,vx,vy, L
                             1.0, 0, 0, 0, 1, 0);
            cells = new Cell[] {a, b};
            a.right = b; b.left = a;
        } else if("orbit".equals(args[0])) {
            Cell a, b;
            printCellCoordinates = (c) -> {
                System.out.format(" %f %f ", c.x, c.y);};
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
            printCellCoordinates = (c) -> {System.out.format(" %f %f ", c.x, c.y);};
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
            System.out.println();
        }
        System.err.println("Number of Calls " + counter);
    }
}
