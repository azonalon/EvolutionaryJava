package physics;
import java.util.*;
import static java.lang.Math.PI;

public class SoftBody {
    Cell[] cells;
    int [] cellStates; // these are used in BFS 'update forces' to check
                       // if bond force was calculated
    static int LEFT = 1 << 0;
    static int RIGHT= 1 << 1;
    static int UP   = 1 << 2;
    static int DOWN = 1 << 3;
    static  double dt = 0.1;

    SoftBody() {
        cells = new Cell[] {
            new Cell(null, null, null, null,
                     1, 1, 1, 1, 1, 0,
                     0, 0, 0, 1, 1, 0)
        };
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
            c.propagate();
            System.out.format("%f ", c.x);
        }
    }

    void visit(Cell a, int side) {
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
        double tTheta = first.theta + second.theta;
        double dTheta = first.theta - second.theta;
        // rest length of the spring
        double l = first.r + second.r;
        double E = (first.E + second.E)/2;
        double D = (first.D + second.D)/2;
        double c = (first.c + second.c)/2;
        double dx  = first.x - second.x;
        double dy  = first.y - second.y;
        double dLength = Math.sqrt(dx*dx + dy*dy);
        double fShear = - 6 * E * l * l * (tTheta);
        dx /= dLength;
        dy /= dLength;

        first.fX += (dLength - l) * averageK * dx + fShear * (-1*dy);
        first.fY += (dLength - l) * averageK * dy + fShear * (   dx);
        second.fX -= first.fX;
        second.fY -= first.fY;

        first.T += 2 * E * l*l*l * (tTheta + first.theta);
        second.T -= 2 * E * l*l*l * (tTheta + second.theta);

        // damping forces
        double vXRelative = first.vx - second.vx;
        double vYRelative = first.vy - second.vy;

        first.fX -= c * vXRelative;
        first.fY -= c * vYRelative;
        second.fX += c * vXRelative;
        second.fY += c * vYRelative;

        first.T -= D * dTheta;
        second.T += D * dTheta;
    }

    static class Cell {
        Cell left, right, up, down;

        Cell(Cell left, Cell right, Cell up, Cell down,
             double m, double zeta, double omega0, double r, double E, int index,
             double x0, double y0, double theta0,
             double vx0, double vy0, double L0) {
            this.m = m;
            this.zeta = zeta;
            this.r = r;
            this.omega0 = omega0;
            this.k = m * omega0 * omega0;
            this.down = down;
            this.left = left;
            this.right = right;
            this.up = up;
            this.E = E;
            this.index = index;
            this.c = 2 * zeta * omega0;
            this.D = 2 * zeta * Math.sqrt(16 * E * r * r * r / I);
            this.x0 = x0;
            this.y0 = y0;
            this.L0 = L0;
            this.theta0 = theta0;
            this.vx0 = vx0;
            this.vy0 = vy0;
        }

        void propagate() {
            x0 = x;
            vx0 = vx;
            vx  = vx0 + dt * fX/m;
            x  = x0 + dt * vx;

            y0 = y;
            vy0 = vy;
            vy  = vy0 + dt * fY/m;
            y  = y0 + dt * vy;

            theta0 = theta;
            L0 = L;
            L  = L0 + dt * T/I;
            theta  = theta0 + dt * L;
        };

        /**
         *
         */
        double x, y, vx, vy, vx0, vy0, x0, y0, L, L0, theta, theta0, m, I, zeta, omega0, r,
               fX, fY, T, E, k, c, D;
        int index;
    }
    public static void main(String[] args) {
        SoftBody b = new SoftBody();
        double t = 0;
        for(int i = 0; i < 10; i++) {
            System.out.format("%f ", t);
            b.propagateCells();
            t += dt;
            System.out.println();
        }
    }
}
