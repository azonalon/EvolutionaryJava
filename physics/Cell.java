package physics;
import static physics.SoftBody.dt;

class Cell {
    Cell left, right, up, down;

    Cell(Cell left, Cell right, Cell up, Cell down,
         double m, double I, double zeta, double omega0, double r, double E,
         int index, double x, double y, double theta, double vx, double vy,
         double L) {
        this.m = m;
        this.I = I;
        this.zeta = zeta;
        this.r = r;
        this.omega0 = omega0;
        this.k = m * omega0 * omega0/2;
        this.down = down;
        this.left = left;
        this.right = right;
        this.up = up;
        this.E = E;
        this.index = index;
        this.c = 2 * zeta * omega0;
        this.D = 2 * zeta * Math.sqrt(16 * E * r * r * r);
        this.x = x;
        this.y = y;
        this.L = L;
        this.theta = theta;
        this.vx = vx;
        this.vy = vy;
        assert(this.m > 0);
        assert(this.I > 0);
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

    double x, y, vx, vy, vx0, vy0, x0, y0, L, L0, theta, theta0, m, I, zeta, omega0, r,
           fX, fY, T, E, k, c, D;
    int index;
}
