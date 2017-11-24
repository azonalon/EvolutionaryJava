package physics;
import static physics.SoftBody.dt;

public class Cell {
    Cell left, right, up, down;

    public Cell(Cell left, Cell right, Cell up, Cell down,
         double m, double I, double zeta, double omega0, double r, double E,
         double x, double y, double vx, double vy,
         double L) {
        this.m = m;
        this.I = I;
        this.zeta = zeta;
        this.r = r;
        this.omega0 = omega0;
        this.k = omega0;
        this.down = down;
        this.left = left;
        this.right = right;
        this.up = up;
        this.E = E;
        this.c = 2 * zeta * Math.sqrt(m * k);
        this.D = 2 * zeta * Math.sqrt(E * I);
        this.x = x;
        this.y = y;
        this.L = L;
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

    public String toString() {
        return String.format("x=%4.4f, y=%4.4f, vx=%4.4f, vy=%4.4f, th=%4.4f, L=%4.4f\n",
                x, y, vx, vy, theta, L) +
               String.format("m=%4.4f, I=%4.4f, D=%4.4f, c=%4.4f, k=%4.4f, index=%d\n",
               m, I, D, c, k, index);
    }

    public double x, y, vx, vy, vx0, vy0, x0, y0, L, L0, theta, theta0, m=1, I=1, zeta, omega0, r=1,
           fX=0, fY=0, T, E, k, c, D;
    public int index;
}
