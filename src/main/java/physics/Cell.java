package physics;


public class Cell {
    public static double dt = 0.1;
    static double t=0;
    public Cell(
         double m, double I, double zeta, double omega0, double r, double E,
         double x, double y, double vx, double vy,
         double L) {
        this.m = m;
        this.I = I;
        this.zeta = zeta;
        this.r = r;
        this.omega0 = omega0;
        this.k = omega0;
        this.E = E;
        this.c = 2 * zeta * Math.sqrt(m * k);
        this.D = 2 * zeta * Math.sqrt(32*r*r*r*E * I);
        this.x = x;
        this.y = y;
        this.L = L;
        this.vx = vx;
        this.vy = vy;
        assert(this.m > 0);
        assert(this.I > 0);
    }
    public void linkToBody(SoftBody body) {
        this.body = body;
    }

    public double getX() {
        return body.X0.get(this.index + 0);
    }
    public double getY() {
        return body.X0.get(this.index + 1);
    }
    public double getAngle() {
        return body.X0.get(this.index + 2);
    }

    public double getVX() {
        return body.v.get(this.index + 0);
    }
    public double getVY() {
        return body.v.get(this.index + 1);
    }
    public double getL() {
        return body.v.get(this.index + 2);
    }
    public double getFX() {
        return body.f.get(this.index + 0);
    }
    public double getFY() {
        return body.f.get(this.index + 1);
    }
    public double getT() {
        return body.f.get(this.index + 2);
    }

    public void setPosition(double x, double y) {
        body.X0.set(this.index + 0, x);
        body.X0.set(this.index + 1, x);
    }
    public void setAngle(double theta) {
        body.X0.set(this.index + 2, x);
    }

    public Cell(double r, double m, double I, double zeta, double omega0, double E) {
        this(m, I, zeta, omega0, r, E, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    public String toString() {
        return String.format("Cell(x=%4.4f, y=%4.4f, vx=%4.4f, vy=%4.4f, th=%4.4f, L=%4.4f\n",
                x, y, vx, vy, theta, L) +
               String.format("m=%4.4f, I=%4.4f, D=%4.4f, c=%4.4f, k=%4.4f, index=%d)\n",
               m, I, D, c, k, index);
    }

    public SoftBody body=null;
    public double x, y, vx, vy, L, theta, m=1, I=1, zeta, omega0, r=1, fX, fY, T, E, k, c, D;
    public int index;
}
