package physics;

public class Cell {
    public static double dt = 0.1;
    static double t=0;
    public Cell(
         double m, double I, double zeta, double k, double r, double E,
         double x, double y, double vx, double vy,
         double L) {
        this(m, I, zeta, zeta, k, E, r);
        this.x = x;
        this.y = y;
        this.L = L;
        this.vx = vx;
        this.vy = vy;
     };
    public Cell(
         double m, double I,
         double zeta, double alpha,
         double k, double E, double r) {
        this.m = m;
        this.I = I;
        this.zeta = zeta;
        this.alpha = zeta;
        this.r = r;
        this.k = k;
        this.E = E;
        assert(this.m > 0);
        assert(this.I > 0);
    }
    /**
     * Deprecated
     */
    public Cell(double r, double m, double I, double zeta, double k, double E) {
        this(m, I, zeta, k, r, E, 0.0, 0.0, 0.0, 0.0, 0.0);
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
        return body.V0.get(this.index + 0);
    }
    public double getVY() {
        return body.V0.get(this.index + 1);
    }
    public double getL() {
        return body.V0.get(this.index + 2);
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
    public void setVelocity(double vx, double vy) {
        body.X1.set(this.index + 0, body.X0.get(this.index + 0) - dt*vx);
        body.X1.set(this.index + 1, body.X0.get(this.index + 1) - dt*vy);
    }
    public void setAngle(double theta) {
        body.X0.set(this.index + 2, x);
    }


    public String toString() {
        return String.format("Cell(x=%4.4f, y=%4.4f, vx=%4.4f, vy=%4.4f, th=%4.4f, L=%4.4f\n",
                x, y, vx, vy, theta, L) +
               String.format("m=%4.4f, I=%4.4f, zeta=%4.4f, alpha=%4.4f, k=%4.4f, index=%d)\n",
               m, I, zeta, alpha, k, index);
    }

    public SoftBody body=null;
    public double x, y, vx, vy, L, theta, m=1, I=1,
                  zeta, // translational damping coefficient
                  alpha,// shear damping coefficient
                  r=1, fX, fY, T, E, k;
    public int index;
}
