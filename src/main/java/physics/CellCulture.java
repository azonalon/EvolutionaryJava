package physics;
import physics.*;

public enum CellCulture {
    STIFF(1, 1, 1, 1, 10, 10),
    SOFT(1, 1, 1, 1, .1, .1),
    NORMAL(1, 1, 1, 1, 1, 1),
    HEAVY(10, 10, 1, 1, 1, 1);
    public Cell grow() {
        return new Cell(r, m, I, zeta, omega0, E);
    }
    CellCulture(double r, double m, double I, double zeta, double omega0, double E) {
        this.r = r;
        this.m = m;
        this.I = I;
        this.zeta = zeta;
        this.omega0 = omega0;
        this.E = E;
    }
    double r, m , I, zeta, omega0, E;
}
