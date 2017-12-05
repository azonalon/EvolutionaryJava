package physics;
import physics.*;

public enum CellCulture {
    STIFF(0.5, 1, 1, 1, 10, 10),
    SOFT(0.5, 1, 1, 1, .1, .1),
    NORMAL(0.5, 1, 1, 1, 1, 1),
    OSCILLATOR(0.5, 1, 1, 0.05, 1, 1),
    HEAVY(0.5, 10, 1, 1, 1, 1);

    public Cell grow() {
        return new Cell(r, m, I, zeta, omega0, E);
    }
    static public Cell[][] growArray(CellCulture[][] types) {
        Cell[][] cells = new Cell[types.length][];
        for(int i = 0; i < types.length; i++) {
            cells[i] = new Cell[types[i].length];
            for(int j = 0; j < types[i].length; j++) {
                cells[i][j] = types[i][j].grow();
            }
        }
        return cells;
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
