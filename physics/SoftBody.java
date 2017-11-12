package physics;

import java.util.*;

public class SoftBody {
    static class Cell {
        List <Cell> connectedCells;
        double x, y, px, py, L, m, zeta, omega0, r;
    }
}
