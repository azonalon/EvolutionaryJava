package physics.test;
import physics.*;
import java.util.*;
import static physics.SoftBody.*;
import static java.lang.Math.PI;

public class TwoCellTests
{
    static double[] energies;
    static int stepCounter=0;

    public static void main(String[] args) {
            Cell[] cells = null;
            Bond[] bonds = null;
            if("damped relative rotation".equals(args[0])) {
                Cell a, b;
                SoftBody.bodyStatusReading = (body) -> ("" + body.totalEnergy() + " " + phi0);
                cellStatusReading = (c) -> { return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
                a = new Cell(
                               //m, I, Z, om0 , r   , E,
                               1, 1, 1, 1, 0.5 , .1,
                               //x, y,vx,vy, L
                               -0.5, 0, 0, 0, 1);
                b = new Cell(
                               //m, I, Z, om0 , r, E,
                                 1, 1, 1, 1, 0.5, .1,
                               //x, y,vx,vy, L
                                 0.5, 0, 0.0, 0, -0.5);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = Double.parseDouble(args[1]);
                int nSteps = Integer.parseInt(args[2]);
                testSimulation(cells, bonds, nSteps);
            } else if("damped symmetric rotation oscillation".equals(args[0])) {
                Cell a, b;
                SoftBody.bodyStatusReading = (body) -> ("" + body.totalEnergy() + " " + phi0);
                cellStatusReading = (c) -> { return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
                a = new Cell(
                               //m, I, Z, om0 , r   , E,
                               1, 1, 1, 1, 0.5 , .1,
                               //x, y,vx,vy, L
                               -0.5, 0, 0, -1, 0);
                b = new Cell(
                               //m, I, Z, om0 , r, E,
                                 1, 1, 1, 1, 0.5, .1,
                               //x, y,vx,vy, L
                                 0.5, 0, 0.0, 1, 0);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = Double.parseDouble(args[1]);
                int nSteps = Integer.parseInt(args[2]);
                testSimulation(cells, bonds, nSteps);
            } else if("symmetric rotation oscillation".equals(args[0])) {
                Cell a, b;
                SoftBody.bodyStatusReading = (body) -> ("" + body.totalEnergy() + " " + phi0);
                cellStatusReading = (c) -> { return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
                a = new Cell(
                               //m, I, Z, om0 , r   , E,
                               1, 1, 0, 1, 0.5 , 1,
                               //x, y,vx,vy, L
                               -0.5, 0, 0, -1, -0.0 * 2 * PI);
                b = new Cell(
                               //m, I, Z, om0 , r, E,
                                 1, 1, 0, 1, 0.5, 0.1,
                               //x, y,vx,vy, L
                                 0.5, 0, 0.0, 1, 0.0 * 2*PI);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = Double.parseDouble(args[1]);
                int nSteps = Integer.parseInt(args[2]);
                testSimulation(cells, bonds, nSteps);
            } else if("beam oscillation rotation".equals(args[0])) {
                Cell a, b;
                bodyStatusReading = (body) -> ("" + body.totalEnergy());
                cellStatusReading = (c) -> {
                    return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
                a = new Cell(
                               //m, I, Z, om0 , r   , E, index
                               10000, 10000, 1, 1, 0.5 , 1,
                               //x, y,vx,vy, L
                               0, 0, 0, 0, -0.0 * 2 * PI);
                b = new Cell(
                               //m, I, Z, om0 , r, E,
                                 1, 1, 1, 1, 0.5, 1,
                               //x, y,vx,vy, L
                                 0, -1, 0.3, 0.0, 0.0 * 2*PI);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = Double.parseDouble(args[1]);
                int nSteps = Integer.parseInt(args[2]);
                testSimulation(cells, bonds, nSteps);
            } else if("beam oscillation".equals(args[0])) {
                Cell a, b;
                bodyStatusReading = (body) -> ("" + body.totalEnergy());
                cellStatusReading = (c) -> {
                    return String.format("%f %f %f ", c.y, c.theta, c.L);};
                a = new Cell(
                               //m, I, Z, om0 , r   , E, index
                             10000, 10000, 0,  2*PI/10000, 0.5 , 1,
                               //x, y,vx,vy, L
                                 0, 0, 0, 0, 0.00);
                b = new Cell(
                               //m, I, Z, om0 , r, E, index
                                 1, 1, 0, 2 * PI, 0.5, 1,
                               //x, y,vx,vy, L
                                 1.0, 0, 0, 0.02, 0);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = Double.parseDouble(args[1]);
                int nSteps = Integer.parseInt(args[2]);
                testSimulation(cells, bonds, nSteps);
            } else if("orbit".equals(args[0])) {
                Cell a, b;
                SoftBody.cellStatusReading = (c) -> {
                    return String.format(" %f %f ", c.x, c.y);};
                a = new Cell(
                               //m, I, Z, om0 , r   , E
                             10000, 1, 0,  2*PI, 0.5 , 0,
                               //x, y,th,vx,vy, L
                                 0, 0, 0, 0, 0);
                b = new Cell(
                               //m, I, Z, om0 , r, E
                                 1, 1, 0, 2 * PI, 0.5, 0,
                               //x, y,th,vx,vy, L
                                 1.0, 0, 0, .1, 0);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = Double.parseDouble(args[1]);
                int nSteps = Integer.parseInt(args[2]);
                testSimulation(cells, bonds, nSteps);
            } else if("relative motion".equals(args[0])) {
                Cell a, b;
                a = new Cell(
                                 1, 1, 0, 2*PI, 0.5, 1,
                                 -0.7, 0, 1, 0, 0);
                b = new Cell(
                             1, 1, 0, 2*PI, 0.5, 1,
                             +0.7, 0, 1, 0, 0);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = Double.parseDouble(args[1]);
                int nSteps = Integer.parseInt(args[2]);
                testSimulation(cells, bonds, nSteps);
            } else if("basic rotation".equals(args[0])) {
                Cell a, b;
                cellStatusReading = (c) -> {
                    return String.format(" %f %f ", c.x, c.y);};
                a = new Cell(
                               //m, I, Z, om0 , r  , E, index
                                 1, 1, 1, 2*PI, 1, 0,
                               //x, y,th,vx,vy, L
                                -2, 0, 0, -2, 0);
                b = new Cell(
                               //m, I, Z, om0 , r  , E, index
                                 1, 1, 1, 2*PI, 1, 0,
                               //x, y,th,vx,vy, L
                                 2, 0, 0, 2, 0);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = Double.parseDouble(args[1]);
                int nSteps = Integer.parseInt(args[2]);
                testSimulation(cells, bonds, nSteps);
            } else if("relative motion damping".equals(args[0])) {
                Cell a, b;
                a = new Cell(
                                 1, 1, 0.1, 1, 1, 1,
                                 -2, 0, 1, 0, 0);
                b = new Cell(
                             1, 1, 0.1, 1, 1, 1,
                             2, 0, 1, 0, 0);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = Double.parseDouble(args[1]);
                int nSteps = Integer.parseInt(args[2]);
                testSimulation(cells, bonds, nSteps);
            } else if("minimize".equals(args[0])) {
                Cell a, b;
                a = new Cell(1, 1, 0, 1, 0.5 , 1, -0.5, 0, 0, -1, 0);
                b = new Cell(1, 1, 0, 1, 0.5, 0.1, 0.5, 0, 0.0, 1, 0);
                cells = new Cell[] {a, b};
                bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
                dt = 0.1;
                int nSteps = 900;
                double[] energies = new double[nSteps];
                ka = Double.parseDouble(args[1]);
                kb = Double.parseDouble(args[2]);
                kc = Double.parseDouble(args[3]);
                kd = Double.parseDouble(args[4]);
                stepCounter = -1;
                bodyStatusReading = (body) -> {
                    stepCounter++;
                    energies[stepCounter] = body.totalEnergy();
                    return "" + energies[stepCounter];
                };
                cellStatusReading = (c) -> { return "";};
                testSimulation(cells, bonds, nSteps);
                System.err.println("Run finished with " + Arrays.toString(energies));
                System.out.println(util.Math.standardDeviation(energies));
            }
        }

    static void testSimulation(Cell[] cells, Bond[] bonds, int nSteps) {
        SoftBody bod = new SoftBody(cells, bonds);
        double t = 0;
        for(int i = 0; i < nSteps; i++) {
            System.out.format("%f ", t);
            bod.updateForces();
            bod.propagateCells();
            t += dt;
        }
    }
}
