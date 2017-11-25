import physics.*;
import java.util.*;
import org.junit.Test;
import java.io.*;
import static util.PrettyPrint.printGrid;
import java.nio.file.*;
import static org.junit.Assert.*;
import static physics.SoftBody.*;
import static java.lang.Math.PI;
import static physics.Cell.dt;
import static util.Math.*;

public class TwoCellTests
{
    // static double[] energies;
    static int stepCounter=0;
    static double t;
    int i=0, j=0;

    static String makeHeader(int nCells) {
        String[] observableNames = {
            "X", "Y", "Theta",
            "VX", "VY", "L",
            "FX", "FY", "T"
        };
        String header = "";
        for(int i = 0; i < nCells; i++) {
            for(String s: observableNames)
                header += s + i + " ";
        }
        header += "Time Energy\n";
        return header;
    }

    @Test
    public void dampedRelativeRotation() throws IOException {
        Cell[] cells = null;
        Bond[] bonds = null;
        int nSteps = 400;
        int nCells = 2;
        dt = 0.01;
        int nCellParams = 9;
        double[][] simulationResults = new double[nSteps][nCells * nCellParams + 2];
        Cell a, b;
        SoftBody.cellStatusCallback = (c) -> {
             simulationResults[i][0+9*(j%nCells)] = c.x;
             simulationResults[i][1+9*(j%nCells)] = c.y;
             simulationResults[i][2+9*(j%nCells)] = c.theta;
             simulationResults[i][3+9*(j%nCells)] = c.vx;
             simulationResults[i][4+9*(j%nCells)] = c.vy;
             simulationResults[i][5+9*(j%nCells)] = c.L;
             simulationResults[i][6+9*(j%nCells)] = c.fX;
             simulationResults[i][7+9*(j%nCells)] = c.fY;
             simulationResults[i][8+9*(j%nCells)] = c.T;
             j++;
        };
        SoftBody.bodyStatusCallback = (bod) -> {
             simulationResults[i][nCellParams * nCells + 0] = t;
             simulationResults[i][nCellParams * nCells + 1] = bod.totalEnergy();
            i++;
            t += dt;
        };
        a = new Cell(
                       //m, I, Z, om0 , r   , E,
                       1, 1, 0, 1, 0.5 , 1,
                       //x, y,vx,vy, L
                       -0.5, 0, 0, -1, 0);
        b = new Cell(
                       //m, I, Z, om0 , r, E,
                         1, 1, 0, 1, 0.5, .1,
                       //x, y,vx,vy, L
                         0.5, 0, 0.0, 1, 0);
        cells = new Cell[] {a, b};
        bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
        testSimulation(cells, bonds, nSteps);
        double[] energies = getColumn(simulationResults, 19);
        String header = makeHeader(nCells);
        Files.write(Paths.get("plot.dat"),
            (header + printGrid(simulationResults)).getBytes()
        );
        double e0 = energies[0];
        double eav = average(energies);
        assertTrue(String.format("Energy is not conserved! E0= %f, Eaverage=%f", e0, eav)
        , eav <= e0);
    }

    // public static void main(String[] args) {
    //         Cell[] cells = null;
    //         Bond[] bonds = null;
    //         if("damped relative rotation".equals(args[0])) {
    //             Cell a, b;
    //             SoftBody.bodyStatusReading = (body) -> ("" + body.totalEnergy() + " " + phi0);
    //             cellStatusReading = (c) -> { return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
    //             a = new Cell(
    //                            //m, I, Z, om0 , r   , E,
    //                            1, 1, 1, 1, 0.5 , .1,
    //                            //x, y,vx,vy, L
    //                            -0.5, 0, 0, 0, 1);
    //             b = new Cell(
    //                            //m, I, Z, om0 , r, E,
    //                              1, 1, 1, 1, 0.5, .1,
    //                            //x, y,vx,vy, L
    //                              0.5, 0, 0.0, 0, -0.5);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = Double.parseDouble(args[1]);
    //             int nSteps = Integer.parseInt(args[2]);
    //             testSimulation(cells, bonds, nSteps);
    //         } else if("damped symmetric rotation oscillation".equals(args[0])) {
    //             Cell a, b;
    //             SoftBody.bodyStatusReading = (body) -> ("" + body.totalEnergy() + " " + phi0);
    //             cellStatusReading = (c) -> { return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
    //             a = new Cell(
    //                            //m, I, Z, om0 , r   , E,
    //                            1, 1, 1, 1, 0.5 , .1,
    //                            //x, y,vx,vy, L
    //                            -0.5, 0, 0, -1, 0);
    //             b = new Cell(
    //                            //m, I, Z, om0 , r, E,
    //                              1, 1, 1, 1, 0.5, .1,
    //                            //x, y,vx,vy, L
    //                              0.5, 0, 0.0, 1, 0);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = Double.parseDouble(args[1]);
    //             int nSteps = Integer.parseInt(args[2]);
    //             testSimulation(cells, bonds, nSteps);
    //         } else if("symmetric rotation oscillation".equals(args[0])) {
    //             Cell a, b;
    //             SoftBody.bodyStatusReading = (body) -> ("" + body.totalEnergy() + " " + phi0);
    //             cellStatusReading = (c) -> { return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
    //             a = new Cell(
    //                            //m, I, Z, om0 , r   , E,
    //                            1, 1, 0, 1, 0.5 , 1,
    //                            //x, y,vx,vy, L
    //                            -0.5, 0, 0, -1, -0.0 * 2 * PI);
    //             b = new Cell(
    //                            //m, I, Z, om0 , r, E,
    //                              1, 1, 0, 1, 0.5, 0.1,
    //                            //x, y,vx,vy, L
    //                              0.5, 0, 0.0, 1, 0.0 * 2*PI);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = Double.parseDouble(args[1]);
    //             int nSteps = Integer.parseInt(args[2]);
    //             testSimulation(cells, bonds, nSteps);
    //         } else if("beam oscillation rotation".equals(args[0])) {
    //             Cell a, b;
    //             bodyStatusReading = (body) -> ("" + body.totalEnergy());
    //             cellStatusReading = (c) -> {
    //                 return String.format("%f %f %f %f %f", c.x, c.y, c.theta, c.fX, c.fY);};
    //             a = new Cell(
    //                            //m, I, Z, om0 , r   , E, index
    //                            10000, 10000, 1, 1, 0.5 , 1,
    //                            //x, y,vx,vy, L
    //                            0, 0, 0, 0, -0.0 * 2 * PI);
    //             b = new Cell(
    //                            //m, I, Z, om0 , r, E,
    //                              1, 1, 1, 1, 0.5, 1,
    //                            //x, y,vx,vy, L
    //                              0, -1, 0.3, 0.0, 0.0 * 2*PI);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = Double.parseDouble(args[1]);
    //             int nSteps = Integer.parseInt(args[2]);
    //             testSimulation(cells, bonds, nSteps);
    //         } else if("beam oscillation".equals(args[0])) {
    //             Cell a, b;
    //             bodyStatusReading = (body) -> ("" + body.totalEnergy());
    //             cellStatusReading = (c) -> {
    //                 return String.format("%f %f %f ", c.y, c.theta, c.L);};
    //             a = new Cell(
    //                            //m, I, Z, om0 , r   , E, index
    //                          10000, 10000, 0,  2*PI/10000, 0.5 , 1,
    //                            //x, y,vx,vy, L
    //                              0, 0, 0, 0, 0.00);
    //             b = new Cell(
    //                            //m, I, Z, om0 , r, E, index
    //                              1, 1, 0, 2 * PI, 0.5, 1,
    //                            //x, y,vx,vy, L
    //                              1.0, 0, 0, 0.02, 0);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = Double.parseDouble(args[1]);
    //             int nSteps = Integer.parseInt(args[2]);
    //             testSimulation(cells, bonds, nSteps);
    //         } else if("orbit".equals(args[0])) {
    //             Cell a, b;
    //             SoftBody.cellStatusReading = (c) -> {
    //                 return String.format(" %f %f ", c.x, c.y);};
    //             a = new Cell(
    //                            //m, I, Z, om0 , r   , E
    //                          10000, 1, 0,  2*PI, 0.5 , 0,
    //                            //x, y,th,vx,vy, L
    //                              0, 0, 0, 0, 0);
    //             b = new Cell(
    //                            //m, I, Z, om0 , r, E
    //                              1, 1, 0, 2 * PI, 0.5, 0,
    //                            //x, y,th,vx,vy, L
    //                              1.0, 0, 0, .1, 0);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = Double.parseDouble(args[1]);
    //             int nSteps = Integer.parseInt(args[2]);
    //             testSimulation(cells, bonds, nSteps);
    //         } else if("relative motion".equals(args[0])) {
    //             Cell a, b;
    //             a = new Cell(
    //                              1, 1, 0, 2*PI, 0.5, 1,
    //                              -0.7, 0, 1, 0, 0);
    //             b = new Cell(
    //                          1, 1, 0, 2*PI, 0.5, 1,
    //                          +0.7, 0, 1, 0, 0);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = Double.parseDouble(args[1]);
    //             int nSteps = Integer.parseInt(args[2]);
    //             testSimulation(cells, bonds, nSteps);
    //         } else if("basic rotation".equals(args[0])) {
    //             Cell a, b;
    //             cellStatusReading = (c) -> {
    //                 return String.format(" %f %f ", c.x, c.y);};
    //             a = new Cell(
    //                            //m, I, Z, om0 , r  , E, index
    //                              1, 1, 1, 2*PI, 1, 0,
    //                            //x, y,th,vx,vy, L
    //                             -2, 0, 0, -2, 0);
    //             b = new Cell(
    //                            //m, I, Z, om0 , r  , E, index
    //                              1, 1, 1, 2*PI, 1, 0,
    //                            //x, y,th,vx,vy, L
    //                              2, 0, 0, 2, 0);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = Double.parseDouble(args[1]);
    //             int nSteps = Integer.parseInt(args[2]);
    //             testSimulation(cells, bonds, nSteps);
    //         } else if("relative motion damping".equals(args[0])) {
    //             Cell a, b;
    //             a = new Cell(
    //                              1, 1, 0.1, 1, 1, 1,
    //                              -2, 0, 1, 0, 0);
    //             b = new Cell(
    //                          1, 1, 0.1, 1, 1, 1,
    //                          2, 0, 1, 0, 0);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = Double.parseDouble(args[1]);
    //             int nSteps = Integer.parseInt(args[2]);
    //             testSimulation(cells, bonds, nSteps);
    //         } else if("minimize".equals(args[0])) {
    //             Cell a, b;
    //             a = new Cell(1, 1, 0, 1, 0.5 , 1, -0.5, 0, 0, -1, 0);
    //             b = new Cell(1, 1, 0, 1, 0.5, 0.1, 0.5, 0, 0.0, 1, 0);
    //             cells = new Cell[] {a, b};
    //             bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
    //             dt = 0.1;
    //             int nSteps = 900;
    //             double[] energies = new double[nSteps];
    //             ka = Double.parseDouble(args[1]);
    //             kb = Double.parseDouble(args[2]);
    //             kc = Double.parseDouble(args[3]);
    //             kd = Double.parseDouble(args[4]);
    //             stepCounter = -1;
    //             bodyStatusReading = (body) -> {
    //                 stepCounter++;
    //                 energies[stepCounter] = body.totalEnergy();
    //                 return "" + energies[stepCounter];
    //             };
    //             cellStatusReading = (c) -> { return "";};
    //             testSimulation(cells, bonds, nSteps);
    //             System.err.println("Run finished with " + Arrays.toString(energies));
    //             System.out.println(util.Math.standardDeviation(energies));
    //         }
    //     }

    static void testSimulation(Cell[] cells, Bond[] bonds, int nSteps) {
        SoftBody bod = new SoftBody(cells, bonds);
        t = 0;
        for(int i = 0; i < nSteps; i++) {
            bod.updateForces();
            bod.propagateCells();
        }
    }
}
