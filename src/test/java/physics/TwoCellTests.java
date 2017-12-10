package physics;
import physics.*;
import org.junit.Test;
import org.junit.Rule;
import java.io.*;
import org.junit.rules.*;
import static util.PrettyPrint.printGrid;
import java.nio.file.*;
import static org.junit.Assert.*;
// import static physics.SoftBody.*;
import static java.lang.Math.PI;
// import static physics.Cell.dt;
import static util.Math.*;

public class TwoCellTests
{
    // static double[] energies;
    static int stepCounter=0;
    static double t;
    int i=0, j=0;

    @Rule
    public TestName name = new TestName();

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

    public void twoCellTestCase(
        double m1, double I1, double zeta1, double k1, double r1, double E1,
        double x1, double y1, double vx1, double vy1, double L1,
        double m2, double I2, double zeta2, double k2, double r2, double E2,
        double x2, double y2, double vx2, double vy2, double L2,
        double dt, int nSteps
    ) {
        Cell a, b;
        a = new Cell(
            m1, I1, zeta1, k1, r1, E1,
            x1, y1, vx1, vy1, L1);
        b = new Cell(
            m2, I2, zeta2, k2, r2, E2,
            x2, y2, vx2, vy2, L2
        );
        runTwoCellTestCase(a, b, dt, nSteps);
    }

    @Test
    public void dampedRelativeRotation (){
        double
        m1= 1.0, I1= 1.0, zeta1= 0.1, k1=  1.0, r1= 0.5, E1= 1.0,
        x1=-0.5, y1= 0.0, vx1=   0.0, vy1= 0.0, L1=-1.0,
        m2= 1.0, I2= 1.0, zeta2= 0.1, k2=  1.0, r2= 0.5, E2= 1.0,
        x2= 0.5, y2= 0.0, vx2=   0.0, vy2= 0.0, L2=+1.0,
        dt= 0.1;
        int nSteps = 100;

        twoCellTestCase(
            m1, I1, zeta1, k1, r1, E1,
            x1, y1, vx1, vy1, L1,
            m2, I2, zeta2, k2, r2, E2,
            x2, y2, vx2, vy2, L2,
            dt, nSteps
        );
    }
    @Test
    public void testotesto (){
        double
        m1= 1.0, I1= 0.1, zeta1= 1.0, k1=  1.0, r1= 0.5, E1= 0.1,
        x1=-0.5, y1= 0.0, vx1=   0.0, vy1=-1.0, L1=+0.0,
        m2= 1.0, I2= 0.1, zeta2= 1.0, k2=  1.0, r2= 0.5, E2= 0.1,
        x2= 0.5, y2= 0.0, vx2=   0.0, vy2= 1.0, L2=+0.0,
        dt= 0.1;
        int nSteps = 100;

        twoCellTestCase(
            m1, I1, zeta1, k1, r1, E1,
            x1, y1, vx1, vy1, L1,
            m2, I2, zeta2, k2, r2, E2,
            x2, y2, vx2, vy2, L2,
            dt, nSteps
        );
    }

    @Test
    public void dampedAbsoluteRotation (){
        double
        m1= 1.0, I1= 1.0, zeta1= 1.0, k1=  1.0, r1= 0.5, E1= 1.0,
        x1=-0.5, y1= 0.0, vx1=   0.0, vy1= 0.0, L1=+1.0,
        m2= 1.0, I2= 1.0, zeta2= 1.0, k2=  1.0, r2= 0.5, E2= 1.0,
        x2= 0.5, y2= 0.0, vx2=   0.0, vy2= 0.0, L2=+1.0,
        dt= 0.001;
        int nSteps = 10000;

        twoCellTestCase(
            m1, I1, zeta1, k1, r1, E1,
            x1, y1, vx1, vy1, L1,
            m2, I2, zeta2, k2, r2, E2,
            x2, y2, vx2, vy2, L2,
            dt, nSteps
        );
    }

    public void runTwoCellTestCase(Cell a, Cell b, double dt, int nSteps) {
        Cell[] cells = null;
        Bond[] bonds = null;
        int nCells = 2;
        int nCellParams = 9;
        double[][] simulationResults = new double[nSteps][nCells * nCellParams + 2];
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
        cells = new Cell[] {a, b};
        bonds = new Bond[] {Bond.harmonicAverageBond(a, b, 0.0)};
        physics.Cell.dt = dt;
        testSimulation(cells, bonds, nSteps);
        double[] energies = getColumn(simulationResults, 19);
        String header = makeHeader(nCells);
        try {
            Path path =  Paths.get("build/test-results/physics/SoftBody/TwoCell");
            if(Files.notExists(path)){
                    Files.createDirectories(path);
                }
            Files.write(path.resolve((name.getMethodName() + ".dat")),
                (header + printGrid(simulationResults)).getBytes()
            );
        } catch(IOException e) {
            throw new RuntimeException("Could not write data");
        }
        double e0 = energies[0];
        double eav = average(energies);
        assertTrue(String.format("Energy is not conserved! E0= %f, Eaverage=%f", e0, eav)
        , eav <= e0);
    }

    @Test
    public void dampedRelativeRotation2() {
        twoCellTestCase(
        1, 1, 1, 1, 0.5 , .1,
        -0.5, 0, 0, 0, 1,
        1, 1, 1, 1, 0.5, .1,
        0.5, 0, 0.0, 0, -1.0,
        0.1, 200
        );
    }
    @Test
    public void dampedSymmetricRotationOscillation() {
                                   twoCellTestCase(
                                   1, 1, 1, 1, 0.5 , .1,
                                   -0.5, 0, 0, -1, 0,
                                     1, 1, 1, 1, 0.5, .1,
                                     0.5, 0, 0.0, 1, 0,
                                     0.01, 5000
                                   );
    }

    @Test
    public void symmetricRotationOscillation() {
                                   twoCellTestCase(
                                   1, 1, 0, 1, 0.5 , 1,
                                   -0.5, 0, 0, -1, 0.0,
                                     1, 1, 0, 1, 0.5, 1,
                                     0.5, 0, 0.0, 1, 0.0,
                                     0.001, 10000
                                   );
    }
    @Test
    public void beamOscillationRotation() {
                                  twoCellTestCase (
                                   10000, 10000, 1, 1, 0.5 , 1,
                                   0, 0, 0, 0, -0.0 * 2 * PI,
                                     1, 1, 1, 1, 0.5, 1,
                                     0, -1, 0.3, 0.0, 0.0 * 2*PI,
                                     0.1, 100
                                   );
    }
    @Test
    public void beamOscillation() {
        twoCellTestCase(
        1000, 1000, 1, 1, 0.5 , .1,
        0, 0, 0, 0, 0,
        1, 1, 1, 1, 0.5, .1,
        1, 0, 0.0, 1, 0,
        0.1, 100
        );
    }
    @Test
    public void orbit() {
                                 twoCellTestCase(
                                 10000, 1, 0,  2*PI, 0.5 , 0,
                                     0, 0, 0, 0, 0,
                                     1, 1, 0, 2 * PI, 0.5, 0,
                                     1.0, 0, 0, 1, 0,
                                     0.1, 100
                                 );
    }
    @Test
    public void relativeMotion() {
                                     twoCellTestCase(
                                         1, 1, 0, 2*PI, 0.5, 1,
                                         -0.7, 0, 1, 0, 0,
                                     1, 1, 0, 2*PI, 0.5, 1,
                                     +0.7, 0, 1, 0, 0,
                                     0.1, 100
                                     );
    }
    @Test
    public void basicRotation() {
        twoCellTestCase(
        1, 1, 1, 2*PI, 1, 0,
        -2, 0, 0, -2, 0,
        1, 1, 1, 2*PI, 1, 0,
        2, 0, 0, 2, 0,
        0.1, 100
        );
    }

    static void testSimulation(Cell[] cells, Bond[] bonds, int nSteps) {
        SoftBody bod = new SoftBody(cells, bonds);
        t = 0;
        for(int i = 0; i < nSteps; i++) {
            bod.updateForces();
            bod.propagateCells();
        }
    }
}
