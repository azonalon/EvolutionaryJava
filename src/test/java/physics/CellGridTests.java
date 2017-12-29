package physics;
import physics.*;
import java.util.*;
import org.junit.Test;
import org.junit.Rule;
import java.nio.file.*;
import java.io.*;
import org.junit.rules.*;
import static util.PrettyPrint.printGrid;
import java.nio.file.*;
import static org.junit.Assert.*;
import static physics.Cell.dt;
import static physics.Cell.t;
import static util.Math.*;
import static physics.CellCulture.*;
import java.util.function.*;

public class CellGridTests
{
    static Consumer<SoftBody> stepFunction;
    static Consumer<SoftBody> implicit = (bod) -> bod.implicitEulerStep();
    // static double[] energies;
    static int stepCounter=0;
    static int nSteps = 3;
    static double dt = 0.1;
    Cell[][] grid;
    SoftBody body=null;
    int i=0, j=0;

    @Rule
    public TestName name = new TestName();

    static String makeHeader(int nCells) {
        String[] observableNames = {
            "     X","     Y"," Theta",
            "    VX","    VY","     L",
            "    FX","    FY","     T"
        };
        String header = "";
        for(int i = 0; i < nCells; i++) {
            for(String s: observableNames)
                header += s + i + " ";
        }
        header += "Time    Energy \n";
        return header;
    }
    @Test
    public void fourOnFourStretching() {
        CellCulture[][] layout = new CellCulture[][] {
            {NORMAL, NORMAL},
            {NORMAL, NORMAL}
        };
        grid = CellCulture.growArray(layout);
        double cellWidth = 1;
        double cellHeight = 1;
        body = new SoftBody(grid, cellWidth, cellHeight);
        body.cellForceCallback = (cell) -> {
            if(cell == grid[0][0])
            cell.fY += +2;
            if(cell == grid[1][0])
            cell.fY += -2;
        };
        executeCellGridTestCase(dt, nSteps);
    }

    @Test
    public void twoOnTwoSimpleMovement() {
        CellCulture[][] layout = new CellCulture[][] {
            {NORMAL, NORMAL},
            {NORMAL, NORMAL}
        };
        grid = CellCulture.growArray(layout);
        double cellWidth = 1;
        double cellHeight = 1;
        body = new SoftBody(grid, cellWidth, cellHeight);
        grid[0][0].setVelocity(1, 0.5);
        grid[1][1].setVelocity(-1, 0.8);
        executeCellGridTestCase(dt, nSteps);
    }

    @Test
    public void gridWaves() {
        CellCulture[][] layout = new CellCulture[][] {
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
        };
        grid = CellCulture.growArray(layout);
        double cellWidth = 1;
        double cellHeight = 1;
        body = new SoftBody(grid, cellWidth, cellHeight);
        body.cellForceCallback = (cell) -> {
            if(cell == grid[0][0]) {
                if(t < 0.03/0.07) {
                        cell.fX += 2*Math.sin(0.03 * 2*Math.PI/0.07*t);
                        cell.fY -= 2*Math.sin(0.03 * 2*Math.PI/0.07*t);
                    }
                }
            // if(cell == grid[4][0])
            //     cell.fX += 2*Math.sin(0.01 * 2*Math.PI/0.07*t);
        };
        executeCellGridTestCase(dt, nSteps);
    }

    @Test
    public void shockForce() {
        CellCulture[][] layout = new CellCulture[][] {
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
        };
        grid = CellCulture.growArray(layout);
        double cellWidth = 1;
        double cellHeight = 1;
        body = new SoftBody(grid, cellWidth, cellHeight);
        body.cellForceCallback = (cell) -> {
            if(cell == grid[0][0]) {
                    cell.fX += 1*Math.sin(0.01 * 2.0*Math.PI/dt*t);
                }
            if(cell == grid[4][0]) {
                    cell.fX += 1*Math.sin(0.01 * 2.0*Math.PI/dt*t);
                }
            // if(cell == grid[4][0])
            //     cell.fX += 2*Math.sin(0.01 * 2*Math.PI/0.07*t);
        };
        executeCellGridTestCase(dt, nSteps);
    }

    @Test
    public void slurp() {
        CellCulture[][] layout = new CellCulture[][] {
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
            {NORMAL, NORMAL, NORMAL, NORMAL, NORMAL},
        };
        grid = CellCulture.growArray(layout);
        double cellWidth = 1;
        double cellHeight = 1;
        body = new SoftBody(grid, cellWidth, cellHeight);
        body.cellForceCallback = (cell) -> {
            if(cell == grid[0][0]) {
                    cell.fX += 1*Math.sin(0.01 * 2*Math.PI/0.07*t);
                }
            if(cell == grid[4][0]) {
                    cell.fX += 1*Math.sin(0.01 * 2*Math.PI/0.07*t);
                }
            // if(cell == grid[4][0])
            //     cell.fX += 2*Math.sin(0.01 * 2*Math.PI/0.07*t);
        };
        executeCellGridTestCase(dt, nSteps);
    }

    @Test
    public void fourOnFourEquilibrium() {
        CellCulture[][] layout = new CellCulture[][] {
            {NORMAL, NORMAL},
            {NORMAL, NORMAL}
        };
        grid = CellCulture.growArray(layout);
        double cellWidth = 1;
        double cellHeight = 1;
        body = new SoftBody(grid, cellWidth, cellHeight);
        body.cellForceCallback = (cell) -> {};
        executeCellGridTestCase(dt, nSteps);
    }

    @Test
    public void fourOnFourExpanded() {
        CellCulture[][] layout = new CellCulture[][] {
            {NORMAL, NORMAL},
            {NORMAL, NORMAL}
        };
        grid = CellCulture.growArray(layout);
        double cellWidth = 2;
        double cellHeight = 2;
        body = new SoftBody(grid, cellWidth, cellHeight);
        body.cellForceCallback = (cell) -> {};
        executeCellGridTestCase(dt, nSteps);
    }


    public void executeCellGridTestCase(double dt, int nSteps) {
        int nCells = grid.length*grid[0].length;
        int nCellParams = 9;
        double[][] simulationResults = new double[nSteps][nCells * nCellParams + 2];
        SoftBody.cellStatusCallback = (c) -> {
             simulationResults[i][0+9*(j%nCells)] = c.getX();
             simulationResults[i][1+9*(j%nCells)] = c.getY();
             simulationResults[i][2+9*(j%nCells)] = c.getAngle();
             simulationResults[i][3+9*(j%nCells)] = c.getVX();
             simulationResults[i][4+9*(j%nCells)] = c.getVY();
             simulationResults[i][5+9*(j%nCells)] = c.getL();
             simulationResults[i][6+9*(j%nCells)] = c.getFX();
             simulationResults[i][7+9*(j%nCells)] = c.getFY();
             simulationResults[i][8+9*(j%nCells)] = c.getT();
             j++;
        };
        SoftBody.bodyStatusCallback = (bod) -> {
             simulationResults[i][nCellParams * nCells + 0] = t;
             simulationResults[i][nCellParams * nCells + 1] = bod.totalEnergy();
            i++;
            t += dt;
        };

        testSimulation(dt, nSteps);
        String header = makeHeader(nCells);
        writeResults(simulationResults, header);

        // TODO: test consistency
        // double[] energies = getColumn(simulationResults, 19);
        // double e0 = energies[0];
        // double eav = average(energies);
        // assertTrue(String.format("Energy is not conserved! E0= %f, Eaverage=%f", e0, eav)
        // , eav <= e0);
    }

    void writeResults(double[][] simulationResults, String header) {
        try {
            Path path =  Paths.get("build/test-results/physics/SoftBody/CellGrid");
            if(Files.notExists(path)){
                    Files.createDirectories(path);
                }
            StringBuilder results = new StringBuilder(header);
            results.append(printGrid(simulationResults));
            Files.write(path.resolve((name.getMethodName() + ".dat")),
                results.toString().getBytes()
            );
        } catch(IOException e) {
            throw new RuntimeException("Could not write data");
        }
    }


    void testSimulation(double deltaT,int nSteps) {
        physics.Cell.dt = deltaT;
        t = 0;
        stepFunction = implicit;
        for(int i = 0; i < nSteps; i++) {
            stepFunction.accept(body);
            // System.out.format("step %d/%d\n", i, nSteps);
        }
    }
}
