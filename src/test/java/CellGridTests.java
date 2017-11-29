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
import static util.Math.*;

public class CellGridTests
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

    public void runTwoCellTestCase(CellCulture[][] grid, double dt, int nSteps) {
        int nCells = grid.length*grid[0].length;
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

        testSimulation(grid, nSteps);
        String header = makeHeader(nCells);
        writeResults(simulationResults, header);

        // TODO: test consistency
        double[] energies = getColumn(simulationResults, 19);
        double e0 = energies[0];
        double eav = average(energies);
        assertTrue(String.format("Energy is not conserved! E0= %f, Eaverage=%f", e0, eav)
        , eav <= e0);
    }

    void writeResults(double[][] simulationResults, String header) {
        try {
            Path path =  Paths.get("build/test-results/cellgrid");
            if(Files.notExists(path)){
                    Files.createDirectory(path);
                }
            Files.write(Paths.get(name.getMethodName() + ".dat").resolve(path),
                (header + printGrid(simulationResults)).getBytes()
            );
        } catch(IOException e) {
            throw new RuntimeException("Could not write data");
        }
    }


    static void testSimulation(CellCulture[][] grid, int nSteps) {
        SoftBody bod = new SoftBody(grid);
        t = 0;
        for(int i = 0; i < nSteps; i++) {
            bod.updateForces();
            bod.propagateCells();
        }
    }
}
