package physics;
import physics.*;
import org.junit.Test;
import org.junit.Rule;
import org.junit.rules.*;
import org.ejml.data.*;
import java.io.*;
import java.util.*;
import java.nio.file.*;
import static org.junit.Assert.*;
import static org.ejml.dense.row.CommonOps_DDRM.*;
import static org.ejml.dense.row.NormOps_DDRM.*;
import static org.ejml.dense.fixed.NormOps_DDF2.*;
import static org.ejml.EjmlUnitTests.*;


public class ElasticModelTests {
    static final Path fDir = Paths.get("build/test-results/physics/ElasticModelTests");
    @Test
    public void compute() {
        double[] vertices1 = new double[] {0, 0, 0, 1, 1, 0};
        double[] vertices2 = new double[] {0, 0, 0, 1, 1, 1};
        double[] k = new double[] {1.0};
        double[] nu = new double[] {0.33};
        int[][] indices = new int[][] {{0,1,2}};
        double[] M = new double[] {1,1,1,1,1,1};
        ElasticModel em1 = new ElasticModel(vertices1, indices, k, nu, M);
        ElasticModel em2 = new ElasticModel(vertices2, indices, k, nu, M);
        DMatrix2x2 Bm1 = new DMatrix2x2(-1, -1,  0, 1);
        DMatrix2x2 Bm2 = new DMatrix2x2( 0, -1, -1, 1);
        assertEquals(em1.Bm[0], Bm1);
        assertEquals(em1.lambda[0], 0.729766, 1e-5);
        assertEquals(em1.mu[0], 0.37594, 1e-5);
        assertEquals(em2.Bm[0], Bm2);
        for(int i=0; i<vertices1.length; i++) {
            em1.x0.set(i, vertices1[i]);
            em1.x1.set(i, vertices1[i]);
        }

        em1.venantPiolaStress(new DMatrix2x2(1, 0, 0, 1), 1, 1, em1.temp2x2P);
        assertEquals(0, normF(em1.temp2x2P), 0);

        em1.computeElasticForces(em1.x0, em1.g);
        assertEquals(em1.g.toString(), 0, normP2(em1.g), 0);
    }

    @Test
    public void elasticForces() {
        double[] vertices1 = new double[] {0, 0, 0, 1, 1, 0};
        double[] vertices2 = new double[] {0, -1, 2, 2, 1, -1};
        double[][] vertices3 = new double[][] {{1}, {0}, {0}, {1}, {1}, {1}};
        double[] k = new double[] {1.0};
        double[] nu = new double[] {0.33};
        int[][] indices = new int[][] {{0,1,2}};
        double[] M = new double[] {1,1,1,1,1,1};
        ElasticModel em = new ElasticModel(vertices1, indices, k, nu, M);
        for(int i=0; i<vertices1.length; i++) {
            em.x0.set(i, vertices2[i]);
            em.x1.set(i, vertices1[i]);
        }

        em.computeElasticForces(em.x1, em.g);
        assertEquals(em.g.toString(), 0, normP2(em.g), 0);

        double psi = em.computeElasticForces(em.x0, em.g);
        DMatrixRMaj f = new DMatrixRMaj(new double[][] {
            {12.207}, {14.4626},
            {-9.26581}, {-13.3348},
            {-2.94118}, {-1.12782},
        });
        assertEquals(em.g, f, 1e-3);
        assertEquals(psi, 27.4215, 1e-4);


        DMatrixRMaj dx = new DMatrixRMaj(vertices3);
        em.computeElasticForceDifferential(em.x0, dx, em.g);
        f = new DMatrixRMaj(new double[][] {
            {-1.84653},{10.7364},
            {2.58735},{-7.04334},
            {-0.740823},{-3.69306}
        });
        assertEquals(em.g, f, 1e-5);
    }

    @Test
    public void venantPiolaStress() {
        double[] vertices1 = new double[] {0, 0, 0, 1, 1, 0};
        double[] k = new double[] {1.0};
        double[] nu = new double[] {0.33};
        int[][] indices = new int[][] {{0,1,2}};
        double[] M = new double[] {1,1,1,1,1,1};
        ElasticModel em1 = new ElasticModel(vertices1, indices, k, nu, M);
        em1.venantPiolaStress(new DMatrix2x2(1, 2, 0, 3), em1.lambda[0], em1.mu[0],
                                em1.temp2x2A);
        assertEquals(new DMatrix2x2(5.88235, 18.5316, 2.25564, 26.6696),
                     em1.temp2x2A, 1e-4);
    }

    @Test public void ejmlTests() {
        assertEquals(normP2(new DMatrixRMaj(new double[][] {
            {4},{3},{9},{2},{-5}
        })), 11.619, 1e-3);
    }
    @Test
    public void backwardEulerTest() {
        double nuV = 0.33;
        double kV =  1.0;
        double[] vertices = new double[] {0, 0, 0, 1, 1, 0};
        double[] k = new double[] {kV};
        double[] nu = new double[] {nuV};
        double[] M = new double[] {1,1,1,1,1,1};
        int[][] indices = new int[][] {{0,1,2}};
        ElasticModel em = new ElasticModel(vertices, indices, k, nu, M);
        Random rand = new Random(0);
        for(int i=0; i<vertices.length; i++) {
            em.x0.set(i, vertices[i] + 0.0*rand.nextDouble());
            em.x1.set(i, vertices[i] + 0.0*rand.nextDouble());
        }
        em.x1.add(1, 0, 0.7);
        em.x0.add(1, 0, 0.7);

        double dt=1, t=0;
        ImplicitODESolver.dt = dt;
        try {
            if(Files.notExists(fDir)) {
                Files.createDirectories(fDir);
            }
            BufferedWriter w = Files.newBufferedWriter(fDir.resolve("backwardEuler1.dat"));
            for(int j=0; j<2; j++) {
                w.write(String.format("% 8.4f", t));
                for(int i=0; i<vertices.length; i++) {
                    w.write(String.format(" % 8.4f ", em.x0.get(i)));
                }
                w.write("\n");
                em.implicitEulerStep();
                t+=dt;
                ImplicitODESolver.t = t;
            }
            w.close();
        }
        catch(IOException e) {
            System.out.println("Error in File IO " + e);
        }
    }

    @Test
    public void simpleForwardEuler() {
        double nuV = 0.33;
        double kV =  1.0;
        double[] vertices = new double[] {0, 0, 0, 1, 1, 0};
        double[] k = new double[] {kV};
        double[] nu = new double[] {nuV};
        int[][] indices = new int[][] {{0,1,2}};
        double[] M = new double[] {1,1,1,1,1,1};
        ElasticModel em = new ElasticModel(vertices, indices, k, nu, M);
        Random rand = new Random(0);
        for(int i=0; i<vertices.length; i++) {
            em.x0.set(i, vertices[i] + 0.0*rand.nextDouble());
            em.x1.set(i, vertices[i] + 0.0*rand.nextDouble());
        }
        em.x1.add(1, 0, 0.7);
        em.x0.add(1, 0, 0.7);

        double dt=0.01, t=0;

        try {
            if(Files.notExists(fDir)) {
                Files.createDirectories(fDir);
            }
            BufferedWriter w = Files.newBufferedWriter(fDir.resolve("simpleEuler1.dat"));
            for(int j=0; j<40; j++) {
                w.write(String.format("% 8.4f", t));
                for(int i=0; i<vertices.length; i++) {
                    w.write(String.format(" % 8.4f ", em.x0.get(i)));
                }
                w.write("\n");
                em.x2.set(em.x1);
                em.x1.set(em.x0);
                em.computeElasticForces(em.x1, em.g);
                // System.out.println("Force " + em.g);
                add(2, em.x1, -1, em.x2, em.x0);
                add(dt*dt, em.g, 1, em.x0, em.x0);
                t+=dt;
            }
            w.close();
        }
        catch(IOException e) {
            System.out.println("Error in File IO " + e);
        }
    }
}
