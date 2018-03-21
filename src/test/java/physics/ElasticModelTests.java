package physics;
import physics.*;
import org.junit.Test;
import org.junit.Rule;
import org.junit.rules.*;
import org.ejml.data.*;
import java.io.*;
import java.util.*;
import java.util.function.*;
import java.nio.file.*;
import static org.junit.Assert.*;
import static org.ejml.dense.row.CommonOps_DDRM.*;
import static org.ejml.dense.row.NormOps_DDRM.*;
import static org.ejml.dense.fixed.CommonOps_DDF2.*;
import static org.ejml.dense.fixed.NormOps_DDF2.*;
import static org.ejml.EjmlUnitTests.*;
import static org.ejml.ops.MatrixIO.*;


public class ElasticModelTests {
    static final Path fDir = Paths.get("build/test-results/physics/ElasticModelTests");
    @Test
    public void conjugateGradient() {
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

        final DMatrixRMaj A = new DMatrixRMaj(new double[][]
        {{10, 0, 4, -2, 1, 0}, {0, 6, 0, 0, 5, 0}, {4, 0, 6, -1, -1, 2}, {-2,
          0, -1, 10, 5, 2}, {1, 5, -1, 5, 10, -1}, {0, 0, 2, 2, -1, 6}}
        );
        DMatrixRMaj b = new DMatrixRMaj(new double[][]
            {{-4}, {-9}, {-9}, {10}, {3}, {-9}}
        );
        DMatrixRMaj x0 = new DMatrixRMaj(new double[][]
            {{0}, {0}, {0}, {0}, {0}, {0}}
        );
        DMatrixRMaj xCorrect = new DMatrixRMaj(new double[][]
            {{0.0223104}, {-2.05611}, {-0.784498}, {0.876219}, {0.667328},
            {-1.41935}}
        );
        DMatrixRMaj x = new DMatrixRMaj(6, 1);
        BiConsumer<DMatrixRMaj, DMatrixRMaj> computeLhs = (m, d) -> {
            mult(A, m, d);
        };
        em.conjugateGradientSolve(b, x0, computeLhs, x);
        assertEquals(x, xCorrect, 1e-3);

        final DMatrixRMaj B = new DMatrixRMaj(new double[][]
          {{10,-6,2,6,-6,-1},{-6,2,3,3,-1,-6},{2,3,8,-2,-3,-6},{6,3,-2,-8,2,-5},{-6,-1,-3,2,-6,-1},{-1,-6,-6,-5,-1,-4}}
        );
        b = new DMatrixRMaj(new double[][]
            {{-10},{-3},{0},{4},{-1},{-5}}
        );
        x0 = new DMatrixRMaj(new double[][]
            {{0}, {0}, {0}, {0}, {0}, {0}}
        );
        xCorrect = new DMatrixRMaj(new double[][]
            {{0.0115413},{-0.631618},{0.193112},{-0.663991},{-1.42397},{-2.29857}}
        );
        x = new DMatrixRMaj(6, 1);
        computeLhs = (m, d) -> {
            mult(B, m, d);
        };
        em.conjugateGradientSolve(b, x0, computeLhs, x);
        assertNotEquals(x, xCorrect);
    }

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
        // double[] vertices2 = new double[] {0, -1, 2, 2, 1, -1};
        double[] vertices2 = new double[] {0, 0.7, 0, 1, 1, 0};
        double[][] vertices3 = new double[][] {{-.05}, {.1}, {-.2}, {.3}, {.4}, {.5}};
        double[] k = new double[] {1.0};
        double[] nu = new double[] {0.33};
        int[][] indices = new int[][] {{0,1,2}};
        double[] M = new double[] {1,1,1,1,1,1};
        ElasticModel em = new ElasticModel(vertices1, indices, k, nu, M, ElasticModel.VENANTKIRCHHOFF);
        for(int i=0; i<vertices1.length; i++) {
            em.x0.set(i, vertices2[i]);
            em.x1.set(i, vertices1[i]);
        }
        DMatrixRMaj dx = new DMatrixRMaj(vertices3);
        double psi;


        em.computeElasticForces(em.x1, em.g);
        assertEquals(em.g.toString(), 0, normP2(em.g), 0);

        psi = em.computeElasticForces(em.x0, em.g);
        DMatrixRMaj f = new DMatrixRMaj(new double[][]
        {{-0.0239938}, {-0.0693498}, {0.0394737}, {0.0466718}, {-0.0154799}, {0.022678}}
        );
        assertEquals(em.g, f, 1e-4);
        assertEquals(psi, 0.0623878, 1e-4);

        em.computeElasticForceDifferential(em.x0, dx, em.g);
        f = new DMatrixRMaj(new double[][]
        {{0.148153}, {-0.125785}, {0.0125663}, {0.0110128}, {-0.16072}, {0.114772}}
        );
        assertEquals(em.g, f, 1e-4);

        em = new ElasticModel(vertices1, indices, k, nu, M, ElasticModel.NEOHOOKEAN);
        for(int i=0; i<vertices1.length; i++) {
            em.x0.set(i, vertices2[i]);
            em.x1.set(i, vertices1[i]);
        }

        psi = em.computeElasticForces(em.x0, em.g);
        f = new DMatrixRMaj(new double[][]
            {{-1.90296}, {-2.16612}, {1.46365}, {2.03454}, {0.439309}, {0.131579}}
        );
        assertEquals(em.g, f, 1e-4);
        assertEquals(psi, 0.451295, 1e-4);

        em.computeElasticForceDifferential(em.x0, dx, em.g);
        f = new DMatrixRMaj(new double[][]
            {{3.0101}, {1.39375}, {-2.58304}, {-1.6322}, {-0.427058}, {0.238451}}
        );
        assertEquals(em.g, f, 1e-4);

        DMatrixRMaj m = new DMatrixRMaj(6,6);
        for(int i=0; i<6; i++) {
            dx = new DMatrixRMaj(6, 1);
            dx.set(i, 0, 1.0);
            em.computeElasticForceDifferential(em.x0, dx, em.g);
            for(int j=0; j<6;j++) {
                // assertEquals(t[i][j], em.g.get(j, 0), 1e-4);
                // System.out.format("%g %g %d %d\n", t[i][j], em.g.get(j, 0), i, j);
                m.set(i, j, em.g.get(j, 0));
            }
        }
    }

    @Test
    public void elasticForcesInvertibleNeoHookean() {
        double[] vertices1 = new double[] {0, 0, 0, 1, 1, 1};
        // double[] vertices2 = new double[] {0, -1, 2, 2, 1, -1};
        double[] vertices2 = new double[] {1.2, 1.1, 0.1, 1.05, 1.06, 1.07};
        double[][] vertices3 = new double[][] {{-.05}, {.1}, {-.2}, {.3}, {.4}, {.5}};
        double[] k = new double[] {1.0};
        double[] nu = new double[] {0.33};
        int[][] indices = new int[][] {{0,1,2}};
        double[] M = new double[] {1,1,1,1,1,1};
        DMatrixRMaj dx = new DMatrixRMaj(vertices3);
        double psi;


        ElasticModel em = new ElasticModel(vertices1, indices, k, nu, M,
                                            ElasticModel.INVERTIBLE_NEOHOOKEAN, 0.3);
        for(int i=0; i<vertices1.length; i++) {
            em.x0.set(i, vertices2[i]);
            em.x1.set(i, vertices1[i]);
        }

        psi = em.computeElasticForces(em.x0, em.g);
        DMatrixRMaj f = new DMatrixRMaj(new double[][]
            {{-0.046988}, {-4.7444}, {0.420897}, {-0.681655}, {-0.373909},  {5.42606}}
        );
        assertEquals(1.68113, psi, 1e-4);
        assertEquals(em.g, f, 1e-4);
        double [][] t = new double [][] {{-0.356551, 0.646852, 0.506009, -4.82607, -0.149458,
                          4.17922}, {0.646852, -6.77856, 4.19242, -0.872466, -4.83927,
                          7.65102}, {0.506009, 4.19242, -0.911141, -0.133457,
                          0.405132, -4.05896}, {-4.82607, -0.872466, -0.133457, -0.417634,
                          4.95953, 1.2901}, {-0.149458, -4.83927, 0.405132,
                          4.95953, -0.255674, -0.120252}, {4.17922, 7.65102, -4.05896,
                          1.2901, -0.120252, -8.94112}};

        DMatrix2x2 F = new DMatrix2x2(0.96,-1.1,0.02,-0.05);
        DMatrix2x2 dF = new DMatrix2x2(0.05,-0.07,-0.07,0.02);
        DMatrix2x2 dP = new DMatrix2x2();
        DMatrix2x2 t1 = new DMatrix2x2(-0.163915, -0.675827, -1.51878, -1.19337);
        InvertibleNeoHookeanModel neo = new InvertibleNeoHookeanModel(0.3);
        neo.computeStressDifferential(F, dF, 0.729766, 0.37594, dP);
        assertEquals(t1, dP, 0.0001);


        DMatrixRMaj m = new DMatrixRMaj(6,6);
        for(int i=0; i<6; i++) {
            dx = new DMatrixRMaj(6, 1);
            dx.set(i, 0, 1.0);
            em.computeElasticForceDifferential(em.x0, dx, em.g);
            for(int j=0; j<6;j++) {
                assertEquals(t[i][j], em.g.get(j, 0), 1e-4);
                // System.out.format("%g %g %d %d\n", t[i][j], em.g.get(j, 0), i, j);
                m.set(i, j, em.g.get(j, 0));
            }
        }
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
    static final double[] flatten(double[][] a) {
        int n= a.length;
        int m= a[0].length;
        double[] f = new double[n*m];
        for(int i=0; i<n; i++) {
            for(int j=0; j<m; j++) {
                f[m*i +j] = a[i][j];
            }
        }
        return f;
    }
    @Test
    public void determinant() {
        DMatrix2x2 m = new DMatrix2x2(1, 0, 0, 1);
        assertEquals(det(m), 1.0, 1e-10);
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

        em.x0.add(1, 0, 0.7);
        em.x1.add(1, 0, 0.7);

        double dt=0.1, t=0;
        em.kDamp = 0.5;
        ImplicitODESolver.dt = dt;
        try {
            if(Files.notExists(fDir)) {
                Files.createDirectories(fDir);
            }
            BufferedWriter w = Files.newBufferedWriter(fDir.resolve("backwardEulerTest.dat"));
            for(int j=0; t<10; j++) {
                w.write(String.format("% 8.4f", t));
                for(int i=0; i<vertices.length; i++) {
                    w.write(String.format(" % 8.4f ", em.x0.get(i)));
                }
                em.fExt.set(0, 0, 1*Math.cos(Math.PI*t));
                w.write(String.format("% 8.4f", em.dampPotential));
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
    public void backwardEulerTestTwoTriangle() {
        double nuV = 0.33;
        double kV =  1.0;
        double[] vertices = new double[] {0, 0, 0, 1, 1, 0, 1, 1};
        double[] k = new double[] {kV, kV};
        double[] nu = new double[] {nuV, nuV};
        double[] M = new double[] {1,1,1,1,1,1,1, 1};
        int[][] indices = new int[][] {{0,1,2}, {1, 2, 3}};
        ElasticModel em = new ElasticModel(vertices, indices, k, nu, M);

        // em.x0.add(1, 0, 0.7);
        // em.x1.add(1, 0, 0.7);

        double dt=0.05, t=0;
        double u = 0.8;
        em.kDamp = 0.0;
        ImplicitODESolver.dt = dt;
        try {
            if(Files.notExists(fDir)) {
                Files.createDirectories(fDir);
            }
            BufferedWriter w = Files.newBufferedWriter(fDir.resolve("backwardEulerTestTwoTriangle.dat"));
            for(int j=0; t<10; j++) {
                w.write(String.format("% 8.4f", t));
                for(int i=0; i<vertices.length; i++) {
                    w.write(String.format(" % 8.4f ", em.x0.get(i)));
                }
                em.fExt.set(0, 0, u);
                em.fExt.set(1, 0, u);
                em.fExt.set(6, 0, -u);
                em.fExt.set(7, 0, -u);
                w.write(String.format("% 8.4f", em.dampPotential));
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

    public static void readObj(Path file, Vector<double[]> vertices,
                                   Vector<int[]> indices) throws IOException {

        Scanner s = new Scanner(Files.newBufferedReader(
            file
        ));
        while(s.hasNextLine()) {
            Scanner ls = new Scanner(s.nextLine());
            if(!ls.hasNext()) {
                continue;
            }
            String token = ls.next();
            if("v".equals(token)) {
                vertices.addElement(new double[]{ls.nextDouble(), ls.nextDouble()});
            } else if("f".equals(token)) {
                indices.addElement(new int[3]);
                for(int i=0; i<3; i++) {
                    Scanner iScan = new Scanner(ls.next());
                    iScan.useDelimiter("/");
                    indices.lastElement()[i] = iScan.nextInt() - 1;
                    iScan.close();
                }
            }
            ls.close();
        }
        s.close();
        // for(int i=0; i<vertices.size(); i++) {
        //     for(int j=0; j< 2; j++) {
        //         System.out.print(vertices.get(i)[j] + " ");
        //     }
        //     System.out.println();
        // }
        // for(int i=0; i<indices.size(); i++) {
        //     for(int j=0; j< 3; j++) {
        //         System.out.print(indices.get(i)[j] + " ");
        //     }
        //     System.out.println();
        // }
    }
    void saveIndices(Path path, int[][] indices) throws IOException {
        BufferedWriter w = Files.newBufferedWriter(path);
        for(int[] face: indices) {
            w.write(String.format("% 4d % 4d % 4d\n", face[0], face[1], face[2]));
        }
        w.close();
    }

    void saveIndices(Path path, Iterable<int[]> indices) throws IOException {
        BufferedWriter w = Files.newBufferedWriter(path);
        for(int[] face: indices) {
            w.write(String.format("% 4d % 4d % 4d\n", face[0], face[1], face[2]));
        }
        w.close();
    }

    @Test
    public void ball() {
        double nuV = 0.4;
        double kV =  50.0;
        ballSimulation(ElasticModel.INVERTIBLE_NEOHOOKEAN, kV, nuV);
    }

    @Test
    public void triangleInversion() {
        double nuV = 0.3;
        double kV =  0.5;
        try {
            double[] vertices = new double[] {0, 0, 1, 0, 0, 1};
            int[][] indices = new int[][] {{0,1,2}};

            double[] k = new double[] {kV} ;
            double[] nu = new double[] {nuV} ;
            double[] M = new double[] {1,1,1,1,1,1};
            ElasticModel em = new ElasticModel(vertices, indices,
                                               k, nu, M,
                                               ElasticModel.INVERTIBLE_NEOHOOKEAN,
                                               0.7);
            em.x1.set(0, 0, 1);
            em.x1.set(1, 0, 1);
            em.x0.set(em.x1);
            em.x2.set(em.x1);
            System.out.println(em.x0);
            System.out.println(em.x1);
            System.out.println(em.x2);

            double dt=0.01, t=0;
            em.kDamp = 0.3;
            ImplicitODESolver.dt = dt;
            if(Files.notExists(fDir)) {
                Files.createDirectories(fDir);
            }
            saveIndices(fDir.resolve("triangle.indices"), indices);
            BufferedWriter w = Files.newBufferedWriter(
                fDir.resolve("triangle.dat")
            );
            for(int j=0; t<50; j++) {
                w.write(String.format("% 8.4f", t));
                for(int i=0; i<em.x0.numRows; i++) {
                    w.write(String.format(" % 8.4f ", em.x0.get(i)));
                }
                if(t < 5) {
                    em.fExt.set(0, 0, 0);
                    em.fExt.set(1, 0, 0);
                    em.fExt.set(2, 0, -0);
                    em.fExt.set(3, 0, 0);
                    em.fExt.set(4, 0, 0);
                    em.fExt.set(5, 0, -0);
                }
                else {
                    em.fExt.set(0, 0, 0);
                    em.fExt.set(1, 0, 0);
                    em.fExt.set(2, 0, 0);
                    em.fExt.set(3, 0, 0);
                    em.fExt.set(4, 0, 0);
                    em.fExt.set(5, 0, 0);
                }
                w.write("\n");
                em.implicitEulerStep();
                t+=dt + 0*j;
                ImplicitODESolver.t = t;
            }
            w.close();
        }
        catch(IOException e) {
            System.out.println("Error in File IO " + e);
        }

    }

    public void ballSimulation(int model, double kV, double nuV) {
        try {
            Vector<double[]> vertices = new Vector<double[]>();
            Vector<int[]> indices = new Vector<int[]>();
            readObj(Paths.get("src/test/resources/physics/ball.obj"),
            vertices, indices);

            double[] k = new double[indices.size()] ;
            double[] nu = new double[indices.size()] ;
            double[] M = new double[vertices.size()*2] ;
            Arrays.fill(k, kV);
            Arrays.fill(nu, nuV);
            Arrays.fill(M, 1);
            ElasticModel em = new ElasticModel(flatten(vertices.toArray(new double[0][])),
                                               indices.toArray(new int[0][]),
                                               k, nu, M, model, 0.7);
            em.x0.set(em.x1);
            em.x2.set(em.x1);

            double dt=0.05, t=0;
            em.kDamp = 0.1;
            ImplicitODESolver.dt = dt;
            if(Files.notExists(fDir)) {
                Files.createDirectories(fDir);
            }
            saveIndices(fDir.resolve("ball.indices"), indices);
            BufferedWriter w = Files.newBufferedWriter(
                fDir.resolve("ball.dat")
            );
            for(int j=0; t<50; j++) {
                w.write(String.format("% 8.4f", t));
                for(int i=0; i<em.x0.numRows; i++) {
                    w.write(String.format(" % 8.4f ", em.x0.get(i)));
                }
                for(int i=0; i<em.x0.numRows; i++) {
                    if(i%2 == 1) {
                        em.fExt.set(i, 0, -0.2);
                        if(em.x0.get(i) < -5) {
                            em.fExt.set(i, 0, 10);
                        }
                    }
                }
                w.write("\n");
                // em.fExt.set(5, 0, -1);
                // em.fExt.set(3, 0, 1);
                em.implicitEulerStep();
                t+=dt + 0*j;
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
        double kV =  0.01;
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

        double dt=0.05, t=0;
        ImplicitODESolver.dt = dt;

        try {
            if(Files.notExists(fDir)) {
                Files.createDirectories(fDir);
            }
            BufferedWriter w = Files.newBufferedWriter(fDir.resolve("simpleForwardEuler.dat"));
            for(int j=0; t<100; j++) {
                w.write(String.format("% 8.4f", t));
                for(int i=0; i<vertices.length; i++) {
                    w.write(String.format(" % 8.4f ", em.x0.get(i)));
                }
                w.write("\n");
                em.x2.set(em.x1);
                em.x1.set(em.x0);
                em.computeForce(em.x1, em.g);
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
