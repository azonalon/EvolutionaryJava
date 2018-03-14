package physics;
import org.ejml.data.*;
import static org.ejml.EjmlUnitTests.*;
import static org.ejml.dense.row.CommonOps_DDRM.*;
import static org.ejml.dense.fixed.CommonOps_DDF2.*;
import static main.Main.DEVEL;
/**
 * Finite Element elasticity model.
 */
public class ElasticModel extends ImplicitODESolver {
    // matrices for triangles
    DMatrix2x2[] Bm; // reference triangle matrices
    DMatrix2x2 temp2x2A, temp2x2B, temp2x2P, temp2x2D, H, P, F, dF, dP;
    final static DMatrix2x2 nId = new DMatrix2x2(-1,0,0,-1);
    final static DMatrix2x2 id = new DMatrix2x2(1,0,0,1);
    static final int NEOHOOKEAN = 0;
    static final int VENANTKIRCHHOFF = 1;
    static int model = NEOHOOKEAN;
    int[][] Te; // nx3 indices for each triangle into point vector
    double[] W; // reference triangle volumes
    double[] mu; // first lame coefficient for each triangle
    double[] lambda; // second lame coefficient
    int n, m; // number of points and triangles

    void timeStepFinished() {

    }

    public static final boolean invertTranspose( DMatrix2x2 a , DMatrix2x2 inv ) {
        double scale = 1.0/elementMaxAbs(a);

        double a11 = a.a11*scale;
        double a12 = a.a12*scale;
        double a21 = a.a21*scale;
        double a22 = a.a22*scale;

        double m11 = a22;
        double m12 = -( a21);
        double m21 = -( a12);
        double m22 = a11;

        double det = (a11*m11 + a12*m12)/scale;

        inv.a11 = m11/det;
        inv.a21 = m21/det;
        inv.a12 = m12/det;
        inv.a22 = m22/det;

        return !Double.isNaN(det) && !Double.isInfinite(det);
    }

    /**
     * Generates a finite element elastic model solver.
     * Vertices are stored in a single vector (x1, y1, x2, y2..., xn, yn)
     * triangle indices are ((t11, t12, t13), (t21, t22, t23), ..., (tm1, tm2, tm3))
     * the tij specify offsets to the x coordinates of the vertex stored in vertices
     */
    public ElasticModel(double[] vertices, int[][] triangles,
                        double[] k, double[] nu, double[] M) {
        super(vertices.length);
        m = triangles.length;
        n = vertices.length; // dimension of parameter space
        this.Te = new int[m][3];
        for(int i=0; i<m; i++) {
            Te[i] = new int[] {
                2*triangles[i][0],
                2*triangles[i][1],
                2*triangles[i][2]
            };
        }
        this.MI = new DMatrixRMaj(n, 1);
        this.M = new DMatrixRMaj(n, 1);
        for(int i=0; i<n; i++) {
            this.M.set(i, 0, M[i]);
            this.MI.set(i, 0, 1/M[i]);
        }
        mu = new double[m];
        lambda = new double[m];
        W = new double[m];
        Bm = new DMatrix2x2[m];
        H  = new DMatrix2x2();
        F  = new DMatrix2x2();
        dF = new DMatrix2x2();
        P  = new DMatrix2x2();
        dP  = new DMatrix2x2();
        temp2x2P = new DMatrix2x2();
        temp2x2A = new DMatrix2x2();
        temp2x2B = new DMatrix2x2();
        temp2x2D = new DMatrix2x2();
        precompute(vertices, k, nu);
        for(int i=0; i<vertices.length; i++) {
            x0.set(i, vertices[i]);
            x1.set(i, vertices[i]);
            x2.set(i, vertices[i]);
        }
    }

    static final void set(DMatrix2x2 x, double a11, double a12, double a21, double a22) {
        x.a11 = a11;
        x.a12 = a12;
        x.a21 = a21;
        x.a22 = a22;
    }

    void precompute(double[] vertices, double[] K, double[] nu) {
        for(int l=0; l<this.m; l++){
            int i=Te[l][0], j=Te[l][1], k=Te[l][2];
            set(temp2x2A,
                  vertices[i + 0] - vertices[k + 0], vertices[j + 0] - vertices[k + 0],
                  vertices[i + 1] - vertices[k + 1], vertices[j + 1] - vertices[k + 1]
            );
            W[l] = Math.abs(det(temp2x2A)/2.0);
            Bm[l] = new DMatrix2x2();
            invert(temp2x2A, Bm[l]);
            assertCountable(Bm[l]);
            lambda[l] = K[l]*nu[l]/(1+nu[l])/(1-2*nu[l]);
            mu[l] = K[l]/2/(1+nu[l]);
        }
    }
    double computeElasticForces(DMatrixRMaj x, DMatrixRMaj dest) {
        dest.zero();
        double stressEnergy=0;
        for(int l=0; l<m; l++){
            int i=Te[l][0], j=Te[l][1], k=Te[l][2];
            set(temp2x2A,
                  x.get(i + 0, 0) - x.get(k + 0, 0), x.get(j + 0, 0) - x.get(k + 0, 0),
                  x.get(i + 1, 0) - x.get(k + 1, 0), x.get(j + 1, 0) - x.get(k + 1, 0)
            );
            mult(temp2x2A, Bm[l], temp2x2B);
            if(model == NEOHOOKEAN)
                stressEnergy += W[l]*neoHookeanStress(temp2x2B, lambda[l], mu[l], temp2x2A);
            else
                stressEnergy += W[l]*venantPiolaStress(temp2x2B, lambda[l], mu[l], temp2x2A);
            multTransB(temp2x2A, Bm[l], temp2x2B);
            addForceMatrixToVector(temp2x2B, dest, i, j, k, l);
        }
        return stressEnergy;
    }

    final static double normP2Squared(DMatrix2x2 m) {
        return
            m.a11*m.a11 +
            m.a12*m.a12 +
            m.a21*m.a21 +
            m.a22*m.a22;
    }

    final void addForceMatrixToVector(DMatrix2x2 H, DMatrixRMaj f,
                                             int i, int j, int k, int l) {
        f.add(i + 0, 0, -W[l]*H.a11);
        f.add(i + 1, 0, -W[l]*H.a21);
        f.add(j + 0, 0, -W[l]*H.a12);
        f.add(j + 1, 0, -W[l]*H.a22);
        f.add(k + 0, 0, +W[l]*H.a11 +W[l]*H.a12);
        f.add(k + 1, 0, +W[l]*H.a21 +W[l]*H.a22);
    }

    final double venantPiolaStress(DMatrix2x2 F,
                                 double lambda, double mu,
                                 DMatrix2x2 dest) {
        double psi=0;
        temp2x2P.set(nId);
        multAddTransA(F, F, temp2x2P);
        psi += mu*normP2Squared(temp2x2P)/4.0;
        double tr = trace(temp2x2P);
        double k = lambda * tr/2.0;
        psi += tr*tr*lambda/8.0;
        add(mu, temp2x2P, k, id, temp2x2P);
        mult(F, temp2x2P, dest);
        return psi;
    }

    static String str(DMatrix2x2 m) {
        return String.format("((%g, %g), (%g, %g))", m.a11, m.a12, m.a21, m.a22);
    }

    final double neoHookeanStress(DMatrix2x2 F,
                                 double lambda, double mu,
                                 DMatrix2x2 dest) {
        multTransA(F, F, temp2x2P);
        double I1 = trace(temp2x2P);
        // double J = det(F);
        double J = Math.abs(det(F));
        assert J > 0: String.format("Negative jacobian not supported for neo hookean! J=%g", J);
        double logJ = Math.log(J);
        invertTranspose(F, temp2x2P);
        add(mu, F, -mu+lambda*logJ, temp2x2P, dest);
        return mu*(I1/2.0 - 1 - logJ) + lambda/2.0*logJ*logJ;
    }

    final void neoHookeanStressDifferential(DMatrix2x2 F,DMatrix2x2 dF,
                                 double lambda, double mu,
                                 DMatrix2x2 dest) {
        // double J = det(F);
        double J = Math.abs(det(F));
        assert J > 0: String.format("Negative jacobian not supported for neo hookean! J=%g\n" + str(F), J);
        double logJ = Math.log(J);

        invertTranspose(F, temp2x2P);
        multTransB(temp2x2P, dF, temp2x2D);
        double trIFdF = trace(temp2x2D);
        mult(temp2x2D, temp2x2P, dest);
        add(mu, dF, mu - lambda*logJ, dest, dest);
        add(lambda*trIFdF, temp2x2P, 1, dest, dest);
    }

    final void venantPiolaStressDifferential(DMatrix2x2 F, DMatrix2x2 dF,
                                 double lambda, double mu,
                                 DMatrix2x2 destDP) {
        temp2x2P.set(nId);
        multAddTransA(F, F, temp2x2P);
        double k = lambda * trace(temp2x2P) / 2;
        add(mu, temp2x2P, k, id, temp2x2P);
        mult(dF, temp2x2P, destDP);

        multTransA(dF, F, temp2x2P);
        multAddTransA(F, dF, temp2x2P);
        double dk = lambda * trace(temp2x2P) / 2;
        add(mu, temp2x2P, dk, id, temp2x2P);
        multAdd(F, temp2x2P, destDP);
    }

    static final void add(double a, DMatrix2x2 A, double b, DMatrix2x2 B,
                     DMatrix2x2 dest) {
        dest.a11 = a * A.a11 + b * B.a11;
        dest.a12 = a * A.a12 + b * B.a12;
        dest.a21 = a * A.a21 + b * B.a21;
        dest.a22 = a * A.a22 + b * B.a22;
    }

    final double computeForce(DMatrixRMaj x, DMatrixRMaj dest) {
        return computeElasticForces(x, dest);
    }
    final void computeForceDifferential(DMatrixRMaj x, DMatrixRMaj dx,
                                               DMatrixRMaj dest) {
        computeElasticForceDifferential(x, dx, dest);
    }

    /**
     * Calculates the change of force df for small displacements around x.
     * f ~ f(x) + K dx = f(x) + df(x, dx)
     * df = K dx; where K ist the Stiffness Matrix.
     * The stiffness Matrix itself is not evaluated. But rather the response
     * of the Piola Kirchhoff Stress due to a small displacement.
     * This is useful for calculating the Newton direction for optimization
     * deltax = K^-1 f
     */
    void computeElasticForceDifferential(DMatrixRMaj x, DMatrixRMaj dx,
                                         DMatrixRMaj dest) {
        dest.zero();
        for(int l=0; l<m; l++){
            int i=Te[l][0], j=Te[l][1], k=Te[l][2];
            set(temp2x2A,
                  x.get(i + 0, 0) - x.get(k + 0, 0), x.get(j + 0, 0) - x.get(k + 0, 0),
                  x.get(i + 1, 0) - x.get(k + 1, 0), x.get(j + 1, 0) - x.get(k + 1, 0)
            );
            set(temp2x2B,
                  dx.get(i + 0, 0) - dx.get(k + 0, 0), dx.get(j + 0, 0) - dx.get(k + 0, 0),
                  dx.get(i + 1, 0) - dx.get(k + 1, 0), dx.get(j + 1, 0) - dx.get(k + 1, 0)
            );
            mult(temp2x2A, Bm[l], F);
            mult(temp2x2B, Bm[l], dF);
            if(model == NEOHOOKEAN)
                neoHookeanStressDifferential(F, dF, lambda[l], mu[l], dP);
            else
                venantPiolaStressDifferential(F, dF, lambda[l], mu[l], dP);
            multTransB(dP, Bm[l], temp2x2B);
            addForceMatrixToVector(temp2x2B, dest, i, j, k, l);
        }
    }
}
