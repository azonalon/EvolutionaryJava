package physics;
import org.ejml.data.*;
import static org.ejml.EjmlUnitTests.*;
import static org.ejml.dense.row.CommonOps_DDRM.*;
import static org.ejml.dense.fixed.CommonOps_DDF2.*;
import static main.Main.DEVEL;
import static util.Math.invertTranspose;
import static util.Math.add;
import static util.Math.str;
import static util.Math.singularValueDecomposition;
import static physics.InvertibleNeoHookeanStressDensity.*;
import java.util.Arrays;
/**
 * Invertible neo-hookean stress tensor computation
 */
public class InvertibleNeoHookeanModel {
    DMatrix2x2 temp2x2A, temp2x2B, temp2x2P, temp2x2D, b;
    double eps = 0.3;
    InvertibleNeoHookeanModel(double eps) {
        temp2x2P = new DMatrix2x2();
        temp2x2A = new DMatrix2x2();
        temp2x2B = new DMatrix2x2();
        temp2x2D = new DMatrix2x2();
        b = new DMatrix2x2();
        assert(eps > 0.01 && eps < 1);
        this.eps = eps;
    }

    final double computeStressTensor(DMatrix2x2 F,
                                 double lambda, double mu,
                                 DMatrix2x2 dest) {
        double J = det(F);
        if(J > eps) {
            // normal neohookean
            multTransA(F, F, temp2x2P);
            double I1 = trace(temp2x2P);
            double logJ = Math.log(J);
            invertTranspose(F, temp2x2P);
            add(mu, F, -mu+lambda*logJ, temp2x2P, dest);
            return mu*(I1/2.0 - 1 - logJ) + lambda/2.0*logJ*logJ;
        } else {
            // extrapolate
            final double[] svd = singularValueDecomposition(F);
            final double u1 = svd[0], u2 = svd[1],
                        sx = svd[2], sy = svd[3],
                        v1 = svd[4], v2 = svd[5];
            final double[] r = computeStressGradientComponents(sx, sy, lambda, mu, eps);
            assert Math.abs(J-sx*sy) < 1e-9;
            final double psi0 = r[0], psi1=r[1], psi2=r[2];
            dest.a11 = psi1*u1*v1 - psi2*u2*v2;
            dest.a12 = -(psi2*u2*v1) - psi1*u1*v2;
            dest.a21 = psi1*u2*v1 + psi2*u1*v2;
            dest.a22 = psi2*u1*v1 - psi1*u2*v2;
            return psi0;
        }
    }

    final void computeStressDifferential(DMatrix2x2 F,DMatrix2x2 dF,
                                 double lambda, double mu,
                                 DMatrix2x2 dest) {
        // double J = det(F);
        double J = det(F);
        if(J>eps) {
            double logJ = Math.log(J);
            // normal neohookean
            invertTranspose(F, temp2x2P);
            multTransB(temp2x2P, dF, temp2x2D);
            double trIFdF = trace(temp2x2D);
            mult(temp2x2D, temp2x2P, dest);
            add(mu, dF, mu - lambda*logJ, dest, dest);
            add(lambda*trIFdF, temp2x2P, 1, dest, dest);
        } else {
            // extrapolate
            double[] svd = singularValueDecomposition(F);
            final double u1 = svd[0], u2 = svd[1],
                    sx = svd[2], sy = svd[3],
                    v1 = svd[4], v2 = svd[5];
            final double[] r = computeStressDifferentialComponents(sx, sy, lambda, mu, eps);
            assert Math.abs(J-sx*sy) < 1e-9;
            final double psi0 = r[0], f1=r[1], f2=r[2], psi11=r[3], psi12=r[4], psi22=r[5];
            temp2x2A.a11=u1;temp2x2A.a12=-u2;
            temp2x2A.a21=u2;temp2x2A.a22=u1;
            temp2x2B.a11=v1;temp2x2B.a12=-v2;
            temp2x2B.a21=v2;temp2x2B.a22=v1;
            multTransB(dF, temp2x2B, b);
            multTransA(temp2x2A, dF, temp2x2D);
            multTransB(temp2x2D, temp2x2B, b);
            temp2x2D.a11 = b.a11*psi11 + b.a22*psi12;
            temp2x2D.a12 = b.a21*f1 + b.a12*f2;
            temp2x2D.a21 = b.a12*f1 + b.a21*f2;
            temp2x2D.a22 = b.a11*psi12 + b.a22*psi22;
            mult(temp2x2A, temp2x2D, temp2x2P);
            mult(temp2x2P, temp2x2B, dest);
        }
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

}
