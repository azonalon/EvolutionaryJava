package physics;
import Math.sqrt;
import Math.signum;
public class InvertibleNeoHookeanStressDensity {
	final static double square(double x) {
		return x*x;
	}
	final static double cube(double x) {
		return x*x*x;
	}
	final static double quart(double x) {
		return x*x*x*x;
	}
	final static public double[] computeStressFunctionComponents(
							double sx, double sy, double l, 
							double mu, double eps) 
	{
		double ux = sx - 1;
		double uy = sy - 1;
		double ux2 = ux*ux;
		double ux3 = ux2*ux;
		double ux4 = ux3*ux;
		double ux5 = ux4*ux;
		double uy2 = uy*uy;
		double uy3 = uy2*uy;
		double uy4 = uy3*uy;
		double uy5 = uy4*uy;
		double vx = sx+1;
		double vx2 = vx*vx;
		double vx3 = vx2*vx;
		double vy = sy+1;
		double vy2 = vy*vy;
		double vy3 = vy2*vy;
		double sy2 = sy*sy;
		double sx2 = sx*sx;
		double eps2 = eps*eps;
		double c0 = square(ux+uy);
		double c1 = c0+4*ux*uy*delta;
		assert c1 > 0;
		assert sx + sy >= 0; // IFE Convention
		double lEps = Math.log(eps);
		double k = mu -(-l+l*lEps);
		double bEps = Math.sqrt(c1);
		double delta = 1+eps;
		double delta2 = delta*delta;
		double delta3 = delta2*delta;
		double d = 2*bEps*ux2*uy2*eps2;
		double psi0=0, f1=0, f2=0, psi11=0, psi12=0, psi22=0;
		
		
		if(Math.abs(ux) < 1e-10) {
			psi0 = (bEps*ux2*uy*signum(uy)*(k*uy*(3 + 6*delta2 + 2*uy + 4*delta*vy) + 
			   signum(uy)*(bEps*k*(3 + 4*delta + (4 + 8*delta)*uy) + 
			     2*uy*(k*(3 + 8*delta - 2*eps2*lEps + 2*delta*uy*(3 + uy) + 
			         uy*(3 + 2*uy) + delta2*(6 + uy*vy)) - 
			       eps*l*(uy2 + delta*(2 + uy*vy) - eps*lEps*(2 + uy*vy) + 
			         eps*square(lEps))))))/(2*d);
		} else if (Math.abs(uy) < 1e-10) {
			psi0 = (bEps*ux*uy2*signum(ux)*(k*ux*(3 + 6*delta2 + 2*ux + 4*delta*vx) + 
			   signum(ux)*(bEps*k*(3 + 4*delta + (4 + 8*delta)*ux) + 
			     2*ux*(k*(3 + 8*delta - 2*eps2*lEps + 2*delta*ux*(3 + ux) + 
			         ux*(3 + 2*ux) + delta2*(6 + ux*vx)) - 
			       eps*l*(ux2 + delta*(2 + ux*vx) - eps*lEps*(2 + ux*vx) + 
			         eps*square(lEps))))))/(2*d);
		
		} else {
			psi0 = (bEps*(-2*eps2*lEps*ux2*uy2*(2*k + l*lEps - l*(2 + ux*vx + uy*vy)) + 
			   k*(bEps*c1*(ux + uy + 2*ux*uy) + 
			     ux4*(1 + 2*uy*(1 + (2 + delta*(2 + delta))*uy)) + 
			     2*ux3*uy*(2 + 3*delta + uy*(3 + 2*delta*(3 + delta + uy))) + 
			     uy4 + 2*ux*uy3*(3*delta + vy) + 2*ux2*uy2*
			      (3 + 8*delta + 2*delta*uy*(3 + uy) + uy*(3 + 2*uy) + 
			       delta2*(6 + uy*vy))) - 2*eps*l*ux2*uy2*
			    (delta*(2 + ux*vx + uy*vy) + square(ux - uy))))/(2*d);
		}
		
		if(Math.abs(ux) < 1e-10) {
			f1 = (ux2*signum(uy)*(bEps*k*(sx + sy)*uy*(-24*(1 + 2*delta) - 
			     uy*(44 + 3*uy*(9 + 2*uy)) + 2*delta2*
			      (-8 + (3 + uy)*(4 + uy)*uy2) + delta*uy*(-74 - 31*uy + 2*uy3) + 
			     2*delta3*sy*vy2) - (2*bEps*eps*l*(sx + sy)*uy3*vy3 + 
			     k*(uy*(-(uy*(80 + uy*(192 + uy*(162 + uy*(58 + 7*uy))))) - 
			         bEps*(32 + uy*(80 + uy*(96 + uy*(67 + uy*(25 + 4*uy))))) + 
			         delta*vy*(sy*uy*(-60 + uy*(-26 + uy*(13 + 8*uy))) + 
			           bEps*(-12 + uy*(-12 + uy*(5 + 2*uy*(4 + uy))*vy))) + 
			         4*delta2*uy*(-2 + uy*(4 + uy*(7 + 2*uy)))*vy2) + 
			       ux*(-2 + uy*(-2 + 5*delta + uy*(-7 + delta*(-3 - 4*delta + 8*
			                eps*uy))) + bEps*(1 - uy - 3*delta*uy + 2*delta*uy3))*
			        vy3))*signum(uy)))/(d*(sx + sy)*uy1*vy3);
		} else if (Math.abs(uy) < 1e-10) {
			f1 = (uy2*(bEps*k*(sx + sy)*ux*(-24*(1 + 2*delta) - 
			     ux*(44 + 3*ux*(9 + 2*ux)) + 2*delta2*
			      (-8 + (3 + ux)*(4 + ux)*ux2) + delta*ux*(-74 - 31*ux + 2*ux3) + 
			     2*delta3*sx*vx2) - 
			   (k*(-((80 + ux*(192 + ux*(162 + ux*(58 + 7*ux))))*ux2) + 
			       4*delta2*ux2*(-2*sy + ux*(4 + 3*uy + ux*(7 + 2*ux + 2*uy)))*
			        vx2 + delta*ux*vx*(sx*ux*(-60 + ux*(-26 + ux*(13 + 8*ux))) + 
			         (5 + ux*(-3 + 8*ux))*uy*vx2) - (2 + ux*(2 + 7*ux))*uy*vx3) + 
			     bEps*(2*eps*l*(sx + sy)*ux3*vx3 + 
			       k*(-(ux*(32 + ux*(80 + ux*(96 + ux*(67 + ux*(25 + 4*ux)))))) + 
			         delta*ux*(-12*sy + ux*(-12*sy + (sx + sy)*ux*(5 + 2*ux*
			                (4 + ux))))*vx - (-1 + ux)*uy*vx3)))/signum(ux))*
			  signum(ux))/(d*(sx + sy)*ux1*vx3);
		
		} else {
			f1 = (sx+sy<1e-10)?(1e10):(c1*k*(sy*ux4 + ux3*(1 + uy*(3 + delta - 2*delta*uy)) + 
			    ux2*uy*(3 + delta - 2*uy*(-2 + delta + delta*uy)) + 
			    ux*(3 + delta + uy*(3 + delta + uy))*uy2 + sy*uy3) + 
			  bEps*(-2*eps*l*(sx + sy)*ux3*uy3 + 
			    k*(sy*ux5 + ux4*(1 + uy*(4 + 3*delta + uy - 2*delta*uy2)) + 
			      sy*uy4 + ux2*uy2*(6 + 7*uy + uy2 + 3*delta*vy) + 
			      ux3*uy*(4 + 3*delta + uy*(7 + 3*delta + 4*uy - 
			          2*delta*uy*vy)) + ux*uy3*(3*delta*sy + vy2))))/
			 (d*(sx + sy)*ux1*uy1);
		}
		
		if(Math.abs(ux) < 1e-10) {
			f2 = (ux2*(4*k*(-(bEps*sy2*(4 + uy*vy*(2 - 3*delta + uy*vy))) - 
			     2*(-1 + delta*uy*vy)*(4 + uy*vy*(6 - delta + 3*uy*vy))) + 
			   4*delta*k*(sx + sy)*uy*(1 + 3*sy*uy)*vy3 + 
			   k*(sx + sy)*uy*(5 + delta + 3*uy*(4 + 3*uy))*vy3 + 
			   bEps*k*(sx + sy)*
			    (1 + uy*(3 + uy*(3 + 2*(2 + delta*(2 + delta))*uy)))*vy3 + 
			   2*bEps*eps2*l*(-1 + lEps)*(sx + sy)*uy3*vy3 - 
			   k*(sx + sy)*(1 + 3*sy*uy)*vy4 - bEps*k*(sx + sy)*
			    (24 + 48*delta + uy*(52 + 3*uy*(15 + uy*(6 + uy))) + 
			     delta*uy*(102 + uy*(89 + 6*uy*(6 + uy))) + 
			     2*delta2*vy*(4 + 3*uy*vy) - 2*delta3*vy2)*Abs(uy)))/
			 (d*(sx + sy)*uy1*vy3);
		} else if (Math.abs(uy) < 1e-10) {
			f2 = (bEps*eps2*ux2*uy2*((9*(-2 + eps)*k)/bEps + (12*delta*(-2 + eps)*k)/
			    bEps + 2*(-1 + eps)*k + 2*l*(-1 + lEps) + 
			   (3*(-2 + eps)*k*(2 + bEps + 4*delta - ux))/(bEps*ux1) + 
			   (3*(-2 + eps)*k)/ux2 + ((-2 + eps)*k*(-1 + 5*delta - 3*ux))/
			    (bEps*ux2) + ((-2 + eps)*k)/ux3 - ((-2 + eps)*k*vx)/(bEps*ux3) + 
			   (4*(-2 + eps)*k*(-(bEps*sx2*(4 + ux*vx*(2 - 3*delta + ux*vx))) - 
			      2*(-1 + delta*ux*vx)*(4 + ux*vx*(6 - delta + 3*ux*vx))))/
			    (bEps*(sx + sy)*ux3*vx3) - 
			   ((-2 + eps)*k*(24 + 48*delta + ux*(52 + 3*ux*(15 + ux*(6 + ux))) + 
			      delta*ux*(102 + ux*(89 + 6*ux*(6 + ux))) + 
			      2*delta2*vx*(4 + 3*ux*vx) - 2*delta3*vx2)*signum(ux))/
			    (ux2*vx3)))/d;
		
		} else {
			f2 = (sx+sy<1e-10)?(-1e10):(c1*k*(ux + uy)*(ux*uy*(2 + delta + 3*uy) + ux2*(1 + 3*sy*uy) + 
			    uy2) + bEps*(2*eps2*l*(-1 + lEps)*(sx + sy)*ux3*uy3 + 
			    k*(ux3*uy*(4 + 3*delta + uy*(9 + 6*delta + 
			          2*uy*(5 + 2*delta*(2 + delta) + (2 + delta*(2 + delta))*
			             uy))) + 
			      ux4*(1 + uy*(3 + uy*(3 + 2*(2 + delta*(2 + delta))*uy))) + 
			      ux*(4 + 3*delta + 3*uy)*uy3 + uy4 + 3*sy*ux2*uy2*
			       (2*delta + vy))))/(d*(sx + sy)*ux1*uy1);
		}
		
		if(Math.abs(ux) < 1e-10) {
			psi11 = (ux2*(k*uy*(uy + 2*uy*(delta + uy + 2*delta*uy) + 
			     bEps*(1 + 2*uy*(1 + (2 + delta*(2 + delta))*uy))) + 
			   2*bEps*eps2*l*(-1 + lEps)*uy3 + 2*bEps*delta*eps*k*
			    (-1 + delta + 3*delta2 - 2*(1 + 2*delta)*uy)*Abs(uy)))/(d*uy1);
		} else if (Math.abs(uy) < 1e-10) {
			psi11 = (2*uy2*((bEps*eps2*l*(-1 + lEps)*ux + 
			     k*(-(delta*eps) + (1 + 2*delta + bEps*(2 + delta*(2 + delta)))*
			        ux))*ux2 + bEps*k*(delta + delta2 - ux - 2*delta*ux)*
			    Abs(ux)))/(d*ux1);
		
		} else {
			psi11 = (-2*bEps*eps2*l*ux4*uy2 + 2*bEps*eps2*l*lEps*ux4*uy2 + 
			  k*(ux5*(1 + 2*uy) + ux4*(bEps + uy*(1 + 2*bEps + 2*delta + 
			        2*(1 + 2*delta + bEps*(2 + delta*(2 + delta)))*uy)) - 
			    2*delta*eps*ux3*(1 + 2*uy)*uy2 + 3*(bEps + uy)*uy4 + 
			    2*ux2*uy3*(delta*(4 + 3*delta) + 2*delta*uy + vy) + 
			    ux*uy3*(uy*(7 + 12*delta + 2*uy) + 2*bEps*(3*delta + vy))))/
			 (d*ux2);
		}
		
		if(Math.abs(ux) < 1e-10) {
			psi12 = (ux2*(k*uy*(uy*(-2 + delta*(-3 - 2*delta + 4*eps*uy)) + 
			     bEps*(-2 + delta*(-3 + 2*uy2))) + 2*bEps*eps*l*uy3 + 
			   2*bEps*delta*eps*k*(delta - 2*uy)*Abs(uy)))/(d*uy1);
		} else if (Math.abs(uy) < 1e-10) {
			psi12 = (uy2*(k*ux*(ux*(-2 + delta*(-3 - 2*delta + 4*eps*ux)) + 
			     bEps*(-2 + delta*(-3 + 2*ux2))) + 2*bEps*eps*l*ux3 + 
			   2*bEps*delta*eps*k*(delta - 2*ux)*Abs(ux)))/(d*ux1);
		
		} else {
			psi12 = (2*bEps*eps*l*ux3*uy3 - 
			  k*(2*sy*ux5 + ux4*(2*bEps*sy + uy*(4 + 7*delta + 
			        (2 + 4*delta)*uy)) + ux3*uy*(bEps*(2 + 3*delta) + 
			      (2 + delta*(3 + 2*delta))*uy - 2*delta*(2 + bEps + 2*delta)*
			       uy2) + ux2*(2 + delta*(3 + 2*delta) + (2 + 4*delta)*uy)*uy3 + 
			    ux*(bEps*(2 + 3*delta + 2*uy) + uy*(4 + 7*delta + 2*uy))*uy3 + 
			    2*(bEps + uy)*uy4))/(d*ux1*uy1);
		}
		
		if(Math.abs(ux) < 1e-10) {
			psi22 = (2*ux2*((bEps*eps2*l*(-1 + lEps)*uy + 
			     k*(-(delta*eps) + (1 + 2*delta + bEps*(2 + delta*(2 + delta)))*
			        uy))*uy2 + bEps*k*(delta + delta2 - uy - 2*delta*uy)*
			    Abs(uy)))/(d*uy1);
		} else if (Math.abs(uy) < 1e-10) {
			psi22 = (uy2*(k*ux*(ux + 2*ux*(delta + ux + 2*delta*ux) + 
			     bEps*(1 + 2*ux*(1 + (2 + delta*(2 + delta))*ux))) + 
			   2*bEps*eps2*l*(-1 + lEps)*ux3 + 2*bEps*delta*eps*k*
			    (-1 + delta + 3*delta2 - 2*(1 + 2*delta)*ux)*Abs(ux)))/(d*ux1);
		
		} else {
			psi22 = (-2*bEps*eps2*l*ux2*uy4 + 2*bEps*eps2*l*lEps*ux2*uy4 + 
			  k*(ux5*(3 + 2*uy) + ux4*(3*bEps + uy*(7 + 2*bEps + 2*uy + 
			        4*delta*(3 + uy))) + 2*ux3*uy*(bEps*(2 + 3*delta) + 
			      uy*(2 + delta*(4 + 3*delta - 2*eps*uy))) + 
			    2*ux2*(-(delta*eps) + (1 + 2*delta + 
			        bEps*(2 + delta*(2 + delta)))*uy)*uy3 + (bEps + uy)*uy4 + 
			    ux*(1 + 2*bEps + 2*delta + 2*uy)*uy4))/(d*uy2);
		}
		
		
		return new double[] {psi0, f1, f2, psi11, psi12, psi22};
	}
}
