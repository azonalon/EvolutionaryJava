package scripts;
public class ArcTan {
    public static void main(String[] args) {
        double s, x, y;
        for(int i=0; i<100; i++) {
            s = (double)i/99.0 * 2 * Math.PI - Math.PI;
            x = Math.cos(s);
            y = Math.sin(s);
            System.out.format("%f %f %f %f\n", s,
            Math.atan(y/x), Math.atan2(y, x), 2*Math.atan(y/(1+x))
            );
        }
        System.out.format("atan2(0, 0)=%g, atan2(0, -1)=%g",
                Math.atan2(0.0, 0.0), Math.atan2(0.0, -1.0));
    }
}
