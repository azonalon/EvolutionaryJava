public class ArcTan {
    public static void main(String[] args) {
        double s, x, y;
        for(int i=0; i<100; i++) {
            s = (double)i/99.0 * 2 * Math.PI - Math.PI/2.0;
            x = Math.cos(s);
            y = Math.sin(s);
            System.out.format("%f %f %f %f\n", s,
            Math.atan(y/x), Math.atan2(y, x), 2*Math.atan(y/(1+x))
            );
        }
    }
}
