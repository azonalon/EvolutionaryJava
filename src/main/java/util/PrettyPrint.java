package util;

public class PrettyPrint {
    public static String printGrid(double[][] a)
    {
        StringBuilder sb = new StringBuilder(a.length * a[0].length * 10);
        for(int i = 0; i < a.length; i++) {
            for(int j = 0; j < a[0].length; j++)
            {
                sb.append(String.format("% 04.4f ", a[i][j]));
            }
            sb.append("\n");
        }
        return sb.toString();
    }
}
