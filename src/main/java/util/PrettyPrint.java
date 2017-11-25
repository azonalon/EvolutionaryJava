package util;

public class PrettyPrint {
    public static String printGrid(double[][] a)
    {
        String s = "";
        for(int i = 0; i < a.length; i++) {
            for(int j = 0; j < a[0].length; j++)
            {
                s += String.format("% 04.4f ", a[i][j]);
            }
            s += "\n";
        }
        return s;
    }
}
