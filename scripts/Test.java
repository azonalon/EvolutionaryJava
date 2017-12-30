package scripts;
import java.util.*;
public class Test
{
    public static void main (String[] args) {
        Random r = new Random(42);
        for(int i=0; i<10; i++) {
            System.out.println(r.nextDouble() + "" +  r.nextDouble());
        }
    }

}
