package renderable;
import shader.ShaderPipeline;
import java.nio.*;
import renderable.*;
import static org.lwjgl.system.MemoryUtil.*;
import org.lwjgl.BufferUtils;

public class Triangle extends Renderable {
    public Triangle (ShaderPipeline pipeline) {
        FloatBuffer vb = BufferUtils.createFloatBuffer(2 * 3);
        float h = (float)Math.sqrt(3.0);
        vb.put(-1.0f).put(-h/3);
        vb.put(+1.0f).put(-h/3);
        vb.put(0.0f).put(2.0f*h/3.0f);
        vb.flip();

        FloatBuffer uvb = BufferUtils.createFloatBuffer(2 * 3);
        uvb.put(0.0f).put(1.0f);
        uvb.put(1.0f).put(1.0f);
        uvb.put(0.5f).put(0.0f);
        uvb.flip();

        IntBuffer ib = BufferUtils.createIntBuffer(6);
        ib.put(0).put(1).put(2);
        ib.flip();
        super.setupVertexArray(pipeline, vb, uvb, ib);
    }
}
