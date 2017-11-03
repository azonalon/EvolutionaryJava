package renderable;
import shader.ShaderPipeline;
import java.nio.*;
import renderable.*;
import static org.lwjgl.system.MemoryUtil.*;
import org.lwjgl.BufferUtils;

public class Square extends Renderable {
    public Square (ShaderPipeline pipeline) {
        FloatBuffer vb = BufferUtils.createFloatBuffer(2 * 6);
        vb.put(-1.0f).put(-1.0f);
        vb.put(1.0f).put(-1.0f);
        vb.put(1.0f).put(1.0f);
        vb.put(1.0f).put(1.0f);
        vb.put(-1.0f).put(1.0f);
        vb.put(-1.0f).put(-1.0f);
        vb.flip();

        FloatBuffer uvb = BufferUtils.createFloatBuffer(2 * 6);
        uvb.put(0.0f).put(1.0f);
        uvb.put(1.0f).put(1.0f);
        uvb.put(1.0f).put(0.0f);
        uvb.put(1.0f).put(0.0f);
        uvb.put(0.0f).put(0.0f);
        uvb.put(0.0f).put(1.0f);
        uvb.flip();

        IntBuffer ib = BufferUtils.createIntBuffer(6);
        ib.put(0).put(1).put(2).put(3).put(4).put(5);
        ib.flip();
        super.setupVertexArray(pipeline, vb, uvb, ib);
    }
}
