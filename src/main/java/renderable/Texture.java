package renderable;
import static org.lwjgl.opengl.GL11.*;
import static org.lwjgl.opengl.GL12.*;
import static org.lwjgl.opengl.GL13.*;
import static org.lwjgl.opengl.GL30.*;
import org.lwjgl.BufferUtils;
import javax.imageio.ImageIO;
import java.awt.image.*;
import java.awt.image.Raster;
import java.io.*;
import java.nio.*;
import java.io.IOException;
import java.net.URL;



public class Texture {
    byte[] image;
    int ID;
    // static int format = GL_RGBA8I;

    void genTexture() {
            //Generally a good idea to enable texturing first
            //generate a texture handle or unique ID for this texture
            ID = glGenTextures();
            //bind the texture
            glBindTexture(GL_TEXTURE_2D, ID);
            glActiveTexture(GL_TEXTURE0);
            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        }

    static Texture testRectangles() {
        Texture texture = new Texture();
        texture.genTexture();
        setTextureParameters();
        short[] tex = new short[] {
             0xFFFF/2, 0xffff/2, 0xFFFF/2, 0xFFFF/2,
             0xFFFF/2, 0xffFF/2, 0xFFFF/2, 0xFFFF/2,
             0xFFFF/2, 0xffFF/2, 0xFFFF/2, 0xFFFF/2,
             0xFFFF/2, 0xffFF/2, 0xFFFF/2, 0xFFFF/2,
        };
        glTexImage2D(GL_TEXTURE_2D,  0, GL_RGBA,
                     2, 2,
                     0, GL_RGBA, GL_SHORT, tex);
        return texture;
    }

    static void setTextureParameters() {
        //Setup filtering, i.e. how OpenGL will interpolate the pixels when scaling up or down
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        //Setup wrap mode, i.e. how OpenGL will handle pixels outside of the expected range
        //Note: GL_CLAMP_TO_EDGE is part of GL12
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    }

    static public void main(String[] args) throws IOException { };

    static public Texture fromResource(String resource) throws IOException {
        URL url = Thread.currentThread().getContextClassLoader().getResource(resource);
        File file = new File(url.getFile());
        return fromFile(file);
    }
    static public Texture fromFile(File f) throws IOException {
        Texture texture = new Texture();
        texture.genTexture();
        setTextureParameters();

        BufferedImage image = ImageIO.read(f);
        int w = image.getWidth(), h = image.getHeight();
        int BYTES_PER_PIXEL=4;
        ByteBuffer buffer = BufferUtils.createByteBuffer(w * h * BYTES_PER_PIXEL);

        for(int y = 0; y < h; y++)
        {
            for(int x = 0; x < w; x++)
            {
                int pixel = image.getRGB(x, y);
                buffer.put((byte) ((pixel >> 16) & 0xFF));    // Red component
                buffer.put((byte) ((pixel >> 8 ) & 0xFF));    // Green component
                buffer.put((byte) ((pixel >> 0 ) & 0xFF));    // Blue component
                buffer.put((byte) ((pixel >> 24) & 0xFF));    // Alpha component. Only for RGBA
            }
        }

        buffer.flip();

        glTexImage2D(GL_TEXTURE_2D,  0, GL_RGBA,
                     image.getWidth(), image.getHeight(),
                     0, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
        glBindTexture(GL_TEXTURE_2D, 0);
        return texture;
    }

}
