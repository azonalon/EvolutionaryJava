package renderable;
import static org.lwjgl.opengl.GL11.*;
import static org.lwjgl.opengl.GL12.*;
import static org.lwjgl.opengl.GL13.*;
import org.lwjgl.BufferUtils;
import javax.imageio.ImageIO;
import java.awt.image.*;
import java.awt.image.Raster;
import java.io.*;
import java.nio.*;
import java.io.IOException;



public class Texture {
    byte[] image;
    int ID;
    // static int format = GL_RGBA8I;
    Texture() {};
    static Texture fromPng() throws IOException {
        BufferedImage im = ImageIO.read(new File("./rsc/textures/triangleEye.png"));
        byte[] pixels = ((DataBufferByte) im.getRaster().getDataBuffer()).getData();
        // byte[] pixels = ((DataBufferByte) im.getRaster().getDataBuffer()).getData();
        Texture texture = new Texture();
        //Generally a good idea to enable texturing first

        //generate a texture handle or unique ID for this texture
        texture.ID = glGenTextures();

        //bind the texture
        glBindTexture(GL_TEXTURE_2D, texture.ID);
        glActiveTexture(GL_TEXTURE0);

        //use an alignment of 1 to be safe
        //this tells OpenGL how to unpack the RGBA bytes we will specify
        // TODO: find out if necessary???
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        setTextureParameters();
        //upload our ByteBuffer to GL
        ByteBuffer buf = BufferUtils.createByteBuffer(pixels.length);
        buf.put(pixels);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, im.getWidth(),
                     im.getHeight(), 0, GL_RGBA8, GL_UNSIGNED_BYTE, buf);

        // glBindTexture(0);
        // System.out.println(pixels.length);

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

    static public void main(String[] args) throws IOException{
        fromPng();
    }

    static void readImageWrong() throws IOException {
        BufferedImage im = ImageIO.read(new File("./rsc/textures/test.png"));
        int w = im.getWidth();
        int h = im.getHeight();
        Raster r = im.getAlphaRaster();
        byte[] a = new byte[w*h*4];
        for(int i = 0; i < w; i++) {
            for(int j = 0; j < h; j++) {
                int rgb = im.getRGB(i, j);
                int red = (rgb >> 16) & 0xFF;
                int green = (rgb >> 8) & 0xFF;
                int blue = rgb & 0xFF;
                float falpha = r.getDataBuffer().getElemFloat(w * i + j);
                byte alpha = (byte)Math.round(falpha * 255 - 127);
                a[w * i * 4 + j * 4 + 0] = (byte)red;
                a[w * i * 4 + j * 4 + 1] = (byte)green;
                a[w * i * 4 + j * 4 + 2] = (byte)blue;
                a[w * i * 4 + j * 4 + 3] = alpha;
                // System.out.format("%d %d %d %d", red, green, blue, alpha);
                // System.out.println("Alpha value:" + falpha);
            }
            // System.out.println();
        }
    }
}
