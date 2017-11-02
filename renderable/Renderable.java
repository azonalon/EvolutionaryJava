package renderable;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import org.lwjgl.*;
import org.lwjgl.BufferUtils;
import org.lwjgl.opengl.*;
import shader.*;
import static org.lwjgl.opengl.GL11.*;
import static org.lwjgl.opengl.GL15.*;
import static org.lwjgl.opengl.GL20.*;
import static org.lwjgl.opengl.GL30.*;
import static org.lwjgl.system.MemoryUtil.*;
import static org.lwjgl.opengl.GL15.*; // glGenBuffers
import static org.lwjgl.opengl.GL20.*; // glVertexAttribPointer
import static org.lwjgl.opengl.GL30.*; // Vertex array objects

import org.joml.Matrix4f;


public class Renderable {
    FloatBuffer vertices; // Edges of polygons
    FloatBuffer colors; // Edges of polygons
    Matrix4f modelMatrix;
    ShaderPipeline pipeline;
    int modelMatrixLocation;
    int uvs; // Texture coordinates
    IntBuffer indices; // indices for indexed drawing
    int positionVbo;
    int textureCoordinateVbo;
    int vaoId; // id of corresponding vertex array object

    public Renderable(ShaderPipeline shader,
               FloatBuffer vertices,  IntBuffer indices) {
        this.vertices = vertices;
        this.indices = indices;
        this.uvs = uvs;
        this.pipeline = shader;
        modelMatrix = new Matrix4f();


        vaoId = glGenVertexArrays();
        glBindVertexArray(vaoId);
        positionVbo = glGenBuffers();
        FloatBuffer fb = BufferUtils.createFloatBuffer(2 * 6);
        fb.put(-1.0f).put(-1.0f);
        fb.put(1.0f).put(-1.0f);
        fb.put(1.0f).put(1.0f);
        fb.put(1.0f).put(1.0f);
        fb.put(-1.0f).put(1.0f);
        fb.put(-1.0f).put(-1.0f);
        fb.flip();
        glBindBuffer(GL_ARRAY_BUFFER, positionVbo);
        glBufferData(GL_ARRAY_BUFFER, fb, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, positionVbo);
        glVertexAttribPointer(shader.positionAttributeId, 2,
                                GL_FLOAT, false, 0, 0L);
        glEnableVertexAttribArray(shader.positionAttributeId);

        textureCoordinateVbo = glGenBuffers();
        fb = BufferUtils.createFloatBuffer(2 * 6);
        fb.put(0.0f).put(1.0f);
        fb.put(1.0f).put(1.0f);
        fb.put(1.0f).put(0.0f);
        fb.put(1.0f).put(0.0f);
        fb.put(0.0f).put(0.0f);
        fb.put(0.0f).put(1.0f);
        fb.flip();
        glBindBuffer(GL_ARRAY_BUFFER, textureCoordinateVbo);
        glBufferData(GL_ARRAY_BUFFER, fb, GL_STATIC_DRAW);
        glVertexAttribPointer(shader.textureCoordinateAttributeId, 2,
                             GL_FLOAT, true, 0, 0L);
        glEnableVertexAttribArray(shader.textureCoordinateAttributeId);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    public void relocate(float x, float y, float z) {
        modelMatrix.translate(x, y, z).get(modelMatrix);
    }
    public void rotateZ(float angle) {
        modelMatrix.rotateZ(angle).get(modelMatrix);
    }

    public void draw() {
        glBindVertexArray(vaoId);
        pipeline.setModelMatrix(modelMatrix);
        glDrawArrays(GL_TRIANGLES, 0, 2*6);
        glBindVertexArray(0);
    }
}
