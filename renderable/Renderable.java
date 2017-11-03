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
    IntBuffer indices; // indices for indexed drawing
    Matrix4f modelMatrix;
    ShaderPipeline pipeline;
    int modelMatrixLocation;
    int uvs; // Texture coordinates
    int positionVbo;
    int textureCoordinateVbo;
    int vaoId; // id of corresponding vertex array object
    int indexVbo;

    protected void setupVertexArray(ShaderPipeline shader,
               FloatBuffer positions,  FloatBuffer textureCoordinates,
               IntBuffer indices) {

        this.vertices = vertices;
        this.indices = indices;
        this.uvs = uvs;
        this.pipeline = shader;
        modelMatrix = new Matrix4f();
        vaoId = glGenVertexArrays();
        glBindVertexArray(vaoId);
        positionVbo = glGenBuffers();
        glBindBuffer(GL_ARRAY_BUFFER, positionVbo);
        glBufferData(GL_ARRAY_BUFFER, positions, GL_STATIC_DRAW);
        glVertexAttribPointer(shader.positionAttributeId, 2,
                                GL_FLOAT, false, 0, 0L);
        glEnableVertexAttribArray(shader.positionAttributeId);


        textureCoordinateVbo = glGenBuffers();
        glBindBuffer(GL_ARRAY_BUFFER, textureCoordinateVbo);
        glBufferData(GL_ARRAY_BUFFER, textureCoordinates, GL_STATIC_DRAW);
        glVertexAttribPointer(shader.textureCoordinateAttributeId, 2,
                             GL_FLOAT, true, 0, 0L);
        glEnableVertexAttribArray(shader.textureCoordinateAttributeId);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        indexVbo = glGenBuffers();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVbo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices, GL_STATIC_DRAW);
        glBindVertexArray(0);
    }

    public void relocate(float x, float y, float z) {
        modelMatrix.translate(x, y, z).get(modelMatrix);
    }
    public void rotateZ(float angle) {
        modelMatrix.rotateZ(angle).get(modelMatrix);
    }

    public void draw() {
        pipeline.setModelMatrix(modelMatrix);
        glBindVertexArray(vaoId);
        glDrawElements(GL_TRIANGLES, indices.capacity(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
    }
}
