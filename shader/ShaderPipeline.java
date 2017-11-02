package shader;
import static org.lwjgl.opengl.GL20.*;
import java.io.*;
import java.nio.*;
import shader.*;
import org.joml.Matrix4f;
import org.lwjgl.BufferUtils;

public class ShaderPipeline {
    public int positionAttributeId;
    public int textureCoordinateAttributeId;
    public int textureLocation;
    public int modelTransformationLocation;
    public int viewTransformationLocation;
    private FloatBuffer matrixBuffer;

    public ShaderPipeline() throws IOException {
        matrixBuffer = BufferUtils.createFloatBuffer(16);
        createShaders();
    }
    void createShaders() throws IOException {
        int program = glCreateProgram();
        int vshader = Shader.createShader("./shader/pipeline.vs", GL_VERTEX_SHADER);
        int fshader = Shader.createShader("./shader/pipeline.fs", GL_FRAGMENT_SHADER);
        glAttachShader(program, vshader);
        glAttachShader(program, fshader);
        glLinkProgram(program);
        int linked = glGetProgrami(program, GL_LINK_STATUS);
        String programLog = glGetProgramInfoLog(program);
        if (programLog != null && programLog.trim().length() > 0)
        System.err.println(programLog);
        if (linked == 0)
        throw new AssertionError("Could not link program");
        glUseProgram(program);
        textureLocation = glGetUniformLocation(program, "tex");
        glUniform1i(textureLocation, 0);
        positionAttributeId = glGetAttribLocation(program, "position");
        textureCoordinateAttributeId = glGetAttribLocation(program, "textureCoordinate");
        textureCoordinateAttributeId = glGetAttribLocation(program, "color");
        modelTransformationLocation = glGetUniformLocation(program, "modelTransformation");
        viewTransformationLocation = glGetUniformLocation(program, "viewTransformation");

        // initialize mvp
        new Matrix4f().get(matrixBuffer);
        glUniformMatrix4fv(modelTransformationLocation, false, matrixBuffer);
        glUniformMatrix4fv(viewTransformationLocation, false, matrixBuffer);
        // glUseProgram(0);
    }
    public void setModelMatrix(Matrix4f m) {
        m.get(matrixBuffer);
        glUniformMatrix4fv(modelTransformationLocation, false, matrixBuffer);
    }
    public void setViewMatrix(Matrix4f m) {
        m.get(matrixBuffer);
        glUniformMatrix4fv(viewTransformationLocation, false, matrixBuffer);
    }
};
