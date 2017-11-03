import java.io.IOException;
import java.nio.*;
import org.joml.Matrix4f;
import org.lwjgl.*;
import org.lwjgl.BufferUtils;
import org.lwjgl.glfw.*;
import org.lwjgl.glfw.GLFWWindowCloseCallbackI;
import org.lwjgl.glfw.GLFWWindowSizeCallback;
import org.lwjgl.opengl.*;
import org.lwjgl.system.*;
import renderable.*;
import shader.*;
import shader.ShaderPipeline;
import view.Window;
import view.View;
import static org.lwjgl.glfw.Callbacks.*;
import static org.lwjgl.glfw.GLFW.*;
import static org.lwjgl.opengl.GL11.*;
import static org.lwjgl.opengl.GL30.*;
import static org.lwjgl.opengl.GL45.*;
import static org.lwjgl.system.MemoryStack.*;
import static org.lwjgl.system.MemoryUtil.*;

public class Main {
    // The window handle
    static float PI = 3.414f;
    private Window window;
    private View view;
    ShaderPipeline pipeline;

	public void run() {
		System.out.println("Hello LWJGL " + Version.getVersion() + "!");

        window = new Window();
		// This line is critical for LWJGL's interoperation with GLFW's
		// OpenGL context, or any context that is managed externally.
		// LWJGL detects the context that is current in the current thread,
		// creates the GLCapabilities instance and makes the OpenGL
		// bindings available for use.
        GL.createCapabilities();
        try {
            initGraphics();
        } catch(IOException exc) {
            System.out.println("Could not load shader");
            System.exit(1);
        }
        view = new View(pipeline, window);
		loop();

	}

    private void initGraphics() throws IOException {
        pipeline = new ShaderPipeline();
        glEnable(GL_TEXTURE_2D);
    }



	private void loop() {

		// Set the clear color
		glClearColor(1.0f, 0.0f, 0.0f, 0.0f);


        Square s = new Square(pipeline);
        Triangle t = new Triangle(pipeline);
        t.relocate(2, 0, 0);

		// Run the rendering loop until the user has attempted to close
		// the window or has pressed the ESCAPE key.
		while ( !glfwWindowShouldClose(window.ID) ) {
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear the framebuffer

            s.draw();
            s.rotateZ(.01f);
            t.draw();
            t.rotateZ(.01f);

            view.updateView();
			glfwSwapBuffers(window.ID); // swap the color buffers
			// Poll for window events. The key callback above will only be
			// invoked during this call.
			glfwPollEvents();
		}
	}

	public static void main(String[] args) {
		new Main().run();
	}

}
