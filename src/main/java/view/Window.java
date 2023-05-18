package view;
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
import static org.lwjgl.glfw.Callbacks.*;
import static org.lwjgl.glfw.GLFW.*;
import static org.lwjgl.opengl.GL11.*;
import static org.lwjgl.opengl.GL30.*;
import static org.lwjgl.system.MemoryStack.*;
import static org.lwjgl.system.MemoryUtil.*;

public class Window {
    public long ID;
    public int width, height;
    public Window() {
        initWindow();
    }
    public float getAspect() {
        return height != 0 ? (float)width/(float)height : 1;
    }
    private void setCoreVersion() {
    	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,4);// use opengl 4.5
    	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,5);
    	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    }
	private void initWindow() {
		// Setup an error callback. The default implementation
		// will print the error message in System.err.
		GLFWErrorCallback.createPrint(System.err).set();

		// Initialize GLFW. Most GLFW functions will not work before doing this.
		if ( !glfwInit() )
			throw new IllegalStateException("Unable to initialize GLFW");
		glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE); // the window will stay hidden after creation
		glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE); // the window will be resizable

		// Configure GLFW
		// glfwDefaultWindowHints(); // optional, the current window hints are already the default

		// Create the window
        setCoreVersion();
		ID = glfwCreateWindow(300, 300, "Hello World!", NULL, NULL);
		if ( ID == NULL )
			throw new RuntimeException("Failed to create the GLFW window");

		// Setup a key callback. It will be called every time a key is pressed, repeated or released.
		glfwSetKeyCallback(ID, (ID, key, scancode, action, mods) -> {
			if ( key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE )
				glfwSetWindowShouldClose(ID, true); // We will detect this in the rendering loop
		});

		// Get the thread stack and push a new frame
		try ( MemoryStack stack = stackPush() ) {
			IntBuffer pWidth = stack.mallocInt(1); // int*
			IntBuffer pHeight = stack.mallocInt(1); // int*

			// Get the window size passed to glfwCreateWindow
			glfwGetWindowSize(ID, pWidth, pHeight);

			// Get the resolution of the primary monitor
			// GLFWVidMode vidmode = glfwGetVideoMode(glfwGetPrimaryMonitor());

			// Center the window
			// glfwSetWindowPos(
			// 	ID,
			// 	(vidmode.width() - pWidth.get(0)) / 2,
			// 	(vidmode.height() - pHeight.get(0)) / 2
			// );
		} // the stack frame is popped automatically

        // GLFWWindowSizeCallback cb = GLFWWindowSizeCallback.create(onWindowSizeChanged);
        glfwSetWindowSizeCallback(ID, new GLFWWindowSizeCallback() {
            @Override
            public void invoke(long window, int w, int h) {
                // System.out.println("Window width and height " + w + "x" + h);
                // glfwSetWindowSize(window, w, h);
                glViewport(0, 0, w, h);
                width = w;
                height = h;
            }
        });

		// Make the OpenGL context current
		glfwMakeContextCurrent(ID);
		// Enable v-sync
		glfwSwapInterval(1);

		// Make the ID visible
		glfwShowWindow(ID);
	}

    // private void terminateWindow() {
	// 	// Free the window callbacks and destroy the window
	// 	glfwFreeCallbacks(ID);
	// 	glfwDestroyWindow(ID);
    //
	// 	// Terminate GLFW and free the error callback
	// 	glfwTerminate();
	// 	glfwSetErrorCallback(null).free();
    // }
}
