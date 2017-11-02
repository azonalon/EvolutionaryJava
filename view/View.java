package view;
import shader.ShaderPipeline;
import view.Window;
import org.joml.Matrix4f;

// import static org.lwjgl.glfw.Callbacks.*;
import static org.lwjgl.glfw.GLFW.*;

// Handler for mouse controls and view matrix etc...
public class View
{
    ShaderPipeline pipeline;
    Window window;
    double zoom = 1;
    double x0=0, y0=0;
    boolean dragMode = false;
    double dragX=0, dragY=0;

    public View(ShaderPipeline pipeline, Window  window) {
        this.pipeline = pipeline;
        this.window = window;
        // glfwSetInputMode(window.ID , GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
		glfwSetScrollCallback(window.ID, (ID, x, y) -> {
            y = -y;
            if(zoom + y > 0) {
                zoom += y;
            }
            System.out.println("Scroll!: " + x +":" +  y);
        });
        glfwSetCursorPosCallback(window.ID, (ID, xpos, ypos) -> {
            if(dragMode == true && glfwGetMouseButton(window.ID, GLFW_MOUSE_BUTTON_LEFT) ==
                                    GLFW_PRESS) {
                double dx = xpos - dragX;
                double dy = ypos - dragY;
                x0 += dx*0.001;
                y0 -= dy*0.001;
                updateDragCoordinates();
            } else {
                dragMode = false;
            }
        });
        glfwSetMouseButtonCallback(window.ID, (ID, button, action, mods) -> {
            if(action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_LEFT) {
                updateDragCoordinates();
            } else if(action == GLFW_RELEASE && button == GLFW_MOUSE_BUTTON_LEFT) {
                dragMode = false;
            }
        });
    }

    void updateDragCoordinates() {
        double xpos[] = new double[1];
        double ypos[] = new double[1];
        glfwGetCursorPos(window.ID, xpos, ypos);
        dragMode = true;
        dragX = xpos[0];
        dragY = ypos[0];
    }

    public void updateView() {
        // glfwSetCursorPos(window.ID , window.width/2, window.height/2);
        updateViewMatrix();
    }
    private void updateViewMatrix() {
        double aspect = window.getAspect();
        double x0 = -zoom  + this.x0;
        double x1 = zoom + this.x0;
        double y0 = -zoom + this.y0;
        double y1 = +zoom + this.y0;
        if(aspect >= 1) {
            x0 *= aspect;
            x1 *= aspect;
        } else {
            y0 *= 1/aspect;
            y1 *= 1/aspect;
        }
        Matrix4f viewTransformation = new Matrix4f()
                // .perspective(PI/8, aspect, 0, -20)
                .ortho((float)x0, (float)x1, (float)y0, (float)y1, 0, 2);
                // .lookAt(0, 0, 1, 0, 0, 0, 1, 0, 0);
        pipeline.setViewMatrix(viewTransformation);
    }
}
