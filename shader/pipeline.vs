#version 130

in vec3 position;
in vec2 textureCoordinate;
in vec4 color;
out vec4 fragmentColor;
uniform mat4 modelTransformation;
uniform mat4 viewTransformation;

void main() {
  gl_Position = viewTransformation * modelTransformation * vec4(position, 1.0);
  fragmentColor  = color;
}
