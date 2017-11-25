#version 330

uniform mat4 modelTransformation;
uniform mat4 viewTransformation;

in vec3 position;
// in vec4 color;
in vec2 textureCoordinate;

out vec4 vertexColor;
out vec2 vertexTextureCoordinate;

void main() {
  gl_Position = viewTransformation * modelTransformation * vec4(position, 1.0);
  // vertexColor = color;
  vertexTextureCoordinate = textureCoordinate;
}
