#version 330

uniform sampler2D textureSampler;

// in vec4 vertexColor;
in vec2 vertexTextureCoordinate;

out vec4 color;

void main() {
  // color = vertexColor;
  // color = vec4(0,1,1,1);
  color = texture(textureSampler, vertexTextureCoordinate) * 0.00001 + vec4(0,1,1,1);
}
