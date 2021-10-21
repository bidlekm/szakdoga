#shader vertex
#version 330 core
layout(location = 0) in vec2 pos;
out vec2 texcoord; // specify a color output to the fragment shader

void main()
{
	texcoord = (pos + vec2(1, 1)) / 2.0;
	gl_Position = vec4(pos.x, pos.y ,0.0, 1.0); 
}



#shader fragment
#version 330
precision highp float;

uniform sampler2D textureUnit;
in  vec2 texcoord;			
out vec4 fragmentColor;	

void main() {
	fragmentColor = texture(textureUnit, texcoord);
}