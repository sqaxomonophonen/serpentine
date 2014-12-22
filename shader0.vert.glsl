#version 130

attribute vec3 a_position;
attribute vec3 a_normal;
attribute vec2 a_uv;

uniform mat4 u_projection;
uniform mat4 u_view;

varying vec3 v_normal;
varying vec2 v_uv;

void main()
{
	v_normal = a_normal;
	v_uv = a_uv;
	gl_Position = u_projection * u_view * vec4(a_position, 1);
}
