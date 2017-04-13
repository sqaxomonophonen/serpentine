precision mediump float;

varying vec3 v_normal;
varying vec2 v_uv;

void main()
{
	vec3 normal_rgb = v_normal * 0.5 + 0.5;
	float border = 0.1;
	vec3 c = normal_rgb * ((fract(v_uv.x) < border || fract(v_uv.x) > (1.0-border) || fract(v_uv.y) < border || fract(v_uv.y) > (1.0-border)) ? 1.0 : 0.2);
	gl_FragColor = vec4(c, 1);
}
