#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "a.h"
#include "m.h"

#include <SDL.h>
#include <GL/glew.h>

void glew_init()
{
	GLenum err = glewInit();
	if (err != GLEW_OK) {
		arghf("glewInit() failed: %s", glewGetErrorString(err));
	}

#define CHECK_GL_EXT(x) { if(!GLEW_ ## x) arghf("OpenGL extension not found: " #x); }
	CHECK_GL_EXT(ARB_shader_objects);
	CHECK_GL_EXT(ARB_vertex_shader);
	CHECK_GL_EXT(ARB_fragment_shader);
	CHECK_GL_EXT(ARB_framebuffer_object);
	CHECK_GL_EXT(ARB_vertex_buffer_object);
#undef CHECK_GL_EXT

	/* to figure out what extension something belongs to, see:
	 * http://www.opengl.org/registry/#specfiles */

	// XXX check that version is at least 1.30?
	printf("GLSL version %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

	GLint value;
#define DUMP_GL_INT(x) { glGetIntegerv(x, &value); printf(#x " = %d\n", value); }
	DUMP_GL_INT(GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS_ARB);
	DUMP_GL_INT(GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS);
	DUMP_GL_INT(GL_MAX_TEXTURE_IMAGE_UNITS);
	DUMP_GL_INT(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS);
	DUMP_GL_INT(GL_MAX_CLIP_PLANES);
	DUMP_GL_INT(GL_MAX_FRAGMENT_UNIFORM_COMPONENTS);
	DUMP_GL_INT(GL_MAX_VERTEX_UNIFORM_COMPONENTS);
	DUMP_GL_INT(GL_MAX_VARYING_FLOATS);
	DUMP_GL_INT(GL_MAX_VERTEX_ATTRIBS);
	DUMP_GL_INT(GL_MAX_DRAW_BUFFERS);
	DUMP_GL_INT(GL_STENCIL_BITS);
#undef DUMP_GL_INT
}


struct vec3i { int32_t s[3]; };

int vec3i_equals(struct vec3i* a, struct vec3i* b)
{
	for (int i = 0; i < 3; i++) {
		if (a->s[i] != b->s[i]) {
			return 0;
		}
	}
	return 1;
}


#define SHADER_MAX_ATTRS (16)

enum shader_attr_type {
	SHADER_ATTR_FLOAT = 1,
	SHADER_ATTR_VEC2,
	SHADER_ATTR_VEC3,
	SHADER_ATTR_VEC4,
};

int shader_attr_type_floats(enum shader_attr_type type)
{
	switch (type) {
		case SHADER_ATTR_FLOAT: return 1;
		case SHADER_ATTR_VEC2: return 2;
		case SHADER_ATTR_VEC3: return 3;
		case SHADER_ATTR_VEC4: return 4;
	}
}

struct shader_attr_spec {
	const char* symbol;
	enum shader_attr_type type;
};

struct shader {
	GLuint program;
	int n_attrs;
	GLuint attr_locations[SHADER_MAX_ATTRS];
	enum shader_attr_type attr_types[SHADER_MAX_ATTRS];
	size_t stride;
};

static GLuint create_shader(GLenum type, const char* src)
{
	GLuint shader = glCreateShader(type); CHKGL;
	glShaderSource(shader, 1, &src, 0);
	glCompileShader(shader);

	GLint status;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE) {
		GLint msglen;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &msglen);
		GLchar* msg = (GLchar*) malloc(msglen + 1);
		glGetShaderInfoLog(shader, msglen, NULL, msg);
		const char* stype = type == GL_VERTEX_SHADER ? "vertex" : type == GL_FRAGMENT_SHADER ? "fragment" : "waaaat";
		arghf("%s shader error: %s -- source:\n%s", stype, msg, src);
	}
	return shader;
}

void shader_init(struct shader* shader, const char* vert_src, const char* frag_src, struct shader_attr_spec* attr_specs)
{
	memset(shader, 0, sizeof(*shader));

	GLuint vertex_shader = create_shader(GL_VERTEX_SHADER, vert_src);
	GLuint fragment_shader = create_shader(GL_FRAGMENT_SHADER, frag_src);

	shader->program = glCreateProgram(); CHKGL;

	glAttachShader(shader->program, vertex_shader);
	glAttachShader(shader->program, fragment_shader);

	glLinkProgram(shader->program);

	GLint status;
	glGetProgramiv(shader->program, GL_LINK_STATUS, &status);
	if (status == GL_FALSE) {
		GLint msglen;
		glGetProgramiv(shader->program, GL_INFO_LOG_LENGTH, &msglen);
		GLchar* msg = (GLchar*) malloc(msglen + 1);
		glGetProgramInfoLog(shader->program, msglen, NULL, msg);
		arghf("shader link error: %s", msg);
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	int i = 0;
	size_t stride = 0;
	for (struct shader_attr_spec* spec = attr_specs; spec->symbol != NULL; spec++, i++) {
		shader->attr_locations[i] = glGetAttribLocation(shader->program, spec->symbol); CHKGL;
		shader->attr_types[i] = spec->type;
		stride += shader_attr_type_floats(shader->attr_types[i]) * sizeof(float);
	}
	shader->n_attrs = i;
	shader->stride = stride;
}

void shader_use(struct shader* shader)
{
	glUseProgram(shader->program);
}

void shader_set_attrib_pointers(struct shader* shader)
{
	char* offset = 0;
	for (int i = 0; i < shader->n_attrs; i++) {
		int n_floats = shader_attr_type_floats(shader->attr_types[i]);
		glVertexAttribPointer(shader->attr_locations[i], n_floats, GL_FLOAT, GL_FALSE, shader->stride, offset); CHKGL;
		offset += n_floats * sizeof(float);
	}
}

void shader_enable_arrays(struct shader* shader)
{
	for (int i = 0; i < shader->n_attrs; i++) {
		glEnableVertexAttribArray(shader->attr_locations[i]); CHKGL;
	}
}

void shader_disable_arrays(struct shader* shader)
{
	for (int i = 0; i < shader->n_attrs; i++) {
		glDisableVertexAttribArray(shader->attr_locations[i]); CHKGL;
	}
}

void shader_uniform_mat44(struct shader* shader, const char* name, struct mat44* m)
{
	GLint location = glGetUniformLocation(shader->program, name); CHKGL;
	glUniformMatrix4fv(location, 1, GL_FALSE, m->s); CHKGL;
}

struct vtxbuf {
	GLuint buffer;
	size_t sz, used;
	float* data;
	struct shader* shader;
	GLenum mode;
};

void vtxbuf_init(struct vtxbuf* vb, size_t sz)
{
	memset(vb, 0, sizeof(*vb));
	vb->sz = sz;
	vb->data = malloc(sz);
	AN(vb->data);
	glGenBuffers(1, &vb->buffer); CHKGL;
	glBindBuffer(GL_ARRAY_BUFFER, vb->buffer); CHKGL;
	glBufferData(GL_ARRAY_BUFFER, sz, vb->data, GL_STREAM_DRAW); CHKGL;
}

void vtxbuf_begin(struct vtxbuf* vb, struct shader* shader, GLenum mode)
{
	vb->shader = shader;
	vb->mode = mode;
	vb->used = 0;
	shader_use(shader);
	shader_enable_arrays(shader);
}

void vtxbuf_flush(struct vtxbuf* vb)
{
	if (vb->used == 0) return;

	glBindBuffer(GL_ARRAY_BUFFER, vb->buffer); CHKGL;
	glBufferSubData(GL_ARRAY_BUFFER, 0, vb->used, vb->data); CHKGL;

	shader_set_attrib_pointers(vb->shader);
	glBindBuffer(GL_ARRAY_BUFFER, vb->buffer); CHKGL;
	glDrawArrays(vb->mode, 0, vb->used / vb->shader->stride);
	vb->used = 0;
}

void vtxbuf_end(struct vtxbuf* vb)
{
	vtxbuf_flush(vb);
	shader_disable_arrays(vb->shader);
}

void vtxbuf_element(struct vtxbuf* vb, float* data, size_t sz)
{
	if ((vb->used + sz) > vb->sz) vtxbuf_flush(vb);
	if ((vb->used + sz) > vb->sz) WRONG("not enough room for even one element");
	memcpy(((uint8_t*)vb->data) + vb->used, data, sz);
	vb->used += sz;
}

static void gl_viewport_from_sdl_window(SDL_Window* window)
{
	int w, h;
	SDL_GetWindowSize(window, &w, &h);
	glViewport(0, 0, w, h); CHKGL;
}

#define WORLD_CHUNK_LENGTH_EXP (3)
#define WORLD_CHUNK_LENGTH (1 << WORLD_CHUNK_LENGTH_EXP)
#define WORLD_CHUNK_VOLUME (1 << (WORLD_CHUNK_LENGTH_EXP * 3))

struct world_cell {
	uint8_t type, _pad0, _pad1, _pad2;
};

struct world_chunk {
	struct world_cell cells[WORLD_CHUNK_VOLUME];
};

static inline int world_chunk_cell_inside(struct world_chunk* chunk, struct vec3i position)
{
	for (int i = 0; i < 3; i++) if (position.s[i] < 0 || position.s[i] >= WORLD_CHUNK_LENGTH) return 0;
	return 1;
}

static inline struct world_cell* world_chunk_cell_atp(struct world_chunk* chunk, struct vec3i position)
{
	int index = 0;
	for (int i = 0; i < 3; i++) index += position.s[i] * (1 << (WORLD_CHUNK_LENGTH_EXP*i));
	return &chunk->cells[index];
}

void world_chunk_init(struct world_chunk* chunk)
{
	memset(chunk, 0, sizeof(*chunk));
	#if 1
	for (int z = 0; z < WORLD_CHUNK_LENGTH; z++)
	for (int y = 0; y < WORLD_CHUNK_LENGTH; y++)
	for (int x = 0; x < WORLD_CHUNK_LENGTH; x++) {
		struct vec3i pos = {{x,y,z}};
		struct world_cell* cell = world_chunk_cell_atp(chunk, pos);
		cell->type = (rand()&7) == 0;
	}
	#endif
}


int world_chunk_draw(struct world_chunk* chunk, struct vtxbuf* vb, struct vec3i offset)
{
	int quads = 0;
	for (int z = 0; z < WORLD_CHUNK_LENGTH; z++)
	for (int y = 0; y < WORLD_CHUNK_LENGTH; y++)
	for (int x = 0; x < WORLD_CHUNK_LENGTH; x++) {
		struct vec3i pos = {{x,y,z}};
		struct world_cell* cell = world_chunk_cell_atp(chunk, pos);
		if (cell->type > 0) {
			for (int axis = 0; axis < 3; axis++) {
				for (int sign = 0; sign < 2; sign++) {
					{
						struct vec3i nbpos = pos;
						nbpos.s[axis] += sign ? 1 : -1;
						if (world_chunk_cell_inside(chunk, nbpos)) {
							struct world_cell* nbcell = world_chunk_cell_atp(chunk, nbpos);
							if (nbcell->type > 0) continue;
						}
					}

					int flip = axis == 1; // XXX wtf is wrong with the y axis?
					float quad[8*4];
					struct vec3 normal = {{0,0,0}};
					normal.s[axis] = sign ? 1 : -1;
					int p = 0;
					for (int i = 0; i < 4; i++) {
						struct vec3 point = {{x,y,z}};
						int ds[2] = {
							(i&1) ^ ((i>>1)&1),
							((i>>1)&1) ^ sign ^ flip ^ 1
						};
						int dsi = 0;
						for (int j = 0; j < 3; j++) {
							float a = 0;
							if (j == axis) {
								a = sign;
							} else {
								a = ds[dsi++];
							}
							quad[p++] = point.s[j] + a + offset.s[j] * WORLD_CHUNK_LENGTH;
						}
						for (int j = 0; j < 3; j++) quad[p++] = normal.s[j];

						float uv[2] = {ds[0], ds[1]};
						for (int j = 0; j < 2; j++) quad[p++] = uv[j];

					}
					vtxbuf_element(vb, quad, sizeof(quad));
					quads++;
				}
			}
		}
	}
	//printf("quads: %d\n", quads);
	return quads;
}

int main(int argc, char** argv)
{
	SAZ(SDL_Init(SDL_INIT_VIDEO));
	atexit(SDL_Quit);

	int bitmask = SDL_WINDOW_FULLSCREEN_DESKTOP | SDL_WINDOW_OPENGL;
	SDL_Window* window = SDL_CreateWindow(
			"??? ???",
			SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
			0, 0,
			bitmask);

	SAN(window);

	SDL_GLContext glctx = SDL_GL_CreateContext(window);
	SAN(glctx);

	SAZ(SDL_GL_SetSwapInterval(1)); // or -1, "late swap tearing"?

	glew_init();

	gl_viewport_from_sdl_window(window);

	struct mat44 projection;
	{
		int width, height;
		SDL_GetWindowSize(window, &width, &height);
		float fovy = 65;
		float aspect = (float)width / (float)height;
		mat44_set_perspective(&projection, fovy, aspect, 0.01f, 409.6f);
	}

	struct vtxbuf vtxbuf;
	vtxbuf_init(&vtxbuf, 65536);

	struct world_chunk chunk;
	world_chunk_init(&chunk);

	struct shader shader0;
	{
		#include "shader0.glsl.inc"
		struct shader_attr_spec specs[] = {
			{"a_position", SHADER_ATTR_VEC3},
			{"a_normal", SHADER_ATTR_VEC3},
			{"a_uv", SHADER_ATTR_VEC2},
			{NULL}
		};
		shader_init(&shader0, shader0_vert_src, shader0_frag_src, specs);
	}

	int ctrl_forward = 0;
	int ctrl_backward = 0;
	int ctrl_left = 0;
	int ctrl_right = 0;

	float yaw = 0;
	float pitch = 0;
	struct vec3 position = {{4,16,4}};

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	SDL_SetRelativeMouseMode(SDL_TRUE);
	int exiting = 0;
	while (!exiting) {
		SDL_Event e;
		int mdx = 0;
		int mdy = 0;
		while (SDL_PollEvent(&e)) {
			 if (e.type == SDL_QUIT) exiting = 1;

			 struct push_key {
				 SDL_Keycode sym;
				 int* intptr;
			 } push_keys[] = {
				 {SDLK_w, &ctrl_forward},
				 {SDLK_s, &ctrl_backward},
				 {SDLK_a, &ctrl_left},
				 {SDLK_d, &ctrl_right},
				 {-1, NULL}
			 };
			 for (struct push_key* tkp = push_keys; tkp->intptr != NULL; tkp++) {
				 if ((e.type == SDL_KEYDOWN || e.type == SDL_KEYUP) && e.key.keysym.sym == tkp->sym) {
					 *(tkp->intptr) = (e.type == SDL_KEYDOWN);
				 }
			 }

			 if (e.type == SDL_KEYDOWN) {
				 if (e.key.keysym.sym == SDLK_ESCAPE) {
					 exiting = 1;
				 }
			 }

			 if (e.type == SDL_MOUSEMOTION) {
				 mdx += e.motion.xrel;
				 mdy += e.motion.yrel;
			 }
		}

		{
			float sensitivity = 0.1f;
			yaw += (float)mdx * sensitivity;
			pitch += (float)mdy * sensitivity;
			float pitch_limit = 90;
			if (pitch > pitch_limit) pitch = pitch_limit;
			if (pitch < -pitch_limit) pitch = -pitch_limit;
		}

		{
			float speed = 0.1f;
			float forward = (float)(ctrl_forward - ctrl_backward) * speed;
			float right = (float)(ctrl_right - ctrl_left) * speed;
			struct vec3 movement;
			vec3_move(&movement, yaw, pitch, forward, right);
			vec3_add_inplace(&position, &movement);
		}

		struct mat44 view;
		mat44_set_identity(&view);
		mat44_rotate_x(&view, pitch);
		mat44_rotate_y(&view, yaw);

		struct vec3 translate;
		vec3_scale(&translate, &position, -1);
		mat44_translate(&view, &translate);

		glClearColor(0.0, 0.0, 0.2, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		vtxbuf_begin(&vtxbuf, &shader0, GL_QUADS);
		shader_uniform_mat44(&shader0, "u_projection", &projection);
		shader_uniform_mat44(&shader0, "u_view", &view);

		int min = -1;
		int max = 1;
		int quads = 0;
		for (int dz = min; dz <= max; dz++)
		for (int dy = min; dy <= max; dy++)
		for (int dx = min; dx <= max; dx++) {
			struct vec3i offset = {{dx,dy,dz}};
			quads += world_chunk_draw(&chunk, &vtxbuf, offset);
		}
		printf("quads: %d\n", quads);

		vtxbuf_end(&vtxbuf);

		SDL_GL_SwapWindow(window);
	}

	SDL_DestroyWindow(window);
	SDL_GL_DeleteContext(glctx);

	return EXIT_SUCCESS;
}

