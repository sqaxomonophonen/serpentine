#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "a.h"
#include "m.h"

#include <SDL.h>
#include <GL/glew.h>

#define GRAVITY_AXIS (1)
#define GRAVITY_SIGN (-1)
#define GRAVITY_FORCE (0.18)
#define JUMP_HEIGHT (1.1)

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

void vec3i_dump(struct vec3i* v)
{
	printf("(%d  %d  %d)\n", v->s[0], v->s[1], v->s[2]);
}

struct vec3i_it {
	struct vec3i current, min, max;
	int more;
};

void vec3i_it_init(struct vec3i_it* it, struct vec3i* min, struct vec3i* max)
{
	it->current = *min;
	it->min = *min;
	it->max = *max;
	it->more = 1;
}

void vec3i_it_next(struct vec3i_it* it)
{
	AN(it->more);
	int axis = 0;
	while (axis < 3) {
		it->current.s[axis]++;
		if (it->current.s[axis] > it->max.s[axis]) {
			it->current.s[axis] = it->min.s[axis];
			axis++;
		} else {
			break;
		}
	}
	if (axis == 3) it->more = 0;
}


struct aabb {
	struct vec3 min, max;
};

void aabb_init_min_max(struct aabb* aabb, struct vec3* min, struct vec3* max)
{
	memcpy(&aabb->min, min, sizeof(*min));
	memcpy(&aabb->max, max, sizeof(*max));
}

void aabb_init_around(struct aabb* aabb, struct vec3* pos, struct vec3* min, struct vec3* max)
{
	vec3_add(&aabb->min, pos, min);
	vec3_add(&aabb->max, pos, max);
}

void aabb_grow_inplace(struct aabb* aabb, float s)
{
	for (int axis = 0; axis < 3; axis++) {
		aabb->min.s[axis] -= s;
		aabb->max.s[axis] += s;
	}
}

void aabb_extend_inplace(struct aabb* aabb, struct vec3* x)
{
	for (int axis = 0; axis < 3; axis++) {
		float v = x->s[axis];
		if (v < 0) {
			aabb->min.s[axis] += v;
		} else {
			aabb->max.s[axis] += v;
		}
	}
}

int aabb_aabb_intersection(struct aabb* a, struct aabb* b)
{
	for (int axis = 0; axis < 3; axis++) {
		if (a->max.s[axis] < b->min.s[axis] || a->min.s[axis] > b->max.s[axis]) return 0;
	}
	return 1;
}

void aabb_cube_at(struct aabb* aabb, struct vec3i* pos)
{
	struct vec3 min = {{pos->s[0], pos->s[1], pos->s[2]}};
	aabb->min = min;
	struct vec3 max = {{min.s[0]+1, min.s[1]+1, min.s[2]+1}};
	aabb->max = max;
}

void aabb_aabb_convolute(struct aabb* dst, struct aabb* a, struct aabb* b)
{
	for (int axis = 0; axis < 3; axis++) {
		dst->min.s[axis] = b->min.s[axis] - a->max.s[axis];
		dst->max.s[axis] = b->max.s[axis] - a->min.s[axis];
	}
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
	uint8_t type;
	uint8_t _pad0, _pad1, _pad2;
};

struct world_chunk {
	struct world_cell cells[WORLD_CHUNK_VOLUME];
};

static inline int world_chunk_cell_inside(struct world_chunk* chunk, struct vec3i position)
{
	for (int i = 0; i < 3; i++) if (position.s[i] < 0 || position.s[i] >= WORLD_CHUNK_LENGTH) return 0;
	return 1;
}

static inline struct world_cell* world_chunk_cell_atp(struct world_chunk* chunk, struct vec3i* position)
{
	int index = 0;
	int shift = 0;
	for (int i = 0; i < 3; i++) {
		index += position->s[i] * (1 << shift);
		shift += WORLD_CHUNK_LENGTH_EXP;
	}
	return &chunk->cells[index];
}

int world_chunk_draw(struct world_chunk* chunk, struct vtxbuf* vb, struct vec3i offset)
{
	int quads = 0;
	for (int z = 0; z < WORLD_CHUNK_LENGTH; z++)
	for (int y = 0; y < WORLD_CHUNK_LENGTH; y++)
	for (int x = 0; x < WORLD_CHUNK_LENGTH; x++) {
		struct vec3i pos = {{x,y,z}};
		struct world_cell* cell = world_chunk_cell_atp(chunk, &pos);
		if (cell->type > 0) {
			for (int axis = 0; axis < 3; axis++) {
				for (int sign = 0; sign < 2; sign++) {
					{
						struct vec3i nbpos = pos;
						nbpos.s[axis] += sign ? 1 : -1;
						if (world_chunk_cell_inside(chunk, nbpos)) {
							struct world_cell* nbcell = world_chunk_cell_atp(chunk, &nbpos);
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
	return quads;
}


struct world {
	struct world_chunk* chunks;
};

struct entity {
	struct vec3 position;
	struct vec3 velocity;
	uint32_t on_ground:1;
	uint32_t ctrl_jump:1;
	uint32_t ctrl_forward:1;
	uint32_t ctrl_backward:1;
	uint32_t ctrl_left:1;
	uint32_t ctrl_right:1;
	float pitch, yaw;
};

void entity_init(struct entity* entity)
{
	memset(entity, 0, sizeof(*entity));
}

void entity_get_aabb(struct entity* entity, struct aabb* aabb)
{
	struct vec3 extents_min = {{-0.4, -0.4, -0.4}};
	struct vec3 extents_max = {{0.4, 0.4, 0.4}};
	aabb_init_around(aabb, &entity->position, &extents_min, &extents_max);
}

struct render {
	SDL_Window* window;
	struct mat44 projection;
	struct shader shader0;
	struct vtxbuf vtxbuf;
};

void render_update_projection(struct render* render)
{
	int width, height;
	SDL_GetWindowSize(render->window, &width, &height);

	gl_viewport_from_sdl_window(render->window);

	float fovy = 65;
	float aspect = (float)width / (float)height;
	mat44_set_perspective(&render->projection, fovy, aspect, 0.01f, 409.6f);
}

void render_init(struct render* render, SDL_Window* window)
{
	render->window = window;
	{
		#include "shader0.glsl.inc"
		struct shader_attr_spec specs[] = {
			{"a_position", SHADER_ATTR_VEC3},
			{"a_normal", SHADER_ATTR_VEC3},
			{"a_uv", SHADER_ATTR_VEC2},
			{NULL}
		};
		shader_init(&render->shader0, shader0_vert_src, shader0_frag_src, specs);
	}
	vtxbuf_init(&render->vtxbuf, 65536);
}

void world_draw(struct world* world, struct render* render, struct mat44* view)
{
	render_update_projection(render);

	vtxbuf_begin(&render->vtxbuf, &render->shader0, GL_QUADS);
	shader_uniform_mat44(&render->shader0, "u_projection", &render->projection);
	shader_uniform_mat44(&render->shader0, "u_view", view);

	int chunk_min = -1;
	int chunk_max = 1;

	int quads = 0;
	int chunk_offset = 0;
	for (int dz = chunk_min; dz <= chunk_max; dz++)
	for (int dy = chunk_min; dy <= chunk_max; dy++)
	for (int dx = chunk_min; dx <= chunk_max; dx++) {
		struct world_chunk* chunk = world->chunks + (chunk_offset++);
		struct vec3i offset = {{dx,dy,dz}};
		quads += world_chunk_draw(chunk, &render->vtxbuf, offset);
	}
	//printf("quads: %d\n", quads);

	vtxbuf_end(&render->vtxbuf);
}

void resolve_cell_coordinates(struct vec3i* out, struct vec3* in)
{
	for (int axis = 0; axis < 3; axis++) {
		out->s[axis] = (int)floorf(in->s[axis]);
	}
}

void resolve_chunk_coordinates(struct vec3i* out, struct vec3i* in)
{
	for (int axis = 0; axis < 3; axis++) {
		out->s[axis] = ((in->s[axis] + (1<<30)) >> WORLD_CHUNK_LENGTH_EXP) - (1<<(30-WORLD_CHUNK_LENGTH_EXP));
	}
}

struct world_chunk* world_chunk_atp(struct world* world, struct vec3i* pos)
{
	int mul = 1;
	int index = 0;
	for (int axis = 0; axis < 3; axis++) {
		int v = pos->s[axis];
		if (v < -1 || v > 1) return NULL;
		index += mul * (v+1);
		mul *= 3;
	}
	return &world->chunks[index];
}


struct world_cell* world_cell_atp(struct world* world, struct vec3i* cellpos)
{
	struct vec3i chunkpos;
	resolve_chunk_coordinates(&chunkpos, cellpos);
	struct world_chunk* chunk = world_chunk_atp(world, &chunkpos);
	if (chunk == NULL) return NULL;
	struct vec3i chunkcellpos;
	for (int axis = 0; axis < 3; axis++) {
		chunkcellpos.s[axis] = cellpos->s[axis] & ((1<<WORLD_CHUNK_LENGTH_EXP)-1);
	}
	return world_chunk_cell_atp(chunk, &chunkcellpos);
}


void world_chunk_init(struct world_chunk* chunk)
{
	memset(chunk, 0, sizeof(*chunk));
	#if 1
	for (int z = 0; z < WORLD_CHUNK_LENGTH; z++)
	for (int y = 0; y < WORLD_CHUNK_LENGTH; y++)
	for (int x = 0; x < WORLD_CHUNK_LENGTH; x++) {
		struct vec3i pos = {{x,y,z}};
		struct world_cell* cell = world_chunk_cell_atp(chunk, &pos);
		cell->type = (rand()&7) == 0;
	}
	#endif
}

void world_init(struct world* world)
{
	srand(3);

	int n = 3*3*3;
	world->chunks = malloc(n * sizeof(struct world_chunk));
	for (int i = 0; i < n; i++) {
		world_chunk_init(&world->chunks[i]);
	}
}

void world_entity_clipmove(struct world* world, struct entity* entity, float dt)
{
	float move_length = vec3_length(&entity->velocity) * dt;
	int nsteps = ceilf(move_length * 32);
	if (nsteps < 1) nsteps = 1;

	float fraction = 1.0 / (float)nsteps;
	struct vec3 move;
	vec3_scale(&move, &entity->velocity, fraction);

	struct aabb selector;
	entity_get_aabb(entity, &selector);
	aabb_extend_inplace(&selector, &entity->velocity);
	aabb_grow_inplace(&selector, 0.01);

	struct vec3i cell_min, cell_max;
	resolve_cell_coordinates(&cell_min, &selector.min);
	resolve_cell_coordinates(&cell_max, &selector.max);

	int ground_collision = 0;

	for (int i = 0; i < nsteps; i++) {
		vec3_add_inplace(&entity->position, &move);

		struct vec3i_it it;
		for (vec3i_it_init(&it, &cell_min, &cell_max); it.more; vec3i_it_next(&it)) {
			struct world_cell* cell = world_cell_atp(world, &it.current);
			if (cell == NULL) continue;
			if (cell->type == 0) continue;

			struct aabb entity_aabb;
			entity_get_aabb(entity, &entity_aabb);

			struct aabb cube;
			aabb_cube_at(&cube, &it.current);

			struct aabb block;
			aabb_aabb_convolute(&block, &entity_aabb, &cube);

			int axis = 0;
			int exit_axis = 0;
			float exit_distance = 1e6;
			int exits = 0;
			for (axis = 0; axis < 3; axis++) {
				if (block.min.s[axis] > 0 || block.max.s[axis] < 0) break;
				for (int sign = 0; sign < 2; sign++) {
					struct vec3i neighbour_pos = it.current;
					neighbour_pos.s[axis] += (sign == 0) ? -1 : 1;
					struct world_cell* neighbour = world_cell_atp(world, &neighbour_pos);
					if (neighbour != NULL && neighbour->type != 0) continue;
					exits++;
					float d = (sign == 0) ? block.min.s[axis] : block.max.s[axis];
					if (fabsf(d) < fabs(exit_distance)) {
						exit_distance = d;
						exit_axis = axis;
					}
				}
			}

			if (axis != 3 || !exits) continue;
			if (fabsf(exit_distance) > 2*vec3_length(&move)) continue; // XXX isn't this a hack?

			entity->position.s[exit_axis] += exit_distance;
			entity->velocity.s[exit_axis] = 0;
			move.s[exit_axis] = 0;

			if (exit_axis == GRAVITY_AXIS && ((exit_distance * GRAVITY_SIGN) <= 0)) {
				ground_collision++;
			}
		}
	}

	entity->on_ground = ground_collision > 0;

	// planar control
	float fwd = (entity->ctrl_forward*-1) + entity->ctrl_backward;
	float sid = (entity->ctrl_left*-1) + entity->ctrl_right;
	float c = cosf(DEG2RAD(entity->yaw));
	float s = sinf(DEG2RAD(entity->yaw));
	float fdt = dt * fraction;

	if (entity->on_ground) {
		if (entity->ctrl_jump) {
			entity->on_ground = 0;
			entity->velocity.s[GRAVITY_AXIS] = -GRAVITY_SIGN * 2.0f * sqrtf(GRAVITY_FORCE * dt * JUMP_HEIGHT);
			//entity->velocity.s[GRAVITY_AXIS] = -0.15 * GRAVITY_SIGN;
			entity->ctrl_jump = 0;
		} else {
			float acceleration = 1.0 * fdt;

			entity->velocity.s[0] += acceleration * (-s * fwd + c * sid);
			entity->velocity.s[2] += acceleration * (c * fwd + s * sid);

			float friction_magnitude = 0.0002;
			float ground_friction = powf(friction_magnitude, fdt);
			for (int axis = 0; axis < 3; axis++) {
				if (axis == GRAVITY_AXIS) continue;
				entity->velocity.s[axis] *= ground_friction;
			}
		}
	} else {
		float acceleration = 1e-1 * fdt;
		entity->velocity.s[0] += acceleration * (-s * fwd + c * sid);
		entity->velocity.s[2] += acceleration * (c * fwd + s * sid);
	}
}

void entity_apply_gravity(struct entity* entity, float gravity)
{
	entity->velocity.s[GRAVITY_AXIS] += GRAVITY_SIGN * gravity;
}

int main(int argc, char** argv)
{
	SAZ(SDL_Init(SDL_INIT_VIDEO));
	atexit(SDL_Quit);

	int bitmask = SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL;
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

	struct render render;
	render_init(&render, window);

	struct world world;
	world_init(&world);

	int ctrl_forward = 0;
	int ctrl_backward = 0;
	int ctrl_left = 0;
	int ctrl_right = 0;
	int ctrl_jump = 0;

	struct entity player;
	entity_init(&player);
	{
		struct vec3 ppos = {{4.5,20,3.5}};
		player.position = ppos;
	}

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	float dt = 1.0 / 60.0;

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
				 {SDLK_SPACE, &ctrl_jump},
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
			entity_apply_gravity(&player, GRAVITY_FORCE * dt);
			player.ctrl_forward = ctrl_forward;
			player.ctrl_backward = ctrl_backward;
			player.ctrl_left = ctrl_left;
			player.ctrl_right = ctrl_right;
			player.ctrl_jump = ctrl_jump;
			world_entity_clipmove(&world, &player, dt);
			//vec3_dump(&player.velocity);
		}

		{
			float sensitivity = 0.1f;
			player.yaw += (float)mdx * sensitivity;
			player.pitch += (float)mdy * sensitivity;
			float pitch_limit = 90;
			if (player.pitch > pitch_limit) player.pitch = pitch_limit;
			if (player.pitch < -pitch_limit) player.pitch = -pitch_limit;
		}

		#if 0 // fly-code
		{
			float speed = 0.1f;
			float forward = (float)(ctrl_forward - ctrl_backward) * speed;
			float right = (float)(ctrl_right - ctrl_left) * speed;
			struct vec3 movement;
			vec3_move(&movement, yaw, pitch, forward, right);
			vec3_add_inplace(&position, &movement);
		}
		#endif

		struct mat44 view;
		mat44_set_identity(&view);
		mat44_rotate_x(&view, player.pitch);
		mat44_rotate_y(&view, player.yaw);

		struct vec3 translate;
		vec3_scale(&translate, &player.position, -1);
		mat44_translate(&view, &translate);

		if (player.on_ground) {
			glClearColor(0.6, 0.0, 0.2, 1.0);
		} else {
			glClearColor(0.0, 0.0, 0.2, 1.0);
		}
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		world_draw(&world, &render, &view);

		SDL_GL_SwapWindow(window);
	}

	SDL_DestroyWindow(window);
	SDL_GL_DeleteContext(glctx);

	return EXIT_SUCCESS;
}

