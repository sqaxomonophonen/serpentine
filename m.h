#ifndef M_H
#define M_H

#include "a.h"

#ifndef M_PI
#define M_PI (3.141592653589793)
#endif

#define DEG2RAD(x) ((x)/180.0f*M_PI)
#define I2RAD(x) ((x)*M_PI*2.0f)

// XXX dangerous due to potential double calculation
//#define MAX(a,b) ((a)>(b) ? (a) : (b))
//#define MIN(a,b) ((a)<(b) ? (a) : (b))

float calc_bezier(float t, float a, float b, float c, float d);
float calc_bezier_deriv(float t, float a, float b, float c, float d);

struct vec3 {
	float s[3];
};

void vec3_dump(struct vec3* x);
void vec3_zero(struct vec3* x);
void vec3_copy(struct vec3* dst, struct vec3* src);
void vec3_scale(struct vec3* dst, struct vec3* src, float scalar);
void vec3_add_inplace(struct vec3* dst, struct vec3* src);
void vec3_add_scaled_inplace(struct vec3* dst, struct vec3* src, float scalar);
void vec3_add(struct vec3* dst, struct vec3* a, struct vec3* b);
void vec3_sub(struct vec3* dst, struct vec3* a, struct vec3* b);
void vec3_scale_inplace(struct vec3* dst, float scalar);
void vec3_lerp(struct vec3* dst, struct vec3* a, struct vec3* b, float t);
float vec3_dot(struct vec3* a, struct vec3* b);
float vec3_length(struct vec3* v);
void vec3_cross(struct vec3* dst, struct vec3* a, struct vec3* b);
void vec3_normalize_inplace(struct vec3* dst);
void vec3_move(struct vec3* move, float yaw, float pitch, float forward, float right);
void vec3_z_up_move(struct vec3* move, float yaw, float pitch, float forward, float right);

void vec3_bezier(struct vec3* dst, float t, struct vec3* a, struct vec3* b, struct vec3* c, struct vec3* d);
void vec3_bezier_deriv(struct vec3* dst, float t, struct vec3* a, struct vec3* b, struct vec3* c, struct vec3* d);

void vec3_complete_basis(struct vec3* normal, struct vec3* a, struct vec3* b);
void vec3_calculate_normal_from_3_points(struct vec3* normal, struct vec3* points3);

struct vec4 {
	float s[4];
};
void vec4_dump(struct vec4* x);
void vec4_copy(struct vec4* dst, struct vec4* src);
float vec4_dot(struct vec4* a, struct vec4* b);


// column-major order
struct mat44 {
	float s[16];
};

static inline int mat44_ati(int col, int row)
{
	ASSERT(col >= 0 && col < 4 && row >= 0 && row < 4);
	return row + col * 4;
}

static inline float* mat44_atp(struct mat44* m, int col, int row)
{
	return &m->s[mat44_ati(col, row)];
}

static inline float mat44_at(struct mat44* m, int col, int row)
{
	return *(mat44_atp(m, col, row));
}

void mat44_dump(struct mat44* x);
void mat44_copy(struct mat44* dst, struct mat44* src);
void mat44_multiply(struct mat44* dst, struct mat44* a, struct mat44* b);
void mat44_multiply_inplace(struct mat44* dst, struct mat44* b);
void mat44_set_identity(struct mat44* m);
void mat44_set_perspective(struct mat44* m, float fovy, float aspect, float znear, float zfar);
float mat44_get_znear(struct mat44* m);
void mat44_set_rotation(struct mat44* m, float angle, struct vec3* axis);
void mat44_rotate(struct mat44* m, float angle, struct vec3* axis);
void mat44_rotate_x(struct mat44* m, float angle);
void mat44_rotate_y(struct mat44* m, float angle);
void mat44_rotate_z(struct mat44* m, float angle);
void mat44_set_translation(struct mat44* m, struct vec3* delta);
void mat44_translate(struct mat44* m, struct vec3* delta);
void mat44_inverse(struct mat44* dst, struct mat44* src);
void mat44_get_bases(struct mat44* m, struct vec3* x, struct vec3* y, struct vec3* z);
void mat44_set_z_up_identity(struct mat44* m);

// combined
void vec3_from_vec4(struct vec3* dst, struct vec4* src);
void vec4_from_vec3(struct vec4* dst, struct vec3* src);
void vec3_apply_mat44(struct vec3* dst, struct vec3* src, struct mat44* m);
void vec3_apply_rotation_mat44(struct vec3* dst, struct vec3* src, struct mat44* m);
void vec4_apply_mat44_to_vec3(struct vec4* dst, struct vec3* src, struct mat44* m);

#endif/*M_H*/
