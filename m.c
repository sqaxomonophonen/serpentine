#include <stdio.h>
#include <math.h>
#include "m.h"

float calc_bezier(float t, float a, float b, float c, float d)
{
	float t2 = t*t;
	float t3 = t*t*t;
	float ti = (1-t);
	float t2i = ti*ti;
	float t3i = ti*ti*ti;
	return a*t3i + 3*b*t2i*t + 3*c*ti*t2 + d*t3;
}

float calc_bezier_deriv(float t, float a, float b, float c, float d)
{
	a = 3*(b-a);
	b = 3*(c-b);
	c = 3*(d-c);
	float t2 = t*t;
	float ti = (1-t);
	float t2i = ti*ti;
	return a*t2i + 2*b*ti*t + 3*c*t2;
}

void vec3_dump(struct vec3* x)
{
	printf("(%.3f  %.3f  %.3f)\n", x->s[0], x->s[1], x->s[2]);
}

void vec3_zero(struct vec3* x)
{
	for (int i = 0; i < 3; i++) x->s[i] = 0;
}

void vec3_copy(struct vec3* dst, struct vec3* src)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = src->s[i];
	}
}

void vec3_scale(struct vec3* dst, struct vec3* src, float scalar)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = src->s[i] * scalar;
	}
}

void vec3_add_inplace(struct vec3* dst, struct vec3* src)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] += src->s[i];
	}
}

void vec3_add_scaled_inplace(struct vec3* dst, struct vec3* src, float scalar)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] += src->s[i] * scalar;
	}
}

void vec3_add(struct vec3* dst, struct vec3* a, struct vec3* b)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = a->s[i] + b->s[i];
	}
}

void vec3_sub(struct vec3* dst, struct vec3* a, struct vec3* b)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = a->s[i] - b->s[i];
	}
}

void vec3_scale_inplace(struct vec3* dst, float scalar)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] *= scalar;
	}
}

void vec3_lerp(struct vec3* dst, struct vec3* a, struct vec3* b, float t)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = a->s[i] + (b->s[i] - a->s[i]) * t;
	}
}

float vec3_dot(struct vec3* a, struct vec3* b)
{
	float r = 0;
	for (int i = 0; i < 3; i++) {
		r += a->s[i] * b->s[i];
	}
	return r;
}

float vec3_length(struct vec3* v)
{
	return sqrtf(vec3_dot(v, v));
}


void vec3_cross(struct vec3* dst, struct vec3* a, struct vec3* b)
{
	for (int i = 0; i < 3; i++) {
		int i1 = (i+1)%3;
		int i2 = (i+2)%3;
		dst->s[i] = a->s[i1]*b->s[i2] - a->s[i2]*b->s[i1];
	}
}

void vec3_normalize_inplace(struct vec3* dst)
{
	vec3_scale_inplace(dst, 1 / vec3_length(dst));
}

void vec3_move(struct vec3* move, float yaw, float pitch, float forward, float right)
{
	vec3_zero(move);

	float yaw_s = sinf(DEG2RAD(yaw));
	float yaw_c = cosf(DEG2RAD(yaw));

	float pitch_s = sinf(DEG2RAD(pitch));
	float pitch_c = cosf(DEG2RAD(pitch));

	struct vec3 fv = {{
		yaw_s * pitch_c,
		-pitch_s,
		-yaw_c * pitch_c
	}};
	vec3_add_scaled_inplace(move, &fv, forward);

	struct vec3 rv = {{
		yaw_c,
		0,
		yaw_s,
	}};
	vec3_add_scaled_inplace(move, &rv, right);
}

void vec3_z_up_move(struct vec3* move, float yaw, float pitch, float forward, float right)
{
	vec3_zero(move);

	float yaw_s = sinf(DEG2RAD(yaw));
	float yaw_c = cosf(DEG2RAD(yaw));

	float pitch_s = sinf(DEG2RAD(pitch));
	float pitch_c = cosf(DEG2RAD(pitch));

	struct vec3 fv = {{
		yaw_s * pitch_c,
		-yaw_c * pitch_c,
		-pitch_s,
	}};
	vec3_add_scaled_inplace(move, &fv, forward);

	struct vec3 rv = {{
		yaw_c,
		yaw_s,
		0,
	}};
	vec3_add_scaled_inplace(move, &rv, right);
}


void vec3_bezier(struct vec3* dst, float t, struct vec3* a, struct vec3* b, struct vec3* c, struct vec3* d)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = calc_bezier(t, a->s[i], b->s[i], c->s[i], d->s[i]);
	}
}

void vec3_bezier_deriv(struct vec3* dst, float t, struct vec3* a, struct vec3* b, struct vec3* c, struct vec3* d)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = calc_bezier_deriv(t, a->s[i], b->s[i], c->s[i], d->s[i]);
	}
}

void vec3_complete_basis(struct vec3* normal, struct vec3* a, struct vec3* b)
{
	struct vec3 q;
	if (normal->s[0] > normal->s[1] && normal->s[0] > normal->s[2]) {
		struct vec3 q0 = {{0,1,0}};
		vec3_copy(&q, &q0);
	} else {
		struct vec3 q1 = {{1,0,0}};
		vec3_copy(&q, &q1);
	}

	vec3_cross(a, normal, &q);
	vec3_cross(b, normal, a);
	vec3_normalize_inplace(a);
	vec3_normalize_inplace(b);
}

void vec3_calculate_normal_from_3_points(struct vec3* normal, struct vec3* points3)
{
	struct vec3 a,b;
	vec3_sub(&a, &points3[1], &points3[0]);
	vec3_sub(&b, &points3[2], &points3[0]);
	vec3_cross(normal, &a, &b);
	vec3_normalize_inplace(normal);
}

void vec4_dump(struct vec4* x)
{
	printf("(%.4f  %.4f  %.4f  %.4f)\n", x->s[0], x->s[1], x->s[2], x->s[3]);
}

void vec4_copy(struct vec4* dst, struct vec4* src)
{
	for (int i = 0; i < 4; i++) dst->s[i] = src->s[i];
}

float vec4_dot(struct vec4* a, struct vec4* b)
{
	float r = 0;
	for (int i = 0; i < 4; i++) {
		r += a->s[i] * b->s[i];
	}
	return r;
}

static inline void mat44_get_row(struct mat44* m, struct vec4* result, int row)
{
	for (int col = 0; col < 4; col++) {
		result->s[col] = mat44_at(m, col, row);
	}
}

static inline void mat44_get_col(struct mat44* m, struct vec4* result, int col)
{
	for (int row = 0; row < 4; row++) {
		result->s[row] = mat44_at(m, col, row);
	}
}

static void mat44_set_zero(struct mat44* m)
{
	for (int i = 0; i < 16; i++) {
		m->s[i] = 0;
	}
}

void mat44_dump(struct mat44* x)
{
	for (int row = 0; row < 4; row++) {
		struct vec4 v;
		mat44_get_row(x, &v, row);
		vec4_dump(&v);
	}
}

void mat44_copy(struct mat44* dst, struct mat44* src)
{
	for (int i = 0; i < 16; i++) {
		dst->s[i] = src->s[i];
	}
}

void mat44_multiply(struct mat44* dst, struct mat44* a, struct mat44* b)
{
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			struct vec4 arow;
			mat44_get_row(a, &arow, row);
			struct vec4 bcol;
			mat44_get_col(b, &bcol, col);
			*(mat44_atp(dst,col,row)) = vec4_dot(&arow, &bcol);
		}
	}
}

void mat44_multiply_inplace(struct mat44* dst, struct mat44* b)
{
	struct mat44 a;
	mat44_copy(&a, dst);
	mat44_multiply(dst, &a, b);
}

void mat44_set_identity(struct mat44* m)
{
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			*(mat44_atp(m,col,row)) = col == row ? 1 : 0;
		}
	}
}

static void mat44_set_frustum(struct mat44* m, float left, float right, float bottom, float top, float znear, float zfar)
{
	// http://docs.gl/gl2/glFrustum

	mat44_set_zero(m);

	*(mat44_atp(m,0,0)) = 2 * znear / (right - left);
	*(mat44_atp(m,1,1)) = 2 * znear / (top - bottom);

	*(mat44_atp(m,2,0)) = (right+left)/(right-left); // A
	*(mat44_atp(m,2,1)) = (top+bottom)/(top-bottom); // B
	*(mat44_atp(m,2,2)) = -((zfar+znear)/(zfar-znear)); // C
	*(mat44_atp(m,2,3)) = -1;

	*(mat44_atp(m,3,2)) = -(2 * zfar * znear / (zfar - znear)); // D
}

void mat44_set_perspective(struct mat44* m, float fovy, float aspect, float znear, float zfar)
{
	float dy = znear * tanf(DEG2RAD(fovy)/2);
	float dx = dy * aspect;
	mat44_set_frustum(m, -dx, dx, -dy, dy, znear, zfar);
}

float mat44_get_znear(struct mat44* m)
{
	// http://shitohichiumaya.blogspot.com/2012/10/how-to-get-camera-parameters-from.html
	float c = mat44_at(m,2,2);
	float k = (c - 1.0f) / (c + 1.0f);
	float d = mat44_at(m,3,2);
	return (d * (1.0f - k)) / (2.0f * k);
}

void mat44_set_rotation(struct mat44* m, float angle, struct vec3* axis)
{
	mat44_set_identity(m);

	float c = cosf(DEG2RAD(angle));
	float s = sinf(DEG2RAD(angle));
	float c1 = 1 - c;

	for (int row = 0; row < 3; row++) {
		for (int col = 0; col < 3; col++) {
			float value = axis->s[col] * axis->s[row] * c1;
			if (col == row) {
				value += c;
			} else {
				int i = 3 - row - col;
				ASSERT(i >= 0 && i < 3);
				float sgn1 = ((col+row)&1) ? 1 : -1;
				float sgn2 = row>col ? 1 : -1;
				value += sgn1 * sgn2 * axis->s[i] * s;

			}
			*(mat44_atp(m,col,row)) = value;
		}
	}
}

void mat44_rotate(struct mat44* m, float angle, struct vec3* axis)
{
	struct mat44 op;
	mat44_set_rotation(&op, angle, axis);
	mat44_multiply_inplace(m, &op);
}

void mat44_rotate_x(struct mat44* m, float angle)
{
	struct vec3 axis = {{1,0,0}};
	mat44_rotate(m, angle, &axis);
}

void mat44_rotate_y(struct mat44* m, float angle)
{
	struct vec3 axis = {{0,1,0}};
	mat44_rotate(m, angle, &axis);
}

void mat44_rotate_z(struct mat44* m, float angle)
{
	struct vec3 axis = {{0,0,1}};
	mat44_rotate(m, angle, &axis);
}

void mat44_set_translation(struct mat44* m, struct vec3* delta)
{
	mat44_set_identity(m);
	for (int row = 0; row < 3; row++) {
		*(mat44_atp(m,3,row)) = delta->s[row];
	}
}

void mat44_translate(struct mat44* m, struct vec3* delta)
{
	struct mat44 op;
	mat44_set_translation(&op, delta);
	mat44_multiply_inplace(m, &op);
}

static float _mat44_sub33_at(struct mat44* m, int delcol, int delrow, int col, int row)
{
	return mat44_at(m, col + (col >= delcol ? 1 : 0), row + (row >= delrow ? 1 : 0));
}

static void mat44_minors(struct mat44* dst, struct mat44* src)
{
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			float det = 0;
			for (int x = 0; x < 3; x++) {
				{
					float s = 1;
					for (int i = 0; i < 3; i++) {
						s *= _mat44_sub33_at(src, col, row, (i+x+6)%3, i);
					}
					det += s;
				}
				{
					float t = 1;
					for (int i = 0; i < 3; i++) {
						t *= _mat44_sub33_at(src, col, row, (-i+x+6)%3, i);
					}
					det -= t;
				}
			}
			*(mat44_atp(dst, col, row)) = det;
		}
	}
}

static void mat44_cofactors(struct mat44* dst, struct mat44* src)
{
	mat44_minors(dst, src);
	for (int col = 0; col < 4; col++) {
		for (int row = 0; row < 4; row++) {
			if ((col+row)&1) {
				*(mat44_atp(dst, col, row)) *= -1;
			}
		}
	}
}

static void mat44_transpose(struct mat44* dst, struct mat44* src)
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			*(mat44_atp(dst, j, i)) = mat44_at(src, i, j);
		}
	}
}

static void mat44_adjugate(struct mat44* dst, struct mat44* src)
{
	struct mat44 cofactor;
	mat44_cofactors(&cofactor, src);
	mat44_transpose(dst, &cofactor);
}

static float mat44_determinant(struct mat44* a)
{
	float det = 0;
	for (int x = 0; x < 4; x++) {
		{
			float s = 1;
			for (int i = 0; i < 4; i++) {
				s *= mat44_at(a, (i+x)&3, i);
			}
			det += s;
		}
		{
			float t = 1;
			for (int i = 0; i < 4; i++) {
				t *= mat44_at(a, (-i+x)&3, i);
			}
			det -= t;
		}
	}
	return det;
}

void mat44_inverse(struct mat44* dst, struct mat44* src)
{
	mat44_adjugate(dst, src);
	float det = mat44_determinant(src);
	AN(det);
	float scalar = 1 / det;
	for (int i = 0; i < 16; i++) {
		dst->s[i] *= scalar;
	}
}

void mat44_get_bases(struct mat44* m, struct vec3* x, struct vec3* y, struct vec3* z)
{
	struct vec3* vs[] = {x,y,z};
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (vs[j] != NULL) vs[j]->s[i] = mat44_at(m, i, j);
		}
	}
}

void mat44_set_z_up_identity(struct mat44* m)
{
	mat44_set_zero(m);
	*(mat44_atp(m,0,0)) = 1;
	*(mat44_atp(m,1,2)) = 1;
	*(mat44_atp(m,2,1)) = 1;
	*(mat44_atp(m,3,3)) = 1;
}

void vec3_from_vec4(struct vec3* dst, struct vec4* src)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = src->s[i] / src->s[3];
	}
}

void vec4_from_vec3(struct vec4* dst, struct vec3* src)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = src->s[i];
	}
	dst->s[3] = 1;
}

void vec3_apply_mat44(struct vec3* dst, struct vec3* src, struct mat44* m)
{
	struct vec4 result;
	struct vec4 src4;
	vec4_from_vec3(&src4, src);
	for (int i = 0; i < 4; i++) {
		struct vec4 o;
		mat44_get_row(m, &o, i);
		result.s[i] = vec4_dot(&src4, &o);
	}
	vec3_from_vec4(dst, &result);
}

void vec3_apply_rotation_mat44(struct vec3* dst, struct vec3* src, struct mat44* m)
{
	for (int i = 0; i < 3; i++) {
		struct vec4 o4;
		mat44_get_row(m, &o4, i);
		o4.s[3] = 1;
		struct vec3 o3;
		vec3_from_vec4(&o3, &o4);
		dst->s[i] = vec3_dot(src, &o3);
	}
}

void vec4_apply_mat44_to_vec3(struct vec4* dst, struct vec3* src, struct mat44* m)
{
	struct vec4 src4;
	vec4_from_vec3(&src4, src);
	for (int i = 0; i < 4; i++) {
		struct vec4 o;
		mat44_get_row(m, &o, i);
		dst->s[i] = vec4_dot(&src4, &o);
	}
}

