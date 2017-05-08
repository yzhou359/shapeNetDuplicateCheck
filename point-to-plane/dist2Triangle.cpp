#include "dist2Triangle.h"

retValues dist2Triangle(const double* p0,
	const double* p1,
	const double* p2,
	const double* p) {
	retValues r;
	double vb[3];
	vec3copy(p0, vb);

	double ve0[3], ve1[3], vd[3];

	vec3sub(p1, vb, ve0);
	vec3sub(p2, vb, ve1);
	vec3sub(vb, p, vd);

	double a, b, c, d, e, f, distance;

	a = vec3dot(ve0, ve0);
	b = vec3dot(ve0, ve1);
	c = vec3dot(ve1, ve1);
	d = vec3dot(ve0, vd);
	e = vec3dot(ve1, vd);
	f = vec3dot(vd, vd);

	double  det = (a * c - b * b);
	double  s(b * e - c * d), t(b * d - a * e);

	if (s + t <= det)
		if (s < (double)0.0)
			if (t < (double)0.0) // region 4
				if (d < (double)0.0) {
					t = (double)0.0;
					if (-d >= a) {
						s = (double)1.0;
						distance = a + 2.0 * d + f;
					}
					else {
						s = -d / a;
						distance = d * s + f;
					}
				}
				else {
					s = (double)0.0;
					if (e >= (double)0.0) {
						t = (double)0.0;
						distance = f;
					}
					else if (-e >= c) {
						t = (double)1.0;
						distance = c + 2.0 * e + f;
					}
					else {
						t = -e / c;
						distance = e * t + f;
					}
				}
			else { // region 3
				s = (double)0.0;
				if (e >= (double)0.0) {
					t = (double)0.0;
					distance = f;
				}
				else if (-e >= c) {
					t = (double)1.0;
					distance = c + 2.0 * e + f;
				}
				else {
					t = -e / c;
					distance = e * t + f;
				}
			}
		else if (t < (double)0.0) { // region 5
			t = (double)0.0;
			if (d >= (double)0.0) {
				s = (double)0.0;
				distance = f;
			}
			else if (-d >= a) {
				s = (double)1.0;
				distance = a + 2.0 * d + f;
			}
			else {
				s = -d / a;
				distance = d * s + f;
			}
		}
		else { // region 0

			const double invDet((double)1.0 / det);

			s *= invDet, t *= invDet;
			distance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
		}
	else {
		double tmp0, tmp1, numer, denom;

		if (s < (double)0.0) // region 2
			if ((tmp1 = c + e) >(tmp0 = b + d))
				if ((numer = tmp1 - tmp0) >= (denom = a - 2.0 * b + c)) {
					s = (double)1.0;
					t = (double)0.0;
					distance = a + 2.0 * d + f;
				}
				else {
					s = numer / denom, t = (double)1.0 - s;
					distance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
				}
			else {
				s = (double)0.0;
				if (tmp1 <= (double)0.0) {
					t = (double)1.0;
					distance = c + 2.0 * e + f;
				}
				else if (e >= (double)0.0) {
					t = (double)0.0;
					distance = f;
				}
				else {
					t = -e / c;
					distance = e * t + f;
				}
			}
		else if (t < (double)0.0) { // region 6 
			if ((tmp1 = a + d) >(tmp0 = b + e))
				if ((numer = tmp1 - tmp0) >= (denom = a - 2.0 * b + c)) {
					t = (double)1.0;
					s = (double)0.0;
					distance = c + 2.0 * e + f;
				}
				else {
					t = numer / denom, s = (double)1.0 - t;
					distance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
				}
			else {
				t = 0.0;
				if (tmp1 <= (double)0.0) {
					s = (double)1.0;
					distance = a + 2.0 * d + f;
				}
				else if (d >= (double)0.0) {
					s = (double)0.0;
					distance = f;
				}
				else {
					s = -d / a;
					distance = d * s + f;;
				}
			}
		}
		else { // region 1
			if ((numer = c + e - b - d) <= (double)0.0) {
				s = (double)0.0;
				t = (double)1.0;
				distance = c + 2.0 * e + f;
			}
			else {
				if (numer >= (denom = a - 2.0 * b + c)) {
					s = (double)1.0;
					t = (double)0.0;
					distance = a + 2.0 * d + f;
				}
				else {
					s = numer / denom, t = (double)1.0 - s;
					distance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
				}
			}
		}
	}

	r.distance = distance;
	r.s = s;
	r.t = t;
	return r;
}