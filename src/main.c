/*
 *	Adapted smallpt to jitter the colour of any walls and 
 *  ceilings uniformly within ±0.05 of fixed RGB value (i.e. 
 *  known to 1dp) using UxHw. The intensity of the light source 
 *  is jittered by ±2% for each ray using standard Monte Carlo.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <uxhw.h>

/*
 *	Vec struct and operations
 */
typedef struct
{
	double	x;
	double	y;
	double	z;
} Vec;

static inline Vec
vecNew(double x, double y, double z)
{
	return (Vec) {x, y, z};
}

static inline Vec
vecAdd(Vec a, Vec b)
{
	return vecNew(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline Vec
vecSub(Vec a, Vec b)
{
	return vecNew(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline Vec
vecScale(Vec a, double b)
{
	return vecNew(a.x * b, a.y * b, a.z * b);
}

static inline Vec
vecMult(Vec a, Vec b)
{
	return vecNew(a.x * b.x, a.y * b.y, a.z * b.z);
}

static inline Vec
vecNorm(Vec a)
{
	double	length = 1.0 / sqrt(a.x * a.x + a.y * a.y + a.z * a.z);

	return vecScale(a, length);
}

static inline double
vecDot(Vec a, Vec b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline Vec
vecCross(Vec a, Vec b)
{
	return vecNew(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}

/*
 *	Ray
 */
typedef struct
{
	Vec	origin;
	Vec	direction;
} Ray;

/*
 *	Reflection type
 */
typedef enum
{
	kReflDiff,
	kReflSpec,
	kReflRefr
} ReflType;

/*
 *	Sphere
 */
typedef struct
{
	double		radius;
	Vec		position;
	Vec		emission;
	Vec		colour;
	ReflType	refl;
} Sphere;

/*
 *	Sphere intersection
 */
double
sphereIntersect(const Sphere *  s, const Ray *  r)
{
	Vec		op = vecSub(s->position, r->origin);
	double		b = vecDot(op, r->direction);
	double		det = b * b - vecDot(op, op) + s->radius * s->radius;
	double		t;
	double		eps = 1e-4;

	if (det < 0)
	{
		return 0;
	}

	det = sqrt(det);
	t = b - det;
	if (t > eps)
	{
		return t;
	}
	t = b + det;

	return (t > eps) ? t : 0;
}

/*
 *	Scene (light source is the last sphere)
 */
Sphere spheres[] =
{
	{1e5, { 1e5 + 1, 40.8, 81.6}, {0, 0, 0}, {.75, .25, .25}, kReflDiff},	/* Left */
	{1e5, {-1e5 + 99, 40.8, 81.6}, {0, 0, 0}, {.25, .25, .75}, kReflDiff},	/* Right */
	{1e5, {50, 40.8,  1e5},        {0, 0, 0}, {.75, .75, .75}, kReflDiff},	/* Back */
	{1e5, {50, 40.8, -1e5 + 170},  {0, 0, 0}, {0, 0, 0},      kReflDiff},	/* Front */
	{1e5, {50,  1e5, 81.6},        {0, 0, 0}, {.75, .75, .75}, kReflDiff},	/* Bottom */
	{1e5, {50, -1e5 + 81.6, 81.6}, {0, 0, 0}, {.75, .75, .75}, kReflDiff},	/* Top */
	{16.5,{27, 16.5, 47},          {0, 0, 0}, {0.999, 0.999, 0.999}, kReflSpec},	/* Mirror */
	{16.5,{73, 16.5, 78},          {0, 0, 0}, {0.999, 0.999, 0.999}, kReflRefr},	/* Glass */
	{600, {50, 681.6 - .27, 81.6}, {12, 12, 12}, {0, 0, 0},   kReflDiff}	/* Light */
};

/*
 *	Utility
 */
static inline double
clamp(double x)
{
	return x < 0 ? 0 : (x > 1 ? 1 : x);
}

static inline int
toInt(double x)
{
	return (int)(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

/*
 *	Intersect all spheres
 */
bool
intersect(const Ray *  r, double *  t, int *  id)
{
	double		n = sizeof(spheres) / sizeof(Sphere);
	double		d;
	double		inf = *t = 1e20;
	int		ii;

	*id = -1;

	for (ii = 0; ii < (int)n; ii++)
	{
		if ((d = sphereIntersect(&spheres[ii], r)) && d < *t)
		{
			*t = d;
			*id = ii;
		}
	}

	return *id != -1;
}

/*
 *	Generate random number in [0,1)
 */
static inline double
rnd(unsigned short *  Xi)
{
	return erand48(Xi);
}

/*
 *	Apply colour jitter using UxHw uniform within ±0.05
 */
static double jitterX = UxHwDoubleUniformDist(-0.05, 0.05);
static double jitterY = UxHwDoubleUniformDist(-0.05, 0.05);
static double jitterZ = UxHwDoubleUniformDist(-0.05, 0.05);

static inline Vec
jitterColour(Vec c)
{
	return vecNew(
		c.x + jitterX,
		c.y + jitterY,
		c.z + jitterZ);
}

/*
 *	Apply small realistic light jitter ±0.02 (2% fluctuations)
 */
static inline Vec
jitterLight(Vec c, unsigned short *  Xi)
{
	return vecScale(c, 1.0 + (rnd(Xi) * 0.04 - 0.02));
}

/*
 *  Radiance
 */
Vec
radiance(Ray r, int depth, unsigned short *  Xi)
{
	double		t;
	int		id = 0;

	if (!intersect(&r, &t, &id))
	{
		return vecNew(0, 0, 0);
	}

	const Sphere *	obj = &spheres[id];
	Vec		x = vecAdd(r.origin, vecScale(r.direction, t));
	Vec		n = vecNorm(vecSub(x, obj->position));
	Vec		nl = vecDot(n, r.direction) < 0 ? n : vecScale(n, -1);

	/*
	 *	Surface colour:
	 *	- Diffuse: apply jitter
	 *	- Specular/Refractive/Light: no jitter
	 */
	Vec		f = obj->colour;
	if (obj->refl == kReflDiff)
	{
		f = jitterColour(f);
	}

	/*
	 *	Emission:
	 *	- Apply jitter only for designated light source
	 */
	Vec		emission = obj->emission;
	if (id == (sizeof(spheres) / sizeof(Sphere)) - 1)
	{
		emission = jitterLight(emission, Xi);
	}

	/*
	 *	Russian roulette
	 */
	double		p = f.x > f.y && f.x > f.z ? f.x : (f.y > f.z ? f.y : f.z);

	if (++depth > 5)
	{
		if (rnd(Xi) < p)
		{
			f = vecScale(f, 1 / p);
		}
		else
		{
			return emission;
		}
	}

	/*
	 *	Reflection models
	 */
	if (obj->refl == kReflDiff)
	{
		double		r1 = 2 * M_PI * rnd(Xi);
		double		r2 = rnd(Xi);
		double		r2s = sqrt(r2);
		Vec		w = nl;
		Vec		u = vecNorm(vecCross(fabs(w.x) > .1 ? vecNew(0, 1, 0) : vecNew(1, 0, 0), w));
		Vec		v = vecCross(w, u);
		Vec		d = vecNorm(vecAdd(
					vecAdd(vecScale(u, cos(r1) * r2s),
						vecScale(v, sin(r1) * r2s)),
					vecScale(w, sqrt(1 - r2))));

		return vecAdd(emission, vecMult(f, radiance((Ray) {x, d}, depth, Xi)));
	}
	else if (obj->refl == kReflSpec)
	{
		Vec		reflDir = vecSub(r.direction, vecScale(n, 2 * vecDot(n, r.direction)));

		return vecAdd(emission, vecMult(f, radiance((Ray) {x, reflDir}, depth, Xi)));
	}

	/*
	 *	Refraction
	 */
	Ray		reflRay = {x, vecSub(r.direction, vecScale(n, 2 * vecDot(n, r.direction)))};
	bool		into = vecDot(n, nl) > 0;
	double		nc = 1;
	double		nt = 1.5;
	double		nnt = into ? nc / nt : nt / nc;
	double		ddn = vecDot(r.direction, nl);
	double		cos2t;

	if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
	{
		return vecAdd(emission, vecMult(f, radiance(reflRay, depth, Xi)));
	}

	Vec		tdir = vecNorm(vecSub(
				vecScale(r.direction, nnt),
				vecScale(n, (into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))));
	double		a = nt - nc;
	double		b = nt + nc;
	double		R0 = a * a / (b * b);
	double		c = 1 - (into ? -ddn : vecDot(tdir, n));
	double		Re = R0 + (1 - R0) * pow(c, 5);
	double		Tr = 1 - Re;
	double		P = .25 + .5 * Re;
	double		RP = Re / P;
	double		TP = Tr / (1 - P);

	if (depth > 2)
	{
		return vecAdd(emission, vecMult(f,
			rnd(Xi) < P
				? vecScale(radiance(reflRay, depth, Xi), RP)
				: vecScale(radiance((Ray) {x, tdir}, depth, Xi), TP)));
	}
	else
	{
		return vecAdd(emission, vecMult(f,
			vecAdd(vecScale(radiance(reflRay, depth, Xi), Re),
				vecScale(radiance((Ray) {x, tdir}, depth, Xi), Tr))));
	}
}


/*
 *	Main
 */
int
main(int argc, char *   argv[])
{
	int		width = 1024;
	int		height = 768;
	int		samples = (argc == 2) ? atoi(argv[1]) / 4 : 1;
	Ray		cam = {vecNew(50, 52, 295.6), vecNorm(vecNew(0, -0.042612, -1))};
	Vec		cx = vecNew(width * .5135 / height, 0, 0);
	Vec		cy = vecScale(vecNorm(vecCross(cx, cam.direction)), .5135);
	Vec *		c = calloc(width * height, sizeof(Vec));
	int		yy;
	int		xx;

	#pragma omp parallel for schedule(dynamic,1)
	for (yy = 0; yy < height; yy++)
	{
		fprintf(stderr, "\rRendering %d spp %5.2f%%", samples * 4, 100. * yy / (height - 1));
		unsigned short	Xi[3] = {0, 0, (unsigned short)(yy * yy * yy)};

		for (xx = 0; xx < width; xx++)
		{
			int		ii = (height - yy - 1) * width + xx;
			Vec		r = vecNew(0, 0, 0);
			int		sy;
			int		sx;

			for (sy = 0; sy < 2; sy++)
			{
				for (sx = 0; sx < 2; sx++)
				{
					Vec		sub = vecNew(0, 0, 0);
					int		ss;

					for (ss = 0; ss < samples; ss++)
					{
						double		r1 = 2 * rnd(Xi);
						double		dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						double		r2 = 2 * rnd(Xi);
						double		dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
						Vec		d = vecAdd(
								vecAdd(
									vecScale(cx, ((sx + .5 + dx) / 2 + xx) / width - .5),
									vecScale(cy, ((sy + .5 + dy) / 2 + yy) / height - .5)),
								cam.direction);
						Ray		ray = {vecAdd(cam.origin, vecScale(d, 140)), vecNorm(d)};
						Vec		rad = radiance(ray, 0, Xi);

						/*
						 *	Apply small light jitter only if hitting light
						 */
						if (ii == width * height - 1)
						{
							rad = jitterLight(rad, Xi);
						}

						sub = vecAdd(sub, vecScale(rad, 1.0 / samples));
					}

					r = vecAdd(r, vecScale(vecNew(clamp(sub.x), clamp(sub.y), clamp(sub.z)), .25));
				}
			}

			c[ii] = r;
		}
	}

	FILE *	f = fopen("mountDir/image.ppm", "w");

	fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);

	int		ii;

	for (ii = 0; ii < width * height; ii++)
	{
		fprintf(f, "%d %d %d ",
			toInt(c[ii].x),
			toInt(c[ii].y),
			toInt(c[ii].z));
	}

	free(c);

	return 0;
}