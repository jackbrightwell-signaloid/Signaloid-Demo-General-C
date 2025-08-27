/*
 *	Adapted smallpt to jitter the colour of walls and ceilings
 *	uniformly within ±0.05 using UxHw, with ±2% light jitter
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <uxhw.h>
#include <string.h>

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
		a.x * b.y - a.y * b.x
	);
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
sphereIntersect(
	const Sphere *	s,
	const Ray *	r
)
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
 *	Scene (light source is last sphere)
 */
Sphere	spheres[] =
{
	{1e5, { 1e5 + 1, 40.8, 81.6}, {0, 0, 0}, {.75, .25, .25}, kReflDiff},	/* Left */
	{1e5, {-1e5 + 99, 40.8, 81.6}, {0, 0, 0}, {.25, .25, .75}, kReflDiff},	/* Right */
	{1e5, {50, 40.8,  1e5},        {0, 0, 0}, {.75, .75, .75}, kReflDiff},	/* Back */
	{1e5, {50, 40.8, -1e5 + 170},  {0, 0, 0}, {0, 0, 0},       kReflDiff},	/* Front */
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
 *	Find first sphere intersection (if any)
 */
bool
intersect(
	const Ray *	r,
	double *	t,
	int *	id
)
{
	double		n = sizeof(spheres) / sizeof(Sphere);
	double		d;
	int		ii;

	*t = 1e20;
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
 *	Random number generator
 */
static inline double
rnd(unsigned short *	Xi)
{
	return erand48(Xi);
}

/*
 *	Initialise RGB jitter (assuming independent jitter 
 *	distributions for each surface, same distributions 
 * 	for each surface)
 */
static double jitterR[sizeof(spheres)/sizeof(Sphere)];
static double jitterG[sizeof(spheres)/sizeof(Sphere)];
static double jitterB[sizeof(spheres)/sizeof(Sphere)];

/*
 *	Initialize jitter
 */
static inline void
initJitter()
{
	int	ii;
	int	n = sizeof(spheres) / sizeof(Sphere);

	for (ii = 0; ii < n; ii++)
	{
		if (spheres[ii].refl == kReflDiff)
		{
			jitterR[ii] = UxHwDoubleUniformDist(-0.05, 0.05);
			jitterG[ii] = UxHwDoubleUniformDist(-0.05, 0.05);
			jitterB[ii] = UxHwDoubleUniformDist(-0.05, 0.05);
		}
		else
		{
			jitterR[ii] = jitterG[ii] = jitterB[ii] = 0.0;
		}
	}
}

/*
 *	Define jitter function (use same distribution within surface)
 */
static inline Vec
jitterColour(
	Vec	c,
	int	id
)
{
	return vecNew(c.x + jitterR[id],
		      c.y + jitterG[id],
		      c.z + jitterB[id]);
}

/*
 *	Apply small realistic light jitter (2% fluctuations)
 */
static inline Vec
jitterLight(
	Vec	c,
	unsigned short *	Xi
)
{
	return vecScale(c, 1.0 + (rnd(Xi) * 0.04 - 0.02));
}

/*
 *	Radiance - core path tracing function, returns RGB radiance 
 *	along ray recursively.
 */
Vec
radiance(
	Ray	r,
	int	depth,
	unsigned short *	Xi
)
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
		f = jitterColour(f, id);
	}

/*
	 *	Emission:
	 *	- Apply jitter only for designated light source
	 */
	Vec		emission = obj->emission;
	if (id == (sizeof(spheres)/sizeof(Sphere)) - 1)
	{
		emission = jitterLight(emission, Xi);
	}

	/*
	 *	Russian roulette - randomly terminate after depth 5
	 */

	double		px = UxHwDoubleNthMoment(f.x, 1);
	double		py = UxHwDoubleNthMoment(f.y, 1);
	double		pz = UxHwDoubleNthMoment(f.z, 1);
	double		p = px > py && px > pz ? px : (py > pz ? py : pz);

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

	/* Diffuse */
	if (obj->refl == kReflDiff)
	{
		double		r1 = 2 * M_PI * rnd(Xi);
		double		r2 = rnd(Xi);
		double		r2s = sqrt(r2);
		Vec		w = nl;
		Vec		u = vecNorm(vecCross(fabs(w.x) > .1 ? vecNew(0,1,0) : vecNew(1,0,0), w));
		Vec		v = vecCross(w, u);
		Vec		d = vecNorm(vecAdd(vecAdd(vecScale(u, cos(r1) * r2s),
						vecScale(v, sin(r1) * r2s)),
					vecScale(w, sqrt(1-r2))));

		return vecAdd(emission, vecMult(f, radiance((Ray){x, d}, depth, Xi)));
	}
	/* Specular */
	else if (obj->refl == kReflSpec)
	{
		Vec		reflDir = vecSub(r.direction, vecScale(n, 2 * vecDot(n, r.direction)));
		return vecAdd(emission, vecMult(f, radiance((Ray){x, reflDir}, depth, Xi)));
	}

	/* Refraction - Snell's law and Fresnel effects */
	Ray		reflRay = {x, vecSub(r.direction, vecScale(n, 2*vecDot(n, r.direction)))};
	bool		into = vecDot(n, nl) > 0;
	double		nc = 1;
	double		nt = 1.5;
	double		nnt = into ? nc/nt : nt/nc;
	double		ddn = vecDot(r.direction, nl);
	double		cos2t;

	if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn)) < 0)
	{
		return vecAdd(emission, vecMult(f, radiance(reflRay, depth, Xi)));
	}

	Vec		tdir = vecNorm(vecSub(vecScale(r.direction, nnt),
					vecScale(n, (into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))));
	double		a = nt - nc;
	double		b = nt + nc;
	double		R0 = a*a/(b*b);
	double		cVal = 1 - (into ? -ddn : vecDot(tdir,n));
	double		Re = R0 + (1 - R0) * pow(cVal, 5);
	double		Tr = 1 - Re;
	double		P = .25 + .5*Re;
	double		RP = Re/P;
	double		TP = Tr/(1-P);

	if (depth > 2)
	{
		return vecAdd(emission, vecMult(f,
				rnd(Xi) < P
					? vecScale(radiance(reflRay, depth, Xi), RP)
					: vecScale(radiance((Ray){x, tdir}, depth, Xi), TP)));
	}
	else
	{
		return vecAdd(emission, vecMult(f,
				vecAdd(vecScale(radiance(reflRay, depth, Xi), Re),
				       vecScale(radiance((Ray){x, tdir}, depth, Xi), Tr))));
	}
}

/*
 *	Main - sets up image and loops over rows. Splits into subpixel for 
 *	anti-aliasing, loops over samples for Monte Carlo rays, accumulates
 *	radiance per subpixel
 */
int
main(int argc, char *	argv[])
{
	int	width = 1024;
	int	height = 768;
	int	samples = (argc >= 2) ? atoi(argv[1])/4 : 1;
	bool	preview = (argc >= 3 &&
			(strcmp(argv[2], "--preview")==0 || strcmp(argv[2], "-p")==0));

	int	startY = 0;
	int	endY = height;
	int	startX = 0;
	int	endX = width;

	if (preview)
	{
		startY = height/4;
		endY   = 3*height/4;
		startX = width/4;
		endX   = 3*width/4;
	}

	Ray	cam = {vecNew(50,52,295.6), vecNorm(vecNew(0,-0.042612,-1))};
	Vec	cx = vecNew(width*.5135/height, 0, 0);
	Vec	cy = vecScale(vecNorm(vecCross(cx, cam.direction)), .5135);
	Vec *	c = calloc(width*height, sizeof(Vec));

	initJitter();

	#pragma omp parallel for schedule(dynamic,1)
	for (int yy = startY; yy < endY; yy++)
	{
		fprintf(stderr,"\rRendering %d spp %5.2f%%", samples*4,
			100.*(yy-startY)/(endY-startY-1));
		unsigned short	Xi[3] = {0,0,(unsigned short)(yy*yy*yy)};

		for (int xx = startX; xx < endX; xx++)
		{
			int	ii = (height-yy-1)*width + xx;
			Vec	r = vecNew(0,0,0);

			for (int sy=0; sy<2; sy++)
			{
				for (int sx=0; sx<2; sx++)
				{
					Vec	sub = vecNew(0,0,0);

					for (int ss=0; ss<samples; ss++)
					{
						double	r1 = 2*rnd(Xi);
						double	dx = r1<1 ? sqrt(r1)-1 : 1-sqrt(2-r1);
						double	r2 = 2*rnd(Xi);
						double	dy = r2<1 ? sqrt(r2)-1 : 1-sqrt(2-r2);
						Vec	d = vecAdd(vecAdd(
									vecScale(cx, ((sx+.5+dx)/2 + xx)/width - .5),
									vecScale(cy, ((sy+.5+dy)/2 + yy)/height - .5)),
								cam.direction);
						Ray	ray = {vecAdd(cam.origin, vecScale(d,140)), vecNorm(d)};
						Vec	rad = radiance(ray,0,Xi);

						sub = vecAdd(sub, vecScale(rad, 1.0/samples));
					}

					r = vecAdd(r, vecScale(sub, .25));
				}
			}

			c[ii] = r;
		}
	}

	FILE *f = fopen("mountDir/image.ppm", "w");
	fprintf(f,"P3\n%d %d\n%d\n", width, height, 255);

	for (int ii=0; ii<width*height; ii++)
	{
		double mx = UxHwDoubleNthMoment(c[ii].x,1);
		double my = UxHwDoubleNthMoment(c[ii].y,1);
		double mz = UxHwDoubleNthMoment(c[ii].z,1);

		if (mx<0) mx=0; else if (mx>1) mx=1;
		if (my<0) my=0; else if (my>1) my=1;
		if (mz<0) mz=0; else if (mz>1) mz=1;

		int r = (int)(pow(mx,1/2.2)*255+.5);
		int g = (int)(pow(my,1/2.2)*255+.5);
		int b = (int)(pow(mz,1/2.2)*255+.5);

		fprintf(f,"%d %d %d ", r,g,b);
	}

	free(c);

	return 0;
}
