//Make : g++ IRTR.cc -o IRTR -O3 -fopenmp

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

const double inf = 1e31;

struct Vector {
  double x, y, z;
  Vector(const double &_x = 0., const double &_y = 0., const double &_z = 0.) {
    x = _x;
    y = _y;
    z = _z;
  }
  Vector operator+(const Vector &rhs) const {
    return Vector(x + rhs.x, y + rhs.y, z + rhs.z);
  }
  Vector operator-(const Vector &rhs) const {
    return Vector(x - rhs.x, y - rhs.y, z - rhs.z);
  }
  Vector operator*(const double &rhs) const {
    return Vector(x * rhs, y * rhs, z * rhs);
  }
  Vector operator*(const Vector &rhs) const { // cross
    return Vector(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x);
  }
  double operator/(const Vector &rhs) const { // dot
    return x * rhs.x + y * rhs.y + z * rhs.z;
  }
  Vector operator%(const Vector &rhs) const { // multiple
    return Vector(x * rhs.x, y * rhs.y, z * rhs.z);
  }
  Vector &normalize() {
    *this = *this * (1. / sqrt(x * x + y * y + z * z));
    return *this;
  }
};

struct Ray {
  Vector ori, dir;
  Ray(const Vector &_ori, const Vector &_dir) {
    ori = _ori;
    dir = _dir;
  }
};

enum ReflectionType {
  Specular, Glass, Diffuse
};
struct Sphere {
  double r;
  Vector p, e, c;
  ReflectionType t;
  static constexpr double eps = 1e-5;
  Sphere(const double &_radius, const Vector &_position,
         const Vector &_emission, const Vector &_colour,
         const ReflectionType &_refl_type) {
    r = _radius;
    p = _position;
    e = _emission;
    c = _colour;
    t = _refl_type;
  }
  double Intersect(const Ray &ray) const { // returns distance, inf if nohit
    const Vector po = ray.ori - p;
    const double a = ray.dir / ray.dir;
    const double b = ray.dir / po * 2;
    const double c = po / po - pow(r, 2);
    double delta = pow(b, 2) - a * c * 4;
    if (delta < 0)
      return inf;
    delta = sqrt(delta);
    const double t_1 = (- b - delta) / (a * 2);
    const double t_2 = (- b + delta) / (a * 2);
    if (t_1 > eps)
      return t_1;
    if (t_2 > eps)
      return t_2;
    return inf;
  }
};

namespace Arguments {
  int samp_num = 10;
  int w = 1024, h = 768;
  FILE *file = fopen("image.ppm", "w");

  bool frog = true;
  double rho = 0.007;
  const Vector frog_c(0.7, 0.7, 0.7); // Colour of frog

  double ipp = 0.07; // inch per pixel in sensor
  const Vector camera_x(1, 0, 0), camera_y(0, 1, 0); // after nomalized
  double lensr = 4, u = 80, v = u / 2.5;
  const Ray camera(Vector(50, 40, 150), Vector(0, 0, -1).normalize());
  const Vector lens_centre = camera.ori + camera.dir * v;
  const Sphere spheres[] = {
    Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(), Vector(.75, .25, .25), Diffuse), // Left
    Sphere(1e5, Vector(-1e5 + 99, 40.8, 81.6), Vector(), Vector(.25, .25, .75), Diffuse), // Rght
    Sphere(1e5, Vector(50, 40.8, 1e5), Vector(), Vector(.75, .75, .75), Specular), // Back
    Sphere(1e5, Vector(50, 40.8, -1e5 + 170), Vector(), Vector(), Diffuse), // Frnt
    Sphere(1e5, Vector(50, 1e5, 81.6), Vector(), Vector(.75, .75, .75), Diffuse), // Botm
    Sphere(1e5, Vector(50, -1e5 + 81.6, 81.6), Vector(.8, .8, .8), Vector(.75, .75, .75), Diffuse), // Top

    Sphere(10, Vector(73, 59, 10), Vector(), Vector(0, .9, .9), Diffuse), // Ball
    Sphere(10, Vector(40, 45, 22), Vector(), Vector(.4, .8, 0), Diffuse), // Ball
    Sphere(10, Vector(27, 30, 37), Vector(), Vector(.8, .8, .1), Diffuse), // Ball
    Sphere(10, Vector(50, 15, 50), Vector(), Vector(1, 1, 1) * 0.999, Glass), // Ball
    Sphere(10, Vector(77, 16.5, 68), Vector(), Vector(.9, .45, .15), Diffuse), // Ball
  };

  void Decode(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
      switch (argv[i][1]) {
        case 'n': {
          samp_num = atoi(argv[++i]);
          break;
        }
        case 'o': {
          file = fopen(argv[++i], "w");
        }
      }
    }
    if (!file) {
      fprintf(stderr, "Failed to open image file.");
      exit(1);
    }
  }
};
using namespace Arguments;

inline double f_atmo(const double &dis) {
  return exp(-pow(rho * dis, 2));
}
inline double clamp(const double &x) {
  if (x < 0)
    return 0;
  if (x > 1)
    return 1;
  return x;
}
inline int Gamma(const double &x) {
  return int(pow(clamp(x), 1 / 2.2) * 255 + .5); // Gamma Correction
}

int Intersect(const Ray &r, double &t) {
  int id;
  t = inf;
  double tmp;
  for (int i = int(sizeof(spheres) / sizeof(Sphere)) - 1; i >= 0; i--) {
    tmp = spheres[i].Intersect(r);
    if (tmp < t) {
      id = i;
      t = tmp;
    }
  }
  if (tmp >= inf)
    return -1;
  return id;
}
Vector Radiance(const Ray &r, int depth, unsigned short *Xi) {
  double t;
  int id = Intersect(r, t);
  if (id < 0)
    return Vector();
  const Sphere &obj = spheres[id]; // the hit object
  Vector col = obj.c;
  Vector point = r.ori + r.dir * t; // intersect point
  Vector n = (point - obj.p).normalize(); // normal
  Vector o_n = n / r.dir < 0 ? n : n * -1; // orinted normal
  double p = std::max(col.x, std::max(col.y, col.z)); // max reflection
  const bool in  = n / r.dir < 0;
  if (++depth > 5) { // Russian Roulette
    if (erand48(Xi) < p)
      col = col * (1 / p);
    else
      return obj.e;
  }
  if (obj.t == Specular) {
    if (frog) {
      const double p_frog = f_atmo(t);
      return obj.e + col % Radiance(Ray(point, r.dir - n * (2 * (n / r.dir))), depth, Xi) * p_frog + frog_c * (1 - p_frog);
    }
    return obj.e + col % Radiance(Ray(point, r.dir - n * (2 * (n / r.dir))), depth, Xi);
  } else if (obj.t == Diffuse) {
    const double ang_ = 2 * M_PI * erand48(Xi); // random angle
    const double dis_ = erand48(Xi), dis_i_ = sqrt(dis_); // random distance
    const Vector u = ((fabs(o_n.x) > .1 ? Vector(0, 1) : Vector(1)) * o_n).normalize(); // u $\perp$ o_n
    const Vector v = o_n * u; // v $\perp$ u && v $\perp$ o_n
    const Vector dir = (u * (cos(ang_) * dis_i_) + v * (sin(ang_) * dis_i_) + o_n * sqrt(1 - dis_)).normalize();
    if (frog) {
      const double p_frog = f_atmo(t);
      return obj.e + col % Radiance(Ray(point, dir), depth, Xi) * p_frog + frog_c * (1 - p_frog);
    }
    return obj.e + col % Radiance(Ray(point, dir), depth, Xi);
  } else if (obj.t == Glass) {
    const Ray refl_ray(point, r.dir - n * (2 * (n / r.dir)));
    const double n_air = 1, n_obj = 1.5, n_relative = in ? n_air / n_obj : n_obj / n_air;
    const double d_d_n = r.dir / o_n;
    const double cosr_2 = 1 - pow(n_relative, 2) * (1 - pow(d_d_n, 2));
    if (cosr_2 < 0) { // total internal reflection
      if (frog && in) {
        const double p_frog = f_atmo(t);
        return obj.e + col % Radiance(refl_ray, depth, Xi) * p_frog + frog_c * (1 - p_frog);
      }
      return obj.e + col % Radiance(refl_ray, depth, Xi);
    } else {
      const Vector t_dir = (r.dir * n_relative - n * ((in ? 1 : -1) * (d_d_n * n_relative + sqrt(cosr_2)))).normalize();
      double a = n_relative - 1, b = n_relative + 1, F_0 = pow(a, 2) / pow(b, 2);
      double Re = F_0 + (1 - F_0) * pow(1 - (in ? -d_d_n : t_dir / n), 5), Tr = 1 - Re, P = .25 + .5 * Re; // Fresnel Reflectance
      const Vector radiance = (depth > 2 ?
                               (erand48(Xi) < P ? Radiance(refl_ray, depth, Xi) * (Re / P)
                                                  : Radiance(Ray(point, t_dir), depth, Xi) * (Tr / (1 - P)))
                               : Radiance(refl_ray, depth, Xi) * Re + Radiance(Ray(point, t_dir), depth, Xi) * Tr);
      if (frog && in) {
        const double p_frog = f_atmo(t);
        return obj.e + col % radiance * p_frog + frog_c * (1 - p_frog);
      }
      return obj.e + col % radiance;
    }
  }
}

int main(int argc, char *argv[]) {
  Decode(argc, argv);
  Vector colour;
  Vector *map = new Vector[w * h];
  #pragma omp parallel for schedule(dynamic, 1) private(colour)
    for (int y = 0; y < h; y++) {
      fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samp_num, 100. * y / (h - 1));
      for (unsigned short x = 0, Xi[3] = {0, 0, (unsigned short)(y * y * y)}; x < w; x++, colour = Vector()) {
        const int id = (h - y - 1) * w + x;
        const Vector sensor_point = camera.ori + camera_x * (w / 2 - x) * ipp + camera_y * (h / 2 - y) * ipp;
        const Vector focus_point = lens_centre + (lens_centre - sensor_point) * (u / v);
        for (int i = 0; i < samp_num; i++) {
          const double theta = 2 * M_PI * erand48(Xi), radius = lensr * erand48(Xi);
          const Vector lens_point = lens_centre + camera_x * (cos(theta) * radius) + camera_y * (sin(theta) * radius);
          colour = colour + Radiance(Ray(lens_point, (focus_point - lens_point).normalize()), 0, Xi) * (1. / samp_num);
        }
        map[id] = Vector(clamp(colour.x), clamp(colour.y), clamp(colour.z));
      }
    }
  fprintf(stderr, "\n");
  fprintf(file, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i = 0; i < w * h; i++)
    fprintf(file, "%d %d %d ", Gamma(map[i].x), Gamma(map[i].y), Gamma(map[i].z));
  return 0;
}