//Make : g++ IRTR.cc -o IRTR -O3 -fopenmp

#include <cmath>
#include <cstdio>
#include <cstdlib>

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
  double Intersect(const Ray &ray) { // returns distance, -1 if nohit
    const Vector o_minus_p = ray.ori - p;
    const double a = ray.dir / ray.dir;
    const double b = ray.dir / o_minus_p * 2;
    const double c = o_minus_p / o_minus_p - r;
    double delta = b * b - a * c * 4;
    if (delta < 0)
      return -1;
    delta = sqrt(delta);
    const double t_1 = (- b - delta) / (a * 2);
    const double t_2 = (- b + delta) / (a * 2);
    if (t_1 > eps)
      return t_1;
    if (t_2 > eps)
      return t_2;
    return -1;
  }
};

// Scene Settings
const Sphere spheres[] = {
  Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(),
         Vector(.75, .25, .25), Specular)
};
const Ray camera(Vector(50, 52, 295.6),
                 Vector(0, -0.042612, -1).normalize());
//

inline double clamp(double x) {
  if (x < 0)
    return 0;
  if (x > 1)
    return 1;
  return x;
}
inline int toInt(double x) {
  return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

Vector Radiance(const Ray &r, int depth, unsigned short *Xi) {

}

int main(int argc, char *argv[]) {
  int w = 1024, h = 768;
  double lens = .5135;
  int samp_num = 10;
  double f = 5.6;
  FILE *file;
  for (int i = 1; i < argc; i++) {
    switch (argv[i][1]) {
      case 's': {
        w = atoi(argv[++i]);
        h = atoi(argv[++i]);
        break;
      }
      case 'l': {
        lens = atof(argv[++i]);
        break;
      }
      case 'n': {
        samp_num = atoi(argv[++i]) / 4;
        break;
      }
      case 'f': {
        f = atof(argv[++i]);
        break;
      }
      case 'o': {
        file = fopen(argv[++i], "w");
      }
    }
  }
  Vector colour;
  Vector *map = new Vector[w * h];
  Vector cx = Vector(w * lens / h);
  Vector cy = (cx * camera.dir).normalize() * lens;
  #pragma omp parallel for schedule(dynamic, 1) private(colour)
    for (int y = 0; y < h; y++) {
    fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samp_num,
            100. * y / (h - 1));
    for (unsigned short x = 0, Xi[3] = {0, 0, (unsigned short)(y * y * y)}; x < w; x++) {
      for (int sy = 0, id = (h - y - 1) * w + x; sy < 2; sy++)
        for (int sx = 0; sx < 2; sx++, colour = Vector()) {
          for (int s = 0; s < samp_num; s++) {
            double r1 = 2 * erand48(Xi);
            double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
            double r2 = 2 * erand48(Xi);
            double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
            Vector d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                       cy * (((sy + .5 + dy) / 2 + y) / h - .5) + camera.dir;
            colour = colour + Radiance(Ray(camera.ori + d * 140, d.normalize()), 0, Xi)
                     * (1. / samp_num);
          }
          map[id] = map[id] + Vector(clamp(colour.x), clamp(colour.y), clamp(colour.z)) * .25;
        }
      }
    }

  fprintf(file, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i = 0; i < w * h; i++)
    fprintf(file, "%d %d %d ", toInt(map[i].x), toInt(map[i].y), toInt(map[i].z));
  return 0;
}