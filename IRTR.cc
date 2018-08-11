//Make : g++ IRTR.cc -o IRTR -O3 -fopenmp

#include <map>
#include <cmath>
#include <cstdio>
#include <cassert>
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
  int w = 1024, h = 768;
  double angle = .5135;
  int samp_num = 10;
  int dof = 10 / 2, l = 30, k = 0.7;
  FILE *file = fopen("image.ppm", "w");

  bool frog = false;
  double rho = 0.007;
  const Vector frog_c(0.7, 0.7, 0.7); // Colour of frog
  const Ray camera(Vector(50, 52, 295.6), Vector(0, -0.042612, -1).normalize());
  Vector tc(0.0588, 0.361, 0.0941);
  Vector sc = Vector(1,1,1)*.7;
  const Sphere spheres[] = {
    // center 50 40.8 62
    // floor 0
    // back  0
    Sphere(1e5, Vector(50, 1e5+130, 0),  Vector(1,1,1)*1.3,Vector(),Diffuse), //lite
    Sphere(1e2, Vector(50, -1e2+2, 47),  Vector(),Vector(1,1,1)*.7,Diffuse), //grnd

    Sphere(1e4, Vector(50, -30, 300)+Vector(-sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4, Vector(), Vector(1,1,1)*.99,Specular),// mirr L
    Sphere(1e4, Vector(50, -30, 300)+Vector(sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4,  Vector(), Vector(1,1,1)*.99,Specular),// mirr R
    Sphere(1e4, Vector(50, -30, -50)+Vector(-sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4,Vector(), Vector(1,1,1)*.99,Specular),// mirr FL
    Sphere(1e4, Vector(50, -30, -50)+Vector(sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4, Vector(), Vector(1,1,1)*.99,Specular),// mirr


    Sphere(4, Vector(50,6*.6,47),   Vector(),Vector(.13,.066,.033), Diffuse),//"tree"
    Sphere(16,Vector(50,6*2+16*.6,47),   Vector(), tc,  Diffuse),//"tree"
    Sphere(11,Vector(50,6*2+16*.6*2+11*.6,47),   Vector(), tc,  Diffuse),//"tree"
    Sphere(7, Vector(50,6*2+16*.6*2+11*.6*2+7*.6,47),   Vector(), tc,  Diffuse),//"tree"

    Sphere(15.5,Vector(50,1.8+6*2+16*.6,47),   Vector(), sc,  Diffuse),//"tree"
    Sphere(10.5,Vector(50,1.8+6*2+16*.6*2+11*.6,47),   Vector(), sc,  Diffuse),//"tree"
    Sphere(6.5, Vector(50,1.8+6*2+16*.6*2+11*.6*2+7*.6,47),   Vector(), sc,  Diffuse),//"tree"

    // Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(), Vector(.75, .25, .25), Diffuse), // Left
    // Sphere(1e5, Vector(-1e5 + 99, 40.8, 81.6), Vector(), Vector(.25, .25, .75), Diffuse), // Rght
    // Sphere(1e5, Vector(50, 40.8, 1e5), Vector(), Vector(.75, .75, .75), Diffuse), // Back
    // Sphere(1e5, Vector(50, 40.8, -1e5 + 170), Vector(), Vector(), Diffuse), // Frnt
    // Sphere(1e5, Vector(50, 1e5, 81.6), Vector(), Vector(.75, .75, .75), Diffuse), // Botm
    // Sphere(1e5, Vector(50, -1e5 + 81.6, 81.6), Vector(), Vector(.75, .75, .75), Diffuse), // Top
    // Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(1, 1, 1) * .999, Specular), // Mirr
    // Sphere(16.5, Vector(73, 16.5, 78), Vector(), Vector(1, 1, 1) * .999, Glass), // Glas
    // Sphere(600, Vector(50, 681.6 - .27, 81.6), Vector(12, 12, 12), Vector(), Diffuse) // Lite
  };

  void Decode(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
      switch (argv[i][1]) {
        case 'n': {
          samp_num = atoi(argv[++i]) / 4;
          break;
        }
        case 'o': {
          file = fopen(argv[++i], "w");
        }
      }
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
inline int GetRadius(const double &dis) {
  if (dis < l - dof)
    return (dof - l - dis) * k;
  else if (dis <= l + dof)
    return 0;
  return (dis - l - dof) * k;
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
Vector Radiance(const Ray &r, int depth, unsigned short *Xi, double &dis) {
  double t;
  int id = Intersect(r, t);
  if (id < 0) {
    dis = inf;
    return Vector();
  }
  if (dis < 0)
    dis = t;
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
      return obj.e + col % Radiance(Ray(point, r.dir - n * 2 * (n / r.dir)), depth, Xi, dis) * p_frog + frog_c * (1 - p_frog);
    }
    return obj.e + col % Radiance(Ray(point, r.dir - n * 2 * (n / r.dir)), depth, Xi, dis);
  } else if (obj.t == Diffuse) {
    const double ang_ = 2 * M_PI * erand48(Xi); // random angle
    const double dis_ = erand48(Xi), dis_i_ = sqrt(dis_); // random distance
    const Vector u = ((fabs(o_n.x) > .1 ? Vector(0, 1) : Vector(1)) * o_n).normalize(); // u $\perp$ o_n
    const Vector v = o_n * u; // v $\perp$ u && v $\perp$ o_n
    const Vector dir = (u * cos(ang_) * dis_i_ + v * sin(ang_) * dis_i_ + o_n * sqrt(1 - dis_)).normalize();
    if (frog) {
      const double p_frog = f_atmo(t);
      return obj.e + col % Radiance(Ray(point, dir), depth, Xi, dis) * p_frog + frog_c * (1 - p_frog);
    }
    return obj.e + col % Radiance(Ray(point, dir), depth, Xi, dis);
  } else if (obj.t == Glass) {
    const Ray refl_ray(point, r.dir - n * 2 * (n / r.dir));
    const double n_air = 1, n_obj = 1.5, n_relative = in ? n_air / n_obj : n_obj / n_air;
    const double d_d_n = r.dir / o_n;
    const double cos_2t = 1 - pow(n_relative, 2) * (1 - pow(d_d_n, 2));
    if (cos_2t < 0) { // total internal reflection
      if (frog && in) {
        const double p_frog = f_atmo(t);
        return obj.e + col % Radiance(refl_ray, depth, Xi, dis) * p_frog + frog_c * (1 - p_frog);
      }
      return obj.e + col % Radiance(refl_ray, depth, Xi, dis);
    } else {
      const Vector t_dir = (r.dir * n_relative - n * (in ? 1 : -1) * (d_d_n * n_relative + sqrt(cos_2t))).normalize();
      double a = n_obj - n_air, b = n_obj + n_air, R0 = pow(a, 2) / pow(b, 2), c = 1 - (in ? -d_d_n : t_dir / n);
      double Re = R0 + (1 - R0) * pow(c, 5), Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
      const Vector radiance = (depth > 2 ?
                               (erand48(Xi) < P ? Radiance(refl_ray, depth, Xi, dis) * RP
                                                  : Radiance(Ray(point, t_dir), depth, Xi, dis) * TP)
                               : Radiance(refl_ray, depth, Xi, dis) * Re + Radiance(Ray(point, t_dir), depth, Xi, dis) * Tr);
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
  double dis;
  Vector colour;
  Vector *map = new Vector[w * h];
  double *depths = new double[w * h];
  Vector *image = new Vector[w * h];
  double *sum_p = new double[w * h];
  Vector cx = Vector(w * angle / h);
  Vector cy = (cx * camera.dir).normalize() * angle;
  std::map<int, double>  table; // the radius of the circle, 1 / the number of the pixels which are covered
  #ifndef DEBUG
    #pragma omp parallel for schedule(dynamic, 1) private(colour, dis)
  #endif // !DEBUG
    for (int y = 0; y < h; y++) {
      fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samp_num, 100. * y / (h - 1));
      for (unsigned short x = 0, Xi[3] = {0, 0, (unsigned short)(y * y * y)}; x < w; x++) {
        const int id = (h - y - 1) * w + x;
        for (int sy = 0; sy < 2; sy++) {
          for (int sx = 0; sx < 2; sx++, colour = Vector()) {
            for (int s = 0; s < samp_num; s++, dis = -1) {
              double r1 = 2 * erand48(Xi);
              double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
              double r2 = 2 * erand48(Xi);
              double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
              Vector d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                        cy * (((sy + .5 + dy) / 2 + y) / h - .5) + camera.dir;
              colour = colour + Radiance(Ray(camera.ori + d * 140, d.normalize()), 0, Xi, dis)
                      * (1. / samp_num);
              depths[id] += dis * (1. / samp_num);
            }
            map[id] = map[id] + Vector(clamp(colour.x), clamp(colour.y), clamp(colour.z)) * .25;
          }
        }
        double p;
        const int r = GetRadius(depths[id]), r_2 = r * r;
        if (table.count(r)) {
          p = table[r];
        } else {
          int cnt = 0;
          for (int i = 1; i <= r; i++) {
            for (int j = 1; j <= i; j++) {
              if (i * i + j * j > r_2)
                continue;
              if (i == j)
                cnt += 4;
              else
                cnt += 8;
            }
          }
          cnt += 4 * r + 1;
          table[r] = p = 1. / cnt;
        }
        for (int i = std::max(0, x - r); i <= std::min(w - 1, x + r); i++) {
          for (int j = std::max(0, y - r); j <= std::min(h - 1, y + r); j++) {
            if (pow(x - i, 2) + pow(y - j, 2) > r_2)
              continue;
            sum_p[(h - j - 1) * w + i] += p;
          }
        }
      }
    }

  fprintf(file, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i = 0; i < w * h; i++)
    fprintf(file, "%d %d %d ", Gamma(map[i].x), Gamma(map[i].y), Gamma(map[i].z));

  fprintf(stderr, "\n");
  #ifndef DEBUG
    #pragma omp parallel for schedule(dynamic, 1)
  #endif // !DEBUG
  for (int y = 0; y < h; y++) {
    fprintf(stderr, "\rCreating %5.2f%%", 100. * y / (h - 1));
    for (int x = 0; x < w; x++) {
      const int id = (h - y - 1) * w + x;
      const int r = GetRadius(depths[id]), r_2 = r * r;
      const double p = table[r];
      for (int i = std::max(0 , x - r); i <= std::min(w - 1, x + r); i++) {
        for (int j = std::max(0, y - r), id_ = (h - j - 1) * w + i; j <= std::min(h - 1, y + r); j++) {
          if (pow(x - i, 2) + pow(y - j, 2) > r_2)
            continue;
          image[id_] = image[id_] + map[id] * (p / sum_p[id_]);
        }
      }
    }
  }

  // fprintf(file, "P3\n%d %d\n%d\n", w, h, 255);
  // for (int i = 0; i < w * h; i++)
  //   fprintf(file, "%d %d %d ", Gamma(image[i].x), Gamma(image[i].y), Gamma(image[i].z));

  return 0;
}