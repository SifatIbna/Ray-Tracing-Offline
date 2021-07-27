#include <bits/stdc++.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #include <windows.h>
// #include <GL/glut.h>

#include "include/glut.h"

#include <vector>

using namespace std;

#define SLICES 100
#define STACKS 100

#define AMB 0
#define DIFF 1
#define SPEC 2
#define REF 3

int recurLevel;

class point
{
public:
  double x, y, z;
  point(double a, double b, double c)
  {
    x = a;
    y = b;
    z = c;
  }
  point()
  {
  }
  double len()
  {
    return sqrt(x * x + y * y + z * z);
  }
  point operator*(point p);
  point operator*(double) const;
  point operator+(point point);
  point operator-(point);
};
point point::operator-(point rightVec)
{
  return point(this->x - rightVec.x, this->y - rightVec.y, this->z - rightVec.z);
}
point point::operator+(point pnt)
{
  return point(this->x + pnt.x, this->y + pnt.y, this->z + pnt.z);
}
point point::operator*(point p)
{
  return point(this->x * p.x, this->y * p.y, this->z * p.z);
}

point point::operator*(double scalar) const
{
  return point(this->x * scalar, this->y * scalar, this->z * scalar);
}

/// For camera part

point Normalized(point &a)
{
  double val = (a.x * a.x + a.y * a.y + a.z * a.z);
  val = sqrt(val);
  a.x /= val;
  a.y /= val;
  a.z /= val;
  return a;
}
void Rotate(point &a, point &b, double angle)
{
  a.x = a.x * cos(angle) + b.x * sin(angle);
  a.y = a.y * cos(angle) + b.y * sin(angle);
  a.z = a.z * cos(angle) + b.z * sin(angle);
  Normalized(a);
}
void Multiplication(point &a, point &b, point &c)
{

  a.x = (b.y * c.z - b.z * c.y);
  a.y = (b.z * c.x - b.x * c.z);
  a.z = (b.x * c.y - b.y * c.x);
  Normalized(a);
}
double Dot(point a, point b)
{
  return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

point Cross(point a, point b)
{
  point c;
  c.x = (a.y * b.z) - (a.z * b.y);
  c.y = (a.z * b.x) - (a.x * b.z);
  c.z = (a.x * b.y) - (a.y * b.x);
  return c;
}
///end for camera part
class Light
{
public:
  point pos;
  point color;
  Light() {}
  Light(point pos, point col)
  {
    this->pos = pos;
    this->color = col;
  }
};

class Ray
{
public:
  point start;
  point dir;
  Ray() {}
  Ray(point s, point d)
  {
    start = s;
    dir = d;
  }
};

class Object
{
public:
  int triangle, floor, sphere, general;
  point reference_point;
  double height, width, length;

  double color[3];
  double coefficients[4];
  int shine;

  Object() {}
  virtual point getNormal(point pointed)
  {
    point p(0, 0, 0);
    return p;
  }
  void setShine(int shine)
  {
    this->shine = shine;
  };

  void setColor(double red, double green, double blue)
  {
    this->color[0] = red;
    this->color[1] = green;
    this->color[2] = blue;
  }

  void setCoEfficients(double amb, double diff, double spec, double ref)
  {
    this->coefficients[0] = amb;
    this->coefficients[1] = diff;
    this->coefficients[2] = spec;
    this->coefficients[3] = ref;
  }

  point getReflection(Ray *r, point normal)
  {
    //* Lighting and shading slide 23
    //* r = a  -  2 ( a . n) n

    point refelection = r->dir - (normal * Dot(r->dir, normal)) * 2.0;
    return Normalized(refelection);
  }

  virtual void draw() {}
  virtual double getIntersecting(Ray *r)
  {
    return -1;
  }
  virtual double intersect(Ray *r, point *Color, int level)
  {
    return -1;
  }
};

vector<Light> lights;
vector<Object *> objects;

double clipValue(double value, double minVal, double maxVal)
{
  if (value < minVal)
  {
    value = minVal;
  }
  if (value > maxVal)
  {
    value = maxVal;
  }

  return value;
}

bool isOutside(double value, double min_val, double max_val)
{
  return (value < min_val) || (value > max_val);
}

void Fill_Color_Obj(Object *obj, Ray *r, point *Color, int level, double t)
{

  //* Ray Tracing and ray Casting Slide 10
  //* P(t) = origin + t * direction
  point normal;
  point intersect = r->start + r->dir * t;

  //* intersectionPointColor = getColorAt(intersectionPoint)
  //* color = intersectionPointColor*coEfficient[AMB]

  Color->x = obj->color[0] * obj->coefficients[AMB];
  Color->y = obj->color[1] * obj->coefficients[AMB];
  Color->z = obj->color[2] * obj->coefficients[AMB];

  if (obj->triangle != 1)
    normal = obj->getNormal(intersect);
  else
    normal = obj->getNormal(r->dir);

  point reflect = obj->getReflection(r, normal);

  for (int i = 0; i < lights.size(); i++)
  {
    point direction = lights[i].pos - intersect;
    point start = intersect + direction * 0.000001;

    double len = direction.len();
    direction = Normalized(direction);

    Ray L_vec(start, direction);

    bool shadow = false;

    for (int j = 0; j < objects.size(); j++)
    {
      double t = objects[j]->getIntersecting(&L_vec);
      if (t > 0)
      {
        shadow = true;
        break;
      }
    }

    if (!shadow)
    {
      //* Lambert = L.N
      //* color and shading slide 20
      double lambert = Dot(L_vec.dir, normal);

      //* Phong = R.V
      //* color and shading slide 20
      double phong = Dot(reflect, r->dir);

      lambert = clipValue(lambert, 0.0, 1.0);
      phong = clipValue(phong, 0.0, 1.0);

      //* color += l.color*coEfficient[DIFF]*lambertValue* intersectionPointColor

      Color->x += lights[i].color.x * obj->color[0] * obj->coefficients[DIFF] * lambert;
      Color->y += lights[i].color.y * obj->color[1] * obj->coefficients[DIFF] * lambert;
      Color->z += lights[i].color.z * obj->color[2] * obj->coefficients[DIFF] * lambert;

      //* color += l.color*coEfficient[SPEC]*phongValueshine * intersectionPointColor

      Color->x += lights[i].color.x * obj->color[0] * obj->coefficients[SPEC] * pow(phong, obj->shine);
      Color->y += lights[i].color.y * obj->color[1] * obj->coefficients[SPEC] * pow(phong, obj->shine);
      Color->z += lights[i].color.z * obj->color[2] * obj->coefficients[SPEC] * pow(phong, obj->shine);
    }

    if (level < recurLevel)
    {
      point start = intersect + reflect * 0.000001; // slight ahead, 0.000001
      Ray reflectionRay(start, reflect);
      int nearest = -1;

      point reflectColor;
      double t_min = 999999999;
      for (int k = 0; k < objects.size(); k++)
      {
        double t = objects[k]->getIntersecting(&reflectionRay);

        if (t <= 0)
          continue;
        else if (t < t_min)
        {
          t_min = t;
          nearest = k;
        }
      }

      if (nearest != -1)
      {
        double t = objects[nearest]->intersect(&reflectionRay, &reflectColor, level + 1);
        if (t != -1)
        {
          //* color += colorreflected * coEfficient[REC_REFFLECTION]

          Color->x += reflectColor.x * obj->coefficients[REF];
          Color->y += reflectColor.y * obj->coefficients[REF];
          Color->z += reflectColor.z * obj->coefficients[REF];
        }
      }
    }

    Color->x = clipValue(Color->x, 0.0, 1.0);
    Color->y = clipValue(Color->y, 0.0, 1.0);
    Color->z = clipValue(Color->z, 0.0, 1.0);
  }
}

class Sphere : public Object
{
public:
  Sphere(point Center, double Radius)
  {
    reference_point = Center;
    height = length = Radius;

    sphere = 1;
    triangle = floor = general = 0;
  }

  void draw()
  {
    glColor3f(color[0], color[1], color[2]);
    glPushMatrix();
    glTranslatef(reference_point.x, reference_point.y, reference_point.z);
    glutSolidSphere(length, SLICES, STACKS);
    glPopMatrix();
  }

  double getIntersecting(Ray *r)
  {
    //* Ray tracing and casting slide 34

    point p = r->start - reference_point;

    double a = Dot(r->dir, r->dir);
    double b = 2.0 * Dot(r->dir, p);
    double c = Dot(p, p) - pow(length, 2);
    double D = b * b - 4 * a * c;

    if (D < 0)
      return -1;

    double t1 = (-b + sqrt(D)) / (2 * a);
    double t2 = (-b - sqrt(D)) / (2 * a);

    return min(t1, t2);
  }

  point getNormal(point pointed)
  {
    //* ray tracing and casting slide 39
    point normal = pointed - reference_point;
    return Normalized(normal);
  }

  double intersect(Ray *r, point *Color, int level)
  {

    double t = getIntersecting(r);

    if (t <= 0)
      return -1;

    if (level == 0)
      return t;

    Fill_Color_Obj(this, r, Color, level, t);

    return t;
  }
};

class Floor : public Object
{
public:
  double floor_width, tile_width;
  Floor() {}
  Floor(point p, double base)
  {
    reference_point = p;
    length = base;
    floor = 1;
    sphere = triangle = general = 0;
  }

  Floor(double floor_width, double tile_width)
  {

    this->floor_width = floor_width;
    this->tile_width = tile_width;

    reference_point.x = -floor_width / 2.0;
    reference_point.y = -floor_width / 2.0;
    reference_point.z = 0;
    length = tile_width;
    height = 0;

    floor = 1;
    sphere = triangle = general = 0;
  }

  point getNormal(point pointed)
  {
    pointed.x = 0;
    pointed.y = 0;
    pointed.z = 1;
    return pointed;
  }

  void draw()
  {

    int nTiles = this->floor_width / this->tile_width;

    glBegin(GL_QUADS);
    for (int i = 0; i < nTiles; i++)
    {
      for (int j = 0; j < nTiles; j++)
      {
        bool b = (i + j) % 2;
        glColor3f(b, b, b);

        glVertex3f(reference_point.x + this->tile_width * i, reference_point.y + this->tile_width * j, reference_point.z);
        glVertex3f(reference_point.x + this->tile_width * (i + 1), reference_point.y + this->tile_width * j, reference_point.z);
        glVertex3f(reference_point.x + this->tile_width * (i + 1), reference_point.y + this->tile_width * (j + 1), reference_point.z);
        glVertex3f(reference_point.x + this->tile_width * i, reference_point.y + this->tile_width * (j + 1), reference_point.z);
      }
    }
    glEnd();
  }

  double getIntersecting(Ray *r)
  {
    //* Ray tracing slide 24
    //* t = -(D+n.Ro)/n.Rd

    point normal = getNormal(reference_point);
    return ((-1.0) * Dot(r->start, normal) / Dot(r->dir, normal));
  }

  double intersect(Ray *r, point *Color, int level)
  {
    double t = getIntersecting(r);

    //* P(t) = Ro + Rd * t
    point intersect = r->start + r->dir * t;

    if (isOutside(intersect.x, reference_point.x, -reference_point.x) || isOutside(intersect.y, reference_point.y, -reference_point.y))
    {
      return -1;
    }
    else if (level != 0)
    {

      int x = (intersect.x - reference_point.x) / length;
      int y = (intersect.y - reference_point.y) / length;

      if ((x + y) % 2 == 0)
      {
        for (int i = 0; i < 3; i++)
        {
          color[i] = 0;
        }
      }
      else
      {
        for (int i = 0; i < 3; i++)
        {
          color[i] = 1;
        }
      }

      Fill_Color_Obj(this, r, Color, level, t);
      return t;
    }
    else
    {
      return t;
    }
  }
};

class Triangle : public Object
{
public:
  point a, b, c;
  Triangle(point a, point b, point c)
  {
    this->a = a;
    this->b = b;
    this->c = c;
    triangle = 1;
    sphere = floor = 0;
  }

  void draw()
  {
    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_TRIANGLES);
    {
      glVertex3f(a.x, a.y, a.z);
      glVertex3f(b.x, b.y, b.z);
      glVertex3f(c.x, c.y, c.z);
    }
    glEnd();
  }

  point getNormal(point pointed)
  {
    //* Normal = cross product of two vectors along the edges e.g. (b-a)X(c-a)
    point normal;
    point u(b.x - a.x, b.y - a.y, b.z - a.z);

    point v(c.x - a.x, c.y - a.y, c.z - a.z);
    normal = Cross(u, v);
    normal = Normalized(normal);

    if (Dot(pointed, normal) > 0)
    {
      normal.x = -normal.x;
      normal.y = -normal.y;
      normal.z = -normal.z;
    }

    return normal;
  }

  double getIntersecting(Ray *r)
  {
    //* Möller–Trumbore intersection algorithm
    //* From Wikipedia https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    //* Link was given in the spec

    const float EPSILON = 0.0000001;

    point edge1, edge2, h, s, q;
    float a1, f, u, v;

    edge1 = b - a;
    edge2 = c - a;

    h = Cross(r->dir, edge2);

    a1 = Dot(edge1, h);
    if (a1 > -EPSILON && a1 < EPSILON)
      return -1; // This ray is parallel to this triangle.
    f = 1.0 / a1;
    s = r->start - a;
    u = f * Dot(s, h);
    if (u < 0.0 || u > 1.0)
      return -1;
    q = Cross(s, edge1);
    v = f * Dot(r->dir, q);
    if (v < 0.0 || u + v > 1.0)
      return -1;
    float t = f * Dot(edge2, q);
    if (t > EPSILON) // ray intersection
    {
      return t;
    }
    else // This means that there is a line intersection but not a ray intersection.
      return -1;
  }
  /*
  double getIntersecting(Ray *r)
  {

    point edge1, edge2, edge3;
    edge1.x = a.x - b.x;
    edge1.y = a.y - b.y;
    edge1.z = a.z - b.z;
    // Normalized(edge1);

    edge2.x = a.x - c.x;
    edge2.y = a.y - c.y;
    edge2.z = a.z - c.z;
    // Normalized(edge2);
    edge3.x = a.x - r->start.x;
    edge3.y = a.y - r->start.y;
    edge3.z = a.z - r->start.z;

    double A = 0.0;

    A += (edge1.x * (edge2.y * r->dir.z - r->dir.y * edge2.z));
    A -= (edge1.y * (edge2.x * r->dir.z - r->dir.x * edge2.z));
    A += (edge1.z * (edge2.x * r->dir.y - r->dir.x * edge2.y));

    if (A < eps)
    {
      return -1;
    }

    double Gamma = 0.0;

    Gamma += (edge1.x * (edge3.y * r->dir.z - r->dir.y * edge3.z));
    Gamma -= (edge1.y * (edge3.x * r->dir.z - r->dir.x * edge3.z));
    Gamma += (edge1.z * (edge3.x * r->dir.y - r->dir.x * edge3.y));
    Gamma /= A;
    if (Gamma < eps || Gamma > 1.0)
    {
      return -1;
    }

    double Beta = 0.0;

    Beta += (edge3.x * (edge2.y * r->dir.z - r->dir.y * edge2.z));
    Beta -= (edge3.y * (edge2.x * r->dir.z - r->dir.x * edge2.z));
    Beta += (edge3.z * (edge2.x * r->dir.y - r->dir.x * edge2.y));
    Beta /= A;
    if (Beta < eps || Beta > 1 - Gamma)
    {
      return -1;
    }

    double t = 0.0;

    t += (edge1.x * (edge2.y * edge3.z - edge3.y * edge2.z));
    t -= (edge1.y * (edge2.x * edge3.z - edge3.x * edge2.z));
    t += (edge1.z * (edge2.x * edge3.y - edge3.x * edge2.y));
    t /= A;

    return t;
  }
*/
  double intersect(Ray *r, point *Color, int level)
  {
    double t = getIntersecting(r);
    if (t <= 0)
    {
      return -1;
    }
    if (level == 0)
    {
      return t;
    }

    Fill_Color_Obj(this, r, Color, level, t);
    return t;
  }
};

class General : public Object
{
public:
  double A, B, C, D, E, F, G, H, I, J;

  General(double coeff[10], point reference_point, double length, double width, double height)
  {
    this->reference_point = reference_point;
    this->height = height;
    this->width = width;
    this->length = length;
    this->A = coeff[0];
    this->B = coeff[1];
    this->C = coeff[2];
    this->D = coeff[3];
    this->E = coeff[4];
    this->F = coeff[5];
    this->G = coeff[6];
    this->H = coeff[7];
    this->I = coeff[8];
    this->J = coeff[9];
    general = 1;
    sphere = triangle = floor = 0;
  }

  void draw() {}

  point getNormal(point p)
  {
    //*Ray tracing and Ray casting slide 55
    //* 2Ax + Dy + Ez + G
    double u = 2 * A * p.x + D * p.y + E * p.z + G;

    //* 2By + Dx + Fz + H
    double v = 2 * B * p.y + D * p.x + F * p.z + H;

    //* 2Cz + Ex + Fy + I
    double z = 2 * C * p.z + E * p.x + F * p.y + I;

    point normal(u, v, z);
    return Normalized(normal);
  }

  double getIntersecting(Ray &ray)
  {
    //*Ray tracing and Ray casting slide 55

    //* General Equation
    //* t^2*a + t* b + c = 0

    //* t^2 * (A * xd^2+ B * yd^2 + C * zd^2 + D* xd* yd + E* yd * zd + F * xd * zd)
    double a = A * pow(ray.dir.x, 2) + B * pow(ray.dir.y, 2) + C * pow(ray.dir.z, 2) + (D * ray.dir.x * ray.dir.y + E * ray.dir.y * ray.dir.z + F * ray.dir.z * ray.dir.x);

    //* t* (2A* xo * xd + 2B * yo * yd + 2C* Zo* Zd + D* xo * yd+ D* yo* xd+ E* yo * zd+ E* yd* zo+ F*zo*xd + F* xo * zd + G*xd + H* yd + I * zd)
    double b = 2 * (A * ray.start.x * ray.dir.x + B * ray.start.y * ray.dir.y + C * ray.start.z * ray.dir.z) + (D * (ray.start.x * ray.dir.y + ray.dir.x * ray.start.y) + E * (ray.start.y * ray.dir.z + ray.dir.y * ray.start.z) + F * (ray.start.z * ray.dir.x + ray.dir.z * ray.start.x)) + (G * ray.dir.x + H * ray.dir.y + I * ray.dir.z);

    //* A * xo ^2 + B * yo ^2 + C * zo^2 + D * xo * yo + E * xo * zo + F * zo * xo + G * xo + H *yo + Izo + J
    double c = A * pow(ray.start.x, 2) + B * pow(ray.start.y, 2) + C * pow(ray.start.z, 2) + (D * ray.start.x * ray.start.y + E * ray.start.y * ray.start.z + F * ray.start.z * ray.start.x) + (G * ray.start.x + H * ray.start.y + I * ray.start.z + J);

    double Des = b * b - 4 * a * c;

    if (Des < 0)
    {
      return -1;
    }

    double t1 = (-b + sqrt(Des)) / (2.0 * a);
    double t2 = (-b - sqrt(Des)) / (2.0 * a);

    point intersectionPoint1 = ray.start + ray.dir * t1;
    point intersectionPoint2 = ray.start + ray.dir * t2;

    double x_min = reference_point.x;
    double x_max = x_min + length;

    double y_min = reference_point.y;
    double y_max = y_min + width;

    double z_min = reference_point.z;
    double z_max = z_min + height;

    double x1 = intersectionPoint1.x;
    double y1 = intersectionPoint1.y;
    double z1 = intersectionPoint1.z;

    double x2 = intersectionPoint2.x;
    double y2 = intersectionPoint2.y;
    double z2 = intersectionPoint2.z;

    //* checking if the t1 and t2 is outside the range

    bool isoutsideT1 = length > 0 && isOutside(x1, x_min, x_max) ||
                       width > 0 && isOutside(y1, y_min, y_max) ||
                       height > 0 && isOutside(z1, z_min, z_max);

    bool isoutsideT2 = length > 0 && isOutside(x2, x_min, x_max) ||
                       width > 0 && isOutside(y2, y_min, y_max) ||
                       height > 0 && isOutside(z2, z_min, z_max);

    if (isoutsideT1 && isoutsideT2)
      return -1;
    else if (isoutsideT1)
      return t2;
    else if (isoutsideT2)
      return t1;
    else
      return min(t1, t2);
  }

  double intersect(Ray *r, point *Color, int level)
  {

    double t = getIntersecting(*r);

    if (t <= 0)
      return -1;

    if (level == 0)
      return t;

    Fill_Color_Obj(this, r, Color, level, t);

    return t;
  }
};