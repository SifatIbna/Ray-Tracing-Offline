#include "1605112_Header.h"
#include "bitmap_image.hpp"

#define pi (acos(-1.0))

#define eps 0.0000001

#define viewAngle 90
#define windowHeight 500
#define windowWidth 500

int imageWidth;
int imageHeight;

int N;
int nLights;

point eye(0, -100, 50), u(0, 0, 1), r(1, 0, 0), l(0, 1, 0);

void capture()
{
  cout << "Capturing" << endl;

  bitmap_image image(imageWidth, imageHeight);

  double planeDistance = (windowHeight / 2.0) / tan((viewAngle / 2.0) * (pi / 180.0));

  point topleft = eye + l * planeDistance - r * (windowWidth / 2.0) + u * (windowHeight / 2.0);

  double du = (1.0 * windowWidth) / imageWidth;
  double dv = (1.0 * windowHeight) / imageHeight;

  topleft = topleft + r * (0.5 * du) - u * (0.5 * dv);

  for (int i = 0; i < imageHeight; i++)
  {
    for (int j = 0; j < imageWidth; j++)
    {
      point curPixel = topleft + r * j * du - u * i * dv;

      point diff = curPixel - eye;
      diff = Normalized(diff);

      Ray ray(eye, diff);

      int nearest = -1;

      point dummyColor;
      double t_min = 999999;
      for (int k = 0; k < objects.size(); k++)
      {
        double t = objects[k]->intersect(&ray, &dummyColor, 0);

        if (t <= 0)
          continue;
        else if (t < t_min)
        {
          t_min = min(t, t_min);
          nearest = k;
        }
      }

      if (nearest != -1)
      {
        double t = objects[nearest]->intersect(&ray, &dummyColor, 1);

        image.set_pixel(j, i, dummyColor.x * 255, dummyColor.y * 255, dummyColor.z * 255);
      }
    }
  }

  image.save_image("image.bmp");
  cout << "Image captured" << endl;
}

void keyboardListener(unsigned char key, int x, int y)
{
  switch (key)
  {
  case '0':
    capture();
    break;

  case '1':
    Rotate(r, l, pi / 60);
    Multiplication(l, u, r);

    break;
  case '2':
    Rotate(r, l, -pi / 60);

    Multiplication(l, u, r);
    break;
  case '3':
    Rotate(l, u, pi / 60);

    Multiplication(u, r, l);

    break;
  case '4':
    Rotate(l, u, -pi / 60);

    Multiplication(u, r, l);

    break;
  case '5':
    Rotate(u, r, pi / 60);

    Multiplication(r, l, u);

    break;
  case '6':
    Rotate(u, r, -pi / 60);

    Multiplication(r, l, u);
    break;

  default:
    break;
  }
}

void specialKeyListener(int key, int x, int y)
{
  switch (key)
  {
  case GLUT_KEY_DOWN: //down arrow key
    eye.x -= l.x * 2;
    eye.y -= l.y * 2;
    eye.z -= l.z * 2;
    break;
  case GLUT_KEY_UP: // up arrow key
    eye.x += l.x * 2;
    eye.y += l.y * 2;
    eye.z += l.z * 2;
    break;

  case GLUT_KEY_RIGHT:
    eye.x += r.x * 2;
    eye.y += r.y * 2;
    eye.z += r.z * 2;
    break;
  case GLUT_KEY_LEFT:
    eye.x -= r.x * 2;
    eye.y -= r.y * 2;
    eye.z -= r.z * 2;
    break;

  case GLUT_KEY_PAGE_UP:
    eye.x += u.x * 2;
    eye.y += u.y * 2;
    eye.z += u.z * 2;
    break;
  case GLUT_KEY_PAGE_DOWN:
    eye.x -= u.x * 2;
    eye.y -= u.y * 2;
    eye.z -= u.z * 2;
    break;

  case GLUT_KEY_INSERT:
    break;

  case GLUT_KEY_HOME:

    break;
  case GLUT_KEY_END:

    break;

  default:
    break;
  }
}

void mouseListener(int button, int state, int x, int y)
{
  switch (button)
  {
  case GLUT_LEFT_BUTTON:
    break;

  case GLUT_RIGHT_BUTTON:
    //........
    break;

  case GLUT_MIDDLE_BUTTON:
    //........
    break;

  default:
    break;
  }
}

void display()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);

  glLoadIdentity();

  gluLookAt(eye.x, eye.y, eye.z, eye.x + l.x, eye.y + l.y, eye.z + l.z, u.x, u.y, u.z);

  glMatrixMode(GL_MODELVIEW);

  for (int i = 0; i < objects.size(); i++)
  {
    objects[i]->draw();
  }

  // for (int i = 0; i < lights.size(); i++)
  // {
  //   glBegin(GL_POINTS);
  //   glColor3f(255, 255, 255);
  //   glVertex3f(lights[i].pos.x, lights[i].pos.y, lights[i].pos.z);
  //   glEnd();
  // }

  glutSwapBuffers();
}

void animate()
{
  glutPostRedisplay();
}

void init()
{

  glClearColor(0, 0, 0, 0);

  glMatrixMode(GL_PROJECTION);

  glLoadIdentity();

  gluPerspective(90, 1, 1, 1000.0);
}

void loadData()
{
  ifstream fin("scene_test.txt");

  fin >> recurLevel;
  fin >> imageWidth;
  imageHeight = imageWidth;

  fin >> N; //* object count
  Object *obj;
  for (int i = 0; i < N; i++)
  {
    string objectName; //* object name
    fin >> objectName;
    if (objectName.compare("sphere") == 0)
    {
      point center;
      double radius;
      double r, g, b;
      double amb, diff, spec, recur;
      int shine;

      fin >> center.x >> center.y >> center.z;
      fin >> radius;
      obj = new Sphere(center, radius);

      fin >> r >> g >> b; //*getting the colors
      obj->setColor(r, g, b);

      fin >> amb >> diff >> spec >> recur; //* getting the co-efficients
      obj->setCoEfficients(amb, diff, spec, recur);

      fin >> shine; //* getting the shine value
      obj->setShine(shine);

      objects.push_back(obj);
    }

    else if (objectName.compare("triangle") == 0)
    {
      point p1, p2, p3;
      double r, g, b;
      double amb, diff, spec, recur;
      int shine;

      fin >> p1.x >> p1.y >> p1.z; //* triangle three points
      fin >> p2.x >> p2.y >> p2.z;
      fin >> p3.x >> p3.y >> p3.z;

      obj = new Triangle(p1, p2, p3);

      fin >> r >> g >> b; //* color parameters
      obj->setColor(r, g, b);

      fin >> amb >> diff >> spec >> recur; //* setting the co-efficients
      obj->setCoEfficients(amb, diff, spec, recur);

      fin >> shine;
      obj->setShine(shine);

      objects.push_back(obj);
    }
    else if (objectName.compare("general") == 0)
    {

      double coeff[10];
      double x, y, z;
      double length, width, height;
      double red, green, blue;
      double ambient, diffuse, specular, reflection;
      double shine;

      point base_point;

      for (int c = 0; c < 10; c++)
      {
        fin >> coeff[c];
      }

      fin >> x >> y >> z;
      base_point = point(x, y, z);

      fin >> length >> width >> height;
      General *obj = new General(coeff, base_point, length, width, height);

      fin >> red >> green >> blue;
      obj->setColor(red, green, blue);

      fin >> ambient >> diffuse >> specular >> reflection;
      obj->setCoEfficients(ambient, diffuse, specular, reflection);

      fin >> shine;
      obj->setShine(shine);

      objects.push_back(obj);
    }
  }

  fin >> nLights;
  // cout << "Light Count: " << nLights;
  point lightPos;
  point lightColor;
  for (int i = 0; i < nLights; i++)
  {
    fin >> lightPos.x >> lightPos.y >> lightPos.z;
    fin >> lightColor.x >> lightColor.y >> lightColor.z;
    Light light(lightPos, lightColor);
    lights.push_back(light);
  }

  obj = new Floor(1000, 30);
  obj->setCoEfficients(0.5, 0.4, 0.3, 0.1);
  obj->setShine(5);
  objects.push_back(obj);
}

void freeMemory()
{
  vector<Light>().swap(lights);
  vector<Object *>().swap(objects);
}

int main(int argc, char **argv)
{
  loadData();
  glutInit(&argc, argv);
  glutInitWindowSize(windowHeight, windowWidth);
  glutInitWindowPosition(0, 0);
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

  glutCreateWindow("Ray Tracing");

  init();

  glEnable(GL_DEPTH_TEST);

  glutDisplayFunc(display);
  glutIdleFunc(animate);

  glutKeyboardFunc(keyboardListener);
  glutSpecialFunc(specialKeyListener);
  glutMouseFunc(mouseListener);

  glutMainLoop();
  freeMemory();
  return 0;
}