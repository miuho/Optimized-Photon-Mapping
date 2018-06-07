#include "ray.hpp"

namespace _462 {

Vector3 dir;
Vector3 up;
real_t AR;

Vector3 cR;
Vector3 cU;
real_t dist;
Vector3 pos;

Ray::Ray(){}

Ray::Ray(Vector3 e, Vector3 d)
{
    this->e = e;
    this->d = d;
}

void Ray::init(const Camera& camera)
{
    dir = camera.get_direction();
    up = camera.get_up();
    AR = camera.get_aspect_ratio();
    cR = cross(dir, up);
    cU = cross(cR, dir);
    pos = camera.get_position();
    dist = tan(camera.get_fov_radians()/2.0);
}

Vector3 Ray::get_pixel_dir(real_t ni, real_t nj)
{
    return normalize(dir + dist*(nj*cU + AR*ni*cR));
}

}
