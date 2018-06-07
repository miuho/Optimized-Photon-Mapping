/**
 * @file raytacer.hpp
 * @brief Raytracer class
 *
 * 
 *
 * @author HingOn Miu (hmiu)
 * 
 */

#ifndef _462_RAYTRACER_HPP_
#define _462_RAYTRACER_HPP_

#define MAX_DEPTH 3

#include "math/color.hpp"
#include "math/random462.hpp"

namespace _462 {

class Scene;
class Ray;
struct Intersection;
struct Photon;
struct Neighbor;
class Raytracer
{
public:

    Raytracer();

    ~Raytracer();

    bool initialize(Scene* scene, size_t num_samples,
                    size_t width, size_t height);

    bool raytrace(unsigned char* buffer, real_t* max_time);

private:

    Color3 trace_pixel(const Scene* scene,
		       size_t x,
		       size_t y,
		       size_t width,
		       size_t height,
                       Neighbor *neighbors);

    // the scene to trace
    Scene* scene;

    // the dimensions of the image to trace
    size_t width, height;

    // the next row to raytrace
    size_t current_row;

    unsigned int num_samples;
};

} /* _462 */

#endif /* _462_RAYTRACER_HPP_ */
