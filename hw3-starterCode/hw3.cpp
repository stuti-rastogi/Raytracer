/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Stuti Rastogi
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

using namespace std;

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

#define PI 3.1415926535

// constants defined to access different planes
#define XY 100
#define YZ 101
#define ZX 102

// constants defined to see which type of object dealing with
#define SPHERE 31
#define TRIANGLE 32

unsigned char buffer[HEIGHT][WIDTH][3];


struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
  // added reflectivity co-effecient for one triangle, not per vertex
  double reflectivity[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
  // added reflectivity co-effecient for the sphere
  double reflectivity[3];
};

struct Light
{
  double position[3];
  double color[3];
};

// 3D vector - use for point and direction
struct Vector3
{
	double x;
	double y;
	double z;
};

// definition of any ray
struct Ray
{
	Vector3 origin;			// point
	Vector3 direction;		// vector
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

// the co-efficient used to calculate final color with reflection
double refl_coeff[3];


void plot_pixel_display(int x, int y, unsigned char r, 
						unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, 
						unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, 
						unsigned char g, unsigned char b);

//Vector Helper Functions

// adds two vectors
Vector3 add(Vector3 a, Vector3 b)
{
	Vector3 result;
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;
	return result;
}

// subtracts one vector from another
Vector3 subtract(Vector3 a, Vector3 b)
{
	Vector3 result;
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;
	return result;
}

// multiply vector with scalar
Vector3 multiplyWithScalar(double t, Vector3 a)
{
	Vector3 result;
	result.x = t * a.x;
	result.y = t * a.y;
	result.z = t * a.z;
	return result;
}

// get the length of any vector
double length(Vector3 a)
{
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

// normalise a vector
Vector3 normalise(Vector3 a)
{
	double magnitude = length(a);
	if (magnitude == 0)
		return a;
	Vector3 result;
	result.x = a.x / magnitude;
	result.y = a.y / magnitude;
	result.z = a.z / magnitude;
	return result;
}

// cross product of two vectors
Vector3 crossProduct(Vector3 a, Vector3 b)
{
	Vector3 result;
	result.x = a.y * b.z - a.z * b.y;
	result.y = a.z * b.x - a.x * b.z;
	result.z = a.x * b.y - a.y * b.x;
	return result;
}

// dot product of two vectors
double dotProduct(Vector3 a, Vector3 b)
{
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

// print a vector in a formatted way - used for debugging purposes
void printVector(Vector3 a)
{
	printf("(%f, %f, %f)\n", a.x, a.y, a.z);
}

// method that generates a ray originating from the camera (eye point)
// passing through the point (x, y)
Ray generateRay(int x, int y)
{
	Ray r;

	// camera at (0, 0, 0)
	r.origin.x = 0.0;
	r.origin.y = 0.0;
	r.origin.z = 0.0;

	Vector3 toPoint;		// pixel in image plane through which ray passes

	double a = (1.0 * WIDTH) / HEIGHT;
	double tanTerm = tan((fov / 2.0) * PI / 180.0);

	// slide 15-10
	toPoint.x = (-a * tanTerm) + (2.0 * a * tanTerm * (1.0 * x / WIDTH));
	toPoint.y = (-tanTerm) + (2.0 * tanTerm * (1.0 * y / HEIGHT));
	toPoint.z = -1.0;		// image plane at -1.0

	r.direction = normalise(subtract(toPoint, r.origin));
	return r;
}

// method that deals with the intersection of a ray with sphere
// returns true if intersection, else false
// the parameter t is assigned to point of closest intersection
bool raySphereIntersection (Ray r, Sphere s, double *t)
{
	Vector3 sphereOrigin = {s.position[0], s.position[1], s.position[2]};
	// (x0 - xc, y0 - yc, z0 - zc)
	Vector3 differenceVector = subtract(r.origin, sphereOrigin);
	
	// slide 16-06
	// quadratic equation constants
	// a = 1, hence not needed
	double b = 2 * (r.direction.x * differenceVector.x +
			r.direction.y * differenceVector.y +
			r.direction.z * differenceVector.z);
	double c = pow(differenceVector.x, 2) + pow(differenceVector.y, 2) + 
				pow(differenceVector.z, 2) - pow(s.radius, 2);

	double discriminant = (b * b) - (4 * c);

	// negative discriminant
	if (discriminant < 1e-10)
	{
		return false;
	}
	else
	{
		double t0 = (-b + sqrt(discriminant))/2.0;
		double t1 = (-b - sqrt(discriminant))/2.0;

		// if both negative, no intersection
		// if only one positive, that is point of intersection
		// if both positive, find closer point by minimum t
		if (t0 < 1e-10)
		{
			if (t1 < 1e-10)
			{
				return false;
			}
			else
			{
				*t = t1;
				return true;
			}
		}
		else
		{
			if (t1 < 1e-10)
			{
				*t = t0;
				return true;
			}
			else
			{
				*t = min(t0, t1);
				return true;
			}
		}

	}
	return false;
}


// method to calculate the Barycentric co-ordinates of a point in a triangle
// the values of alpha, beta, gamma are assigned
// point - the point at which Barycentric co-rodinates are needed
// triangleNormal - needed to find orientation of triangle
void findBarycentricCoordinates(Vector3 point, Triangle t, 
								Vector3 triangleNormal, 
								double *alpha, double *beta, double *gamma)
{
	// check which plane to project on
	Vector3 normalToXY = {0.0, 0.0, 1.0};
	Vector3 normalToYZ = {1.0, 0.0, 0.0};
	Vector3 normalToZX = {0.0, 1.0, 0.0};

	// take absolute value because triangle normal could be either side
	double XY_dotProduct = fabs(dotProduct(triangleNormal, normalToXY));
	double YZ_dotProduct = fabs(dotProduct(triangleNormal, normalToYZ));
	double ZX_dotProduct = fabs(dotProduct(triangleNormal, normalToZX));

	int planeToProject;

	// take the plane with which angle is minimum (maximum dot product)
	if (XY_dotProduct > ZX_dotProduct)
	{
		if (XY_dotProduct > YZ_dotProduct)
			planeToProject = XY;
		else
			planeToProject = YZ;
	}
	else
	{
		if (ZX_dotProduct > YZ_dotProduct)
			planeToProject = ZX;
		else
			planeToProject = YZ;
	}

	// calculate 2D projected vertices and intersection point
	double c0_x, c0_y, c1_x, c1_y, c2_x, c2_y;
	double c_x, c_y;
	switch(planeToProject)
	{
		case XY:
			c0_x = t.v[0].position[0];
			c0_y = t.v[0].position[1];

			c1_x = t.v[1].position[0];
			c1_y = t.v[1].position[1];

			c2_x = t.v[2].position[0];
			c2_y = t.v[2].position[1];

			c_x = point.x;
			c_y = point.y;
			break;

		case YZ:
			c0_x = t.v[0].position[1];
			c0_y = t.v[0].position[2];

			c1_x = t.v[1].position[1];
			c1_y = t.v[1].position[2];

			c2_x = t.v[2].position[1];
			c2_y = t.v[2].position[2];

			c_x = point.y;
			c_y = point.z;
			break;

		case ZX:
			c0_x = t.v[0].position[2];
			c0_y = t.v[0].position[0];

			c1_x = t.v[1].position[2];
			c1_y = t.v[1].position[0];

			c2_x = t.v[2].position[2];
			c2_y = t.v[2].position[0];

			c_x = point.z;
			c_y = point.x;
			break;
	}

	// calculate the four areas
	double area_c_c1_c2, area_c0_c_c2, area_c0_c1_c, area_c0_c1_c2;

	// slide 16-20
	// A(ABC) = (1/2) ((bx – ax)(cy – ay) – (cx – ax) (by – ay))
	area_c0_c1_c2 = (1.0/2.0) * ((c1_x - c0_x) * (c2_y - c0_y) - 		
								(c2_x - c0_x) * (c1_y - c0_y));
	area_c_c1_c2 = (1.0/2.0) * ((c1_x - c_x) * (c2_y - c_y) - 
								(c2_x - c_x) * (c1_y - c_y));
	area_c0_c_c2 = (1.0/2.0) * ((c_x - c0_x) * (c2_y - c0_y) - 
								(c2_x - c0_x) * (c_y - c0_y));
	area_c0_c1_c = (1.0/2.0) * ((c1_x - c0_x) * (c_y - c0_y) - 
								(c_x - c0_x) * (c1_y - c0_y));

	// slide 16-18
	// compute barycentric co-ordinates
	*alpha = area_c_c1_c2 / area_c0_c1_c2;
	*beta = area_c0_c_c2 / area_c0_c1_c2;
	*gamma = area_c0_c1_c / area_c0_c1_c2;
}

// method that finds the normal of the given triangle
Vector3 findTriangleNormal(Triangle triangle)
{
	// assign the 3 vertices to Vector3 variables
	// find 2 edge vectors using the difference between vertices
	// compute cross product and normalise

	// vertex vectors
	Vector3 v0 = {triangle.v[0].position[0], triangle.v[0].position[1], 
												triangle.v[0].position[2]};
	Vector3 v1 = {triangle.v[1].position[0], triangle.v[1].position[1], 
												triangle.v[1].position[2]};
	Vector3 v2 = {triangle.v[2].position[0], triangle.v[2].position[1], 
												triangle.v[2].position[2]};

	// edge vectors
	Vector3 e01 = subtract(v1, v0);
	Vector3 e02 = subtract(v2, v0);

	Vector3 triangleNormal = normalise(crossProduct(e01, e02));
	return triangleNormal;
}

// method that computes the intersection of a ray with a triangle
// returns true if intersection, false otherwise
// t is the parameter of the ray at intersection point
bool rayTriangleIntersection(Ray r, Triangle triangle, double *t)
{
	// first find intersection of ray with the plane of the triangle
	// if it intersects the plane, check if point of intersection is inside
	// the triangle using Barycentric co-ordinates

	// plane of triangle
	Vector3 planeNormal = findTriangleNormal(triangle);

	//plane equation co-efficients
	// ax + by + cz + d = 0
	double a = planeNormal.x;
	double b = planeNormal.y;
	double c = planeNormal.z;

	// d = -a(v0_x) - b(v0_y) - c(v0_z) since v0 lies on the plane
	double d = -1.0 * (a * triangle.v[0].position[0] + 
						b * triangle.v[0].position[1] + 
						c * triangle.v[0].position[2]);

	// find intersection point with plane
	// slide 16-10
	double denominator = (a * r.direction.x + b * r.direction.y + 
												c * r.direction.z);

	// ray parallel to plane
	if (denominator == 0)
		return false;
	else
	{
		double t_intersection = -1.0 * (a * r.origin.x + b * r.origin.y + 
											c * r.origin.z + d) / denominator;
		
		// ray away from the plane, t < 0
		if (t_intersection <= 1e-10)
			return false;

		// check if intersection point inside triangle
		Vector3 intersectionPoint = add(r.origin, 
							multiplyWithScalar(t_intersection, r.direction));
		
		double alpha, beta, gamma;
		findBarycentricCoordinates(intersectionPoint, triangle, planeNormal, 
									&alpha, &beta, &gamma);

		// point inside triangle iff 0 <= alpha, beta, gamma <= 1
		if (alpha >= 0 && alpha <= 1 && 
				beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1)
		{
			*t = t_intersection;
			return true;
		}
		return false;
	}
}

// method that computes the color by the Phong lighting model at a point
// N - normal at the point
// objectHit - sphere or triangle, needed to use appropriate ks, kd, shininess
// objectIndex - access the correct sphere or triangle from the array
// red, green, blue - values set by the method
void computePhongColor(Light light, Vector3 point, Vector3 N, 
						int objectHit, int objectIndex,
						double *red, double *green, double *blue)
{
	double colors[3];
	Vector3 light_pos = {light.position[0], 
							light.position[1], light.position[2]};
	// vector from point to light
	Vector3 L = normalise(subtract(light_pos, point));

	Vector3 eye = {0.0, 0.0, 0.0};
	// view vector
	Vector3 V = normalise(subtract(eye, point));

	// R = 2(N.L)N - L
	Vector3 R = normalise(subtract(multiplyWithScalar(
									(2 * dotProduct(L, N)), N), L));

	double kd[3], ks[3], shi;			// at the point, the one to be used
	if (objectHit == SPHERE)
	{
		Sphere s = spheres[objectIndex];

		for (int i = 0; i < 3; i++)
		{
			kd[i] = s.color_diffuse[i];
			ks[i] = s.color_specular[i];
		}

		shi = s.shininess;
	}
	else if (objectHit == TRIANGLE)
	{
		// we need to interpolate ks, kd and shininess for the triangle
		Triangle t = triangles[objectIndex];
		double alpha, beta, gamma;
		findBarycentricCoordinates(point, t, findTriangleNormal(t), 
									&alpha, &beta, &gamma);

		// for each channel - R, G, B
		for (int i = 0; i < 3; i++)
		{
			kd[i] = (alpha * t.v[0].color_diffuse[i]) + 
					(beta * t.v[1].color_diffuse[i]) + 
					(gamma * t.v[2].color_diffuse[i]);
			ks[i] = (alpha * t.v[0].color_specular[i]) + 
					(beta * t.v[1].color_specular[i]) + 
					(gamma * t.v[2].color_specular[i]);
		}

		shi = (alpha * t.v[0].shininess) + 
			(beta * t.v[1].shininess) + 
			(gamma * t.v[2].shininess);
	}

	// negative LdotN, RdotV => clamp to 0
	double LdotN = dotProduct(L, N);
	if (LdotN <= 1e-10)
		LdotN = 0.0;
	double RdotV = dotProduct(R, V);
	if (RdotV <= 1e-10)
		RdotV = 0.0;

	// R, G, B
	for (int i = 0; i < 3; i++)
	{
		// Phong equation
		colors[i] = light.color[i] * ((kd[i] * LdotN) + 
					(ks[i] * pow(RdotV, shi)));
	}

	*red = colors[0];
	*green = colors[1];
	*blue = colors[2];
}

// method that looks for intersection by the ray with the scene and computes
// the color, by shadows and the Phong model
// objectHit and objectIndex are returned as they are needed for correct
// reflectivity values
bool colorComputedByRay(Ray primaryRay, double *redValue, double *greenValue, 
					double *blueValue, Vector3 *point, Vector3 *normal, 
					int *objHitReturn, int *objIndReturn)
{
	double r = 0.0, g = 0.0, b = 0.0;
	bool hit;
	Vector3 intersectionPoint;
	Vector3 intersectionNormal;

	// store 't's of intersection with objects to find minimum
	double triangleIntersections[MAX_TRIANGLES];
	double sphereIntersections[MAX_SPHERES];
	double t_fetch;

	// if no intersection - set t to some high value
	for (int sphereIndex = 0; sphereIndex < num_spheres; sphereIndex++)
	{			
		if (raySphereIntersection(primaryRay, spheres[sphereIndex], &t_fetch))
			sphereIntersections[sphereIndex] = t_fetch;
		else
			sphereIntersections[sphereIndex] = 1e30;
	}

	for (int triangleIndex = 0; triangleIndex < num_triangles; triangleIndex++)
	{
		if (rayTriangleIntersection(primaryRay, 
					triangles[triangleIndex], &t_fetch))
			triangleIntersections[triangleIndex] = t_fetch;
		else
			triangleIntersections[triangleIndex] = 1e30;
	}

	// check first point of intersection - minimum
	double t_intersection = 1e25;
	int objectHit;			// sphere or triangle
	int objectIndex;
	for (int sphereIndex = 0; sphereIndex < num_spheres; sphereIndex++)
	{
		if (sphereIntersections[sphereIndex] < t_intersection)
		{
			objectHit = SPHERE;
			objectIndex = sphereIndex;
			t_intersection = sphereIntersections[sphereIndex];
		}
	}
	for (int triangleIndex = 0; triangleIndex < num_triangles; triangleIndex++)
	{
		if (triangleIntersections[triangleIndex] < t_intersection)
		{
			objectHit = TRIANGLE;
			objectIndex = triangleIndex;
			t_intersection = triangleIntersections[triangleIndex];
		}
	}

	// no intersection with anything - background color
	if (t_intersection > 1e20)
	{
		r = 1.0;
		g = 1.0;
		b = 1.0;
		hit = false;
	}
	else
	{
		hit = true;
		intersectionPoint = add(primaryRay.origin, 
			multiplyWithScalar(t_intersection, primaryRay.direction));

		// launch shadow ray from hit point for each light
		for (int lightIndex = 0; lightIndex < num_lights; lightIndex++)
		{
			bool inShadow = false;
			Ray shadowRay;
	
			// intersectionPoint is origin
			// toPoint is light position
			shadowRay.origin = intersectionPoint;
			Vector3 toPoint;

			toPoint.x = lights[lightIndex].position[0];
			toPoint.y = lights[lightIndex].position[1];
			toPoint.z = lights[lightIndex].position[2];
			
			shadowRay.direction = normalise(
								subtract(toPoint, intersectionPoint));
			double t_max = length(subtract(toPoint, intersectionPoint));

			// check if ray intersects with anything => shadow 
			double t;

			// intersection should be found before light, hence t < tmax
			// and t should be positive
			for (int sphereIndex = 0; sphereIndex < num_spheres; sphereIndex++)
			{
				if (raySphereIntersection(shadowRay, spheres[sphereIndex], &t) 
										&& t > 1e-4 && t < t_max)
				{
					inShadow = true;
					break;
				}
			}

			// keep going only if not in shadow yet, otherwise avoid this
			// check with every triangle
			if (inShadow == false)
			{
				for (int triangleIndex = 0; triangleIndex < num_triangles; 
												triangleIndex++)
				{
					if (rayTriangleIntersection(shadowRay, 
						triangles[triangleIndex], &t) && t > 1e-4 && t < t_max)
					{
						inShadow = true;
						break;
					}
				}
			}

			// if point in shadow, contribution for this light is 0, so don't
			// do anything
			if (inShadow == false)
			{
				// compute color for this light
				// first compute the normal depending on kind of object
				double redComputed = 0.0, greenComputed = 0.0, 
											blueComputed = 0.0;
				if (objectHit == SPHERE)
				{
					Sphere sphereHit = spheres[objectIndex];
					Vector3 sphereHitOrigin = {sphereHit.position[0], 
											sphereHit.position[1], 
											sphereHit.position[2]};
					// slide 16-07
					intersectionNormal = normalise(multiplyWithScalar(
						(1.0/sphereHit.radius), subtract(intersectionPoint, 
														sphereHitOrigin)));
				}
				else if (objectHit == TRIANGLE)
				{
					Triangle triangleHit = triangles[objectIndex];
					double alpha, beta, gamma;
					Vector3 normal = findTriangleNormal(triangleHit);
					findBarycentricCoordinates(intersectionPoint, triangleHit, 
						normal, &alpha, &beta, &gamma);

					// make vertex normal vectors
					Vector3 n0 = {triangleHit.v[0].normal[0], 
						triangleHit.v[0].normal[1], triangleHit.v[0].normal[2]};
					Vector3 n1 = {triangleHit.v[1].normal[0], 
						triangleHit.v[1].normal[1], triangleHit.v[1].normal[2]};
					Vector3 n2 = {triangleHit.v[2].normal[0], 
						triangleHit.v[2].normal[1], triangleHit.v[2].normal[2]};

					// Barycentric interpolation
					Vector3 temp1 = multiplyWithScalar(alpha, n0);
					Vector3 temp2 = multiplyWithScalar(beta, n1);
					Vector3 temp3 = multiplyWithScalar(gamma, n2);

					intersectionNormal = normalise(add(temp3, 
												add(temp1, temp2)));
				}

				computePhongColor(lights[lightIndex], intersectionPoint, 
					intersectionNormal, objectHit, objectIndex, 
					&redComputed, &greenComputed, &blueComputed);
				
				// accumulate colors
				r = r + redComputed;
				g = g + greenComputed;
				b = b + blueComputed;
			}
		}
	}

	// set the values to be returned
	*point = intersectionPoint;
	*normal = intersectionNormal;
	*objHitReturn = objectHit;
	*objIndReturn = objectIndex;
	*redValue = r;
	*greenValue = g;
	*blueValue = b;

	return hit;
}

// the method that does ray tracing for each pixel
void draw_scene()
{
	for (unsigned int x = 0; x < WIDTH; x++)
	{
	    glPointSize(2.0);  
	    glBegin(GL_POINTS);
		for (unsigned int y = 0; y < HEIGHT; y++)
		{
			double red = 0, green = 0, blue = 0;

			// add ambient light to each pixel
			red = ambient_light[0];
			green = ambient_light[1];
			blue = ambient_light[2];

			// get ray through this pixel
			Ray primaryRay = generateRay(x, y);

			double redPhong = 0.0, greenPhong = 0.0, 
					bluePhong = 0.0;
			double redReflected = 0.0, greenReflected = 0.0, 
					blueReflected = 0.0;
			
			Vector3 P;				// intersection point
			Vector3 N;				// intersection normal
			int objectHit, objectIndex;
			bool hit = colorComputedByRay(primaryRay, &redPhong, &greenPhong, 
				&bluePhong, &P, &N, &objectHit, &objectIndex);
			
			// if not hitting background, do recursion
			if (hit)
			{
				// fetch corresponding reflection co-efficients
				if (objectHit == SPHERE)
				{
					refl_coeff[0] = spheres[objectIndex].reflectivity[0];
					refl_coeff[1] = spheres[objectIndex].reflectivity[1];
					refl_coeff[2] = spheres[objectIndex].reflectivity[2];
				}
				else
				{
					refl_coeff[0] = triangles[objectIndex].reflectivity[0];
					refl_coeff[1] = triangles[objectIndex].reflectivity[1];
					refl_coeff[2] = triangles[objectIndex].reflectivity[2];
				}

				// if reflectivity is 0, no need to send secondary ray
				if (refl_coeff[0] != 0 && refl_coeff[1] != 0 && 
					refl_coeff[2] != 0)
				{
					// compute reflected ray
					Vector3 L = normalise(subtract(primaryRay.origin, P));
					
					// R = 2(N.L)N - L
					Vector3 R = normalise(subtract(multiplyWithScalar(
									(2 * dotProduct(L, N)), N), L));
	
					Ray reflectedRay;
					reflectedRay.origin = P;
					reflectedRay.direction = R;
					
					colorComputedByRay(reflectedRay, &redReflected, 
						&greenReflected, &blueReflected, &P, &N, 
						&objectHit, &objectIndex);
				}

				// formula given on HW page
				red = red + (1-refl_coeff[0])*redPhong + 
							refl_coeff[0]*redReflected;
				green = green + (1-refl_coeff[1])*greenPhong + 
							refl_coeff[1]*greenReflected;
				blue = blue + (1-refl_coeff[2])*bluePhong + 
							refl_coeff[2]*blueReflected;
			}
			else
			{
				red = red + redPhong;
				green = green + greenPhong;
				blue = blue + bluePhong;				
			}


			// clamp values if exceed 1
			if (red >= 1.0)
				red = 1.0;
			if (green >= 1.0)
				green = 1.0;
			if (blue >= 1.0)
				blue = 1.0;

			// convert float to char to plot
			plot_pixel(x, y, (unsigned char)ceil(red*255.0), 
								(unsigned char)ceil(green*255.0), 
								(unsigned char)ceil(blue*255.0));
		}

    	glEnd();
    	glFlush();
  	}
  	printf("Done!\n"); 
  	fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, 
					unsigned char g, unsigned char b)
{
	glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
	glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, 
					unsigned char g, unsigned char b)
{
	buffer[y][x][0] = r;
	buffer[y][x][1] = g;
	buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
	plot_pixel_display(x,y,r,g,b);
	if(mode == MODE_JPEG)
		plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
	printf("Saving JPEG file: %s\n", filename);

	ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
	if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
		printf("Error in Saving\n");
	else 
		printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
	if(strcasecmp(expected,found))
	{
		printf("Expected '%s ' found '%s '\n", expected, found);
		printf("Parse error, abnormal abortion\n");
		exit(0);
	}
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check(check,str);
	fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
	printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check("rad:",str);
	fscanf(file,"%lf",r);
	printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
	char s[100];
	fscanf(file,"%s",s);
	parse_check("shi:",s);
	fscanf(file,"%lf",shi);
	printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
	FILE * file = fopen(argv,"r");
	int number_of_objects;
	char type[50];
	Triangle t;
	Sphere s;
	Light l;
	fscanf(file,"%i", &number_of_objects);

	printf("number of objects: %i\n",number_of_objects);

	parse_doubles(file,"amb:",ambient_light);

	for(int i=0; i<number_of_objects; i++)
	{
		fscanf(file,"%s\n",type);
		printf("%s\n",type);
		if(strcasecmp(type,"triangle")==0)
		{
			printf("found triangle\n");
			for(int j=0;j < 3;j++)
			{
		        parse_doubles(file,"pos:",t.v[j].position);
		        parse_doubles(file,"nor:",t.v[j].normal);
		        parse_doubles(file,"dif:",t.v[j].color_diffuse);
		        parse_doubles(file,"spe:",t.v[j].color_specular);
		        parse_shi(file,&t.v[j].shininess);
      		}

      		if(num_triangles == MAX_TRIANGLES)
      		{
		    	printf("too many triangles, increase MAX_TRIANGLES!\n");
		        exit(0);
      		}
      		triangles[num_triangles++] = t;
    	}

		else if(strcasecmp(type,"sphere")==0)
		{
			printf("found sphere\n");

			parse_doubles(file,"pos:",s.position);
			parse_rad(file,&s.radius);
			parse_doubles(file,"dif:",s.color_diffuse);
			parse_doubles(file,"spe:",s.color_specular);
			parse_shi(file,&s.shininess);

			if(num_spheres == MAX_SPHERES)
			{
				printf("too many spheres, increase MAX_SPHERES!\n");
				exit(0);
      		}
      		spheres[num_spheres++] = s;
    	}

		else if(strcasecmp(type,"light")==0)
		{
			printf("found light\n");
			parse_doubles(file,"pos:",l.position);
			parse_doubles(file,"col:",l.color);

			if(num_lights == MAX_LIGHTS)
			{
				printf("too many lights, increase MAX_LIGHTS!\n");
				exit(0);
      		}
      		lights[num_lights++] = l;
    	}
    
    	else
    	{
			printf("unknown type in scene description:\n%s\n",type);
			exit(0);
    	}
  	}

  	// set reflectivity values for each object
  	// hard coding since not allowed to change scene descriptions
  	for (int i = 0; i < num_spheres; i++)
  	{
  		spheres[i].reflectivity[0] = 0.2;
  		spheres[i].reflectivity[1] = 0.2;
  		spheres[i].reflectivity[2] = 0.2;
  	} 

  	for (int i = 0; i < num_triangles; i++)
  	{
  		triangles[i].reflectivity[0] = 0.0;
		triangles[i].reflectivity[1] = 0.0;
  		triangles[i].reflectivity[2] = 0.0;
  	}

  	return 0;
}

void display()
{
}

void init()
{
	glMatrixMode(GL_PROJECTION);
	glOrtho(0,WIDTH,0,HEIGHT,1,-1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);
}

// this is the method for keyboard callback
void keyboardFunc(unsigned char key, int x, int y)
{
	switch (key)
	{
	  	// exit on ESC key
	    case 27:
			exit(0); // exit the program
	    break;
  	}
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
	if ((argc < 2) || (argc > 3))
	{  
		printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
		exit(0);
	}
	
	if(argc == 3)
	{
		mode = MODE_JPEG;
		filename = argv[2];
	}
	else if(argc == 2)
		mode = MODE_DISPLAY;

	glutInit(&argc,argv);
	loadScene(argv[1]);

	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowPosition(0,0);
	glutInitWindowSize(WIDTH,HEIGHT);
	int window = glutCreateWindow("Ray Tracer");
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboardFunc);
	init();
	glutMainLoop();
}