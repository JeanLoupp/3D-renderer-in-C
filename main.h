#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <libgen.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <time.h>

typedef struct Point3D
{
    float x;
    float y;
    float z;
    SDL_FPoint point2D;
    float dz;
} Point3D;

typedef struct Vertex
{
    float x;
    float y;
    float z;
} Vertex;

typedef struct VertexN
{
    float x;
    float y;
    float z;
    float lightD;
    float lightS;
} VertexN;

typedef struct Face
{
    Point3D *p[3];
    VertexN *n[3];
    SDL_FPoint *t[3];
    float distance;
    short show;
    Vertex normal;
} Face;

typedef struct ColorFloat
{
    float r,g,b;
} ColorFloat;

typedef struct PhongModel
{
    ColorFloat ambiant;
    ColorFloat diffuse;
    ColorFloat specular;
} PhongModel;

typedef struct Polygon
{
    Face *faces;
    SDL_Surface *surfaceT;
    int size;
    char name[64];
    PhongModel phongModel;
} Polygon;

typedef struct Camera
{
    float x;
    float y;
    float z;
    float cx;
    float cy;
    float sx;
    float sy;
    float Rx;
    float Ry;
    Vertex center;
    float distance;
} Camera;

void projectPoints(Point3D points[], Camera *camera);
void projectPoint(Point3D *point, Camera *camera);
void createCube(int size, Point3D points[]);
void RenderPoints3D(SDL_Renderer *renderer, Point3D points[]);
void readFile(char *filePath, Point3D **points, SDL_FPoint **textures, VertexN **normals, Polygon **polygons, Camera *camera);
void readMTL(char file[], Polygon **polygons);
void calculateDistance(Face faces[], int size);
void calculateLightD(VertexN norms[], Vertex *light);
void calculateLightS(VertexN norms[], Vertex *light, Camera *camera);
void calculateNorms(Face faces[], int size);
void sortFaces(Face faces[], int n);
void RenderFaces(SDL_Renderer *renderer, Face faces[]);
void RenderPolygons(SDL_Renderer *renderer, Polygon polygons[]);
void RenderPolygonsTextured(SDL_Renderer *renderer, SDL_Surface *screenSurface, Polygon polygons[], float zbuffer[], Camera *camera);
float Q_rsqrt(float number);
void setCameraAngle(Camera *camera, float Rx, float Ry);
void addCameraAngle(Camera *camera, float Rx, float Ry);
void calculateCameraPos(Camera * camera);
void updateMin(Vertex *point, float x, float y, float z);
void updateMax(Vertex *point, float x, float y, float z);
float min(float a, float b);
float max(float a, float b);
void drawTriangle(Point3D *triangle[], VertexN *normals[], SDL_FPoint *texture[], Uint32 *screenPixels, SDL_Surface *Tsurface, PhongModel *phongModel, float zbuffer[], Camera *camera);
void setVisibleFaces(Face faces[], int size, Camera *camera);
Polygon* searchPoly(Polygon polygons[], char name[]);
void interpolatePoint(Point3D *p1, Point3D *p2, Point3D *dest, float alpha);
void interpolate2DPoint(SDL_FPoint *p1, SDL_FPoint *p2, SDL_FPoint *dest, float alpha);
#endif