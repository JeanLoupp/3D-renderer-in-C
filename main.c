#include "main.h"

#define SQUARE(x) ((x) * (x))

int WIDTH;
int HEIGHT;
float focal;
int totalPoints = 0;
int totalNorms = 0;
int totalTextures = 0;
int totalPolys = 0;
char hasTextures = 0;
char hasNorms = 0;
char uniformLight = 0;
float near = 1;
short orbit = 0;
SDL_PixelFormat *pixelFormat;
const float Ka = 1;
const float Kd = 1;
const float Ks = 1;

int main(int argc, char *argv[])
{
    if (argc < 6){
        perror("Donner comme arguments: FILEPATH WIDTH HEIGHT FOCAL DISTANCE");
        return 1;
    }

    /* Lecture des entrées */

    char *filepath = argv[1];
    WIDTH = atoi(argv[2]);
    HEIGHT = atoi(argv[3]);
    focal = atoi(argv[4]);
    float distance = atoi(argv[5]);

    /* Initialisations */
    SDL_Window *window = NULL;
    SDL_Renderer *renderer = NULL;
    Point3D *points = NULL;
    VertexN *normals = NULL;
    SDL_FPoint *textures = NULL;
    SDL_Texture *screenTexture = NULL;
    SDL_Surface *screenSurface = NULL;
    Polygon *polygons = NULL;
    float zbuffer[WIDTH*HEIGHT];
    pixelFormat = SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888);

    for (int i=0; i<WIDTH*HEIGHT; i++){
        zbuffer[i] = 100000;
    }

    Camera camera;
    camera.distance = distance;

    Vertex light = {1/SDL_sqrt(3), 1/SDL_sqrt(3), 1/SDL_sqrt(3)};

    readFile(filepath, &points, &textures, &normals, &polygons, &camera);

    SDL_Init(SDL_INIT_VIDEO);

    window = SDL_CreateWindow("SDL2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                              WIDTH, HEIGHT, SDL_WINDOW_SHOWN);

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    /* Debut du code */

    screenTexture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STATIC, WIDTH, HEIGHT);
    screenSurface = SDL_CreateRGBSurfaceWithFormat(0, WIDTH, HEIGHT, 32, SDL_PIXELFORMAT_RGBA8888);
    SDL_LockSurface(screenSurface);

    int isClicking=0;

    clock_t debut;
    SDL_Keymod keymod;
    Uint8 *clavier;

    SDL_Event event;
    SDL_bool quit = SDL_FALSE;
    while (!quit)
    {
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
                quit = SDL_TRUE;
            else if (event.type == SDL_MOUSEBUTTONDOWN)
                isClicking = 1;
            else if (event.type == SDL_MOUSEBUTTONUP)
                isClicking = 0;
            else if (event.type == SDL_MOUSEWHEEL)
            {
                if(event.wheel.y > 0) // scroll up
                {
                    camera.distance += -event.wheel.y*3;
                    if (camera.distance < 0)
                        camera.distance = 0;
                    addCameraAngle(&camera, 0, 0);
                }
                else if(event.wheel.y < 0) // scroll down
                {
                    camera.distance += -event.wheel.y*3;
                    addCameraAngle(&camera, 0, 0);
                }
            }
            else if (event.type == SDL_MOUSEMOTION){
                if (isClicking)
                    addCameraAngle(&camera, (float) event.motion.yrel/100, (float) event.motion.xrel/100);
            }
        }

        // Gestion évènements
        SDL_PumpEvents();
        clavier = SDL_GetKeyboardState(NULL);
        keymod = SDL_GetModState();

        if (!orbit){
            if (clavier[SDL_SCANCODE_W]){
                camera.z += 100*camera.cy;
                camera.x += 100*camera.sy;
            }if (clavier[SDL_SCANCODE_S]){
                camera.z += -100*camera.cy;
                camera.x += -100*camera.sy;
            }if (clavier[SDL_SCANCODE_D]){
                camera.z += -100*camera.sy;
                camera.x += 100*camera.cy;
            }if (clavier[SDL_SCANCODE_A]){
                camera.z += 100*camera.sy;
                camera.x += -100*camera.cy;
            }if (clavier[SDL_SCANCODE_SPACE])
                camera.y += 100;
            if (clavier[SDL_SCANCODE_LSHIFT])
                camera.y += -100;
        }

        calculateLightS(normals, &light, &camera);

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        debut = clock();
        projectPoints(points, &camera);

        //printf("Project %lf \n", ((double)(clock() - debut)) / CLOCKS_PER_SEC);
        for (int i=0; i<totalPolys; i++){
            //setVisibleFaces(polygons[i].faces, polygons[i].size, &camera);
        }
        debut = clock();

        RenderPolygonsTextured(renderer, screenSurface, polygons, zbuffer, &camera);
        SDL_UpdateTexture(screenTexture, NULL, screenSurface->pixels, screenSurface->pitch);
        SDL_RenderCopy(renderer, screenTexture, NULL, NULL);
        //printf("Render %lf\n", ((double)(clock() - debut)) / CLOCKS_PER_SEC);

        SDL_RenderPresent(renderer);

        //SDL_Delay(20);
    }

    SDL_DestroyWindow(window);
    SDL_DestroyRenderer(renderer);
    SDL_FreeSurface(screenSurface);
    SDL_DestroyTexture(screenTexture);
    SDL_FreeFormat(pixelFormat);
    SDL_Quit();
    free(points);
    free(normals);
    free(textures);
    for (int i=0; i<totalPolys; i++){
        free(polygons[i].faces);
    }
    free(polygons);
    return 0;
}

void projectPoints(Point3D points[], Camera *camera)
{
    float x, y, z, cx, cy, sx, sy, dx, dy, dz;
    cx = camera->cx;
    cy = camera->cy;
    sx = camera->sx;
    sy = camera->sy;

    for (int i = 0; i < totalPoints; i++)
    {
        x = points[i].x - camera->x;
        y = points[i].y - camera->y;
        z = points[i].z - camera->z;

        dx = cy*x - sy*z;
        dy = sx*(cy*z + sy*x) + cx*y;
        dz = cx*(cy*z + sy*x) - sx*y;

        points[i].point2D.x = focal * dx / dz + WIDTH / 2;
        points[i].point2D.y = HEIGHT/2 - focal * dy / dz;
        points[i].dz = (float) dz;
    }
}

void projectPoint(Point3D *point, Camera *camera)
{
    // Attention dz pas attribué
    float x, y, z, cx, cy, sx, sy, dx, dy, dz;
    cx = camera->cx;
    cy = camera->cy;
    sx = camera->sx;
    sy = camera->sy;

    x = point->x - camera->x;
    y = point->y - camera->y;
    z = point->z - camera->z;

    dx = cy*x - sy*z;
    dy = sx*(cy*z + sy*x) + cx*y;
    dz = cx*(cy*z + sy*x) - sx*y;

    point->point2D.x = focal * dx / dz + WIDTH / 2;
    point->point2D.y = HEIGHT/2 - focal * dy / dz;
}

void interpolatePoint(Point3D *p1, Point3D *p2, Point3D *dest, float alpha){
    dest->x = (1-alpha)*p1->x + alpha*p2->x;
    dest->y = (1-alpha)*p1->y + alpha*p2->y;
    dest->z = (1-alpha)*p1->z + alpha*p2->z;
}

void interpolateNormal(VertexN *n1, VertexN *n2, VertexN *dest, float alpha){
    dest->x = (1-alpha)*n1->x + alpha*n2->x;
    dest->y = (1-alpha)*n1->y + alpha*n2->y;
    dest->z = (1-alpha)*n1->z + alpha*n2->z;
    dest->lightD = (1-alpha)*n1->lightD + alpha*n2->lightD;
}

void interpolate2DPoint(SDL_FPoint *p1, SDL_FPoint *p2, SDL_FPoint *dest, float alpha){
    dest->x = (1-alpha)*p1->x + alpha*p2->x;
    dest->y = (1-alpha)*p1->y + alpha*p2->y;
}

void createCube(int size, Point3D points[])
{
    size = size / 2;
    int count = 0;
    for (int i = -1; i < 2; i += 2)
    {
        for (int j = -1; j < 2; j += 2)
        {
            for (int k = -1; k < 2; k += 2)
            {
                points[count].x = size * i;
                points[count].y = size * j;
                points[count].z = size * k;
                count++;
            }
        }
    }
}

void RenderPoints3D(SDL_Renderer *renderer, Point3D points[])
{
    for (int i = 0; i < totalPoints; i++)
    {
        SDL_RenderDrawPoint(renderer, points[i].point2D.x, points[i].point2D.x);
    }
}

void RenderPolygons(SDL_Renderer *renderer, Polygon polygons[]){
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_Color color = {255, 255, 255, 255};
    Face *faces;

    for (int j=0; j<totalPolys; j++){
        // TOFIX attention à l'ordre des polygons
        faces = polygons[j].faces;

        sortFaces(faces, polygons[j].size);

        SDL_Vertex vertex_1 = {{0, 0}, color, {1, 1}};
        SDL_Vertex vertex_2 = {{0, 0}, color, {1, 1}};
        SDL_Vertex vertex_3 = {{0, 0}, color, {1, 1}};

        SDL_Vertex vertices[3];

        for (int i = polygons[j].size - 1; i >= 0; i--)
        {
            vertex_1.position = faces[i].p[0]->point2D;
            vertex_1.color = (SDL_Color) {faces[i].n[0]->lightD, faces[i].n[0]->lightD, faces[i].n[0]->lightD};
            
            vertex_2.position = faces[i].p[1]->point2D;
            vertex_2.color = (SDL_Color) {faces[i].n[1]->lightD, faces[i].n[1]->lightD, faces[i].n[1]->lightD};
            
            vertex_3.position = faces[i].p[2]->point2D;
            vertex_3.color = (SDL_Color) {faces[i].n[2]->lightD, faces[i].n[2]->lightD, faces[i].n[2]->lightD};

            vertices[0] = vertex_1;
            vertices[1] = vertex_2;
            vertices[2] = vertex_3;

            SDL_RenderGeometry(renderer, NULL, vertices, 3, NULL, 0);
        }
    }
}

void RenderPolygonsTextured(SDL_Renderer *renderer, SDL_Surface *screenSurface, Polygon polygons[], float zbuffer[], Camera *camera){

    SDL_Rect rect = { 0, 0, screenSurface->w, screenSurface->h };
    SDL_FillRect(screenSurface, &rect, SDL_MapRGB(screenSurface->format, 0, 0, 0));
    for (int i=0; i<WIDTH*HEIGHT; i++){
        zbuffer[i] = 1000000000;  // TOFIX
    }

    int i, j;
    for (i=0; i<totalPolys; i++){
        for (j=0; j<polygons[i].size; j++){
            if (polygons[i].faces[j].show)
                drawTriangle(polygons[i].faces[j].p, polygons[i].faces[j].n, polygons[i].faces[j].t, screenSurface->pixels, polygons[i].surfaceT, &polygons[i].phongModel, zbuffer, camera);
        }
    }
}

void readFile(char *filePath, Point3D **points, SDL_FPoint **textures, VertexN **normals, Polygon **polygons, Camera *camera)
{
    FILE *fichier = fopen(filePath, "r");
    char* repertoire = dirname(strdup(filePath));
    chdir(repertoire);
    free(repertoire);

    if (fichier == NULL)
    {
        perror("Erreur lors de l'ouverture du fichier");
        return;
    }

    hasNorms = 0;
    hasTextures = 0;
    int totalFaces = 0;
    totalPoints = 0;
    totalNorms = 0;
    totalTextures = 0;

    char ligne[128];

    free(*points);
    free(*normals);
    free(*textures);
    for (int i=0; i<totalPolys; i++){
        free((*polygons)[i].faces);
    }
    free(*polygons);
    totalPolys = 0;

    char mtlFile[64];
    char polyName[64];
    Polygon *currentPoly = NULL;

    while (fgets(ligne, sizeof(ligne), fichier) != NULL)
    {
        if (strncmp(ligne, "mtllib", 6) == 0){
            sscanf(ligne, "mtllib %s", mtlFile);
            readMTL(mtlFile, polygons);
        }
        if (ligne[0] == 'v'){
            if (ligne[1]==' ')
                totalPoints++;
            else if (ligne[1]=='n'){
                totalNorms++;
                hasNorms = 1;
            } else if (ligne[1] == 't'){
                hasTextures = 1;
                totalTextures ++;
            }
        }else if (ligne[0] == 'f'){
            totalFaces++;
        }else if ((strncmp(ligne, "usemtl", 6) == 0)){
            sscanf(ligne, "usemtl %s", polyName);
            if (currentPoly == NULL){
                currentPoly = searchPoly(*polygons, polyName);
                continue;
            }
            currentPoly->size = totalFaces;
            currentPoly->faces = malloc(sizeof(Face) * totalFaces);
            if (!hasNorms) totalNorms += totalFaces;
            totalFaces = 0;
            currentPoly = searchPoly(*polygons, polyName);
            if (currentPoly == NULL){
                fprintf(stderr, "Unknown material in obj file: %s\n", ligne);
                return;
            }
        }
    }

    if (currentPoly != NULL){
        currentPoly->size = totalFaces;
        if (!hasNorms) totalNorms += totalFaces;
        currentPoly->faces = malloc(sizeof(Face) * totalFaces);
        currentPoly = searchPoly(*polygons, polyName);
        if (currentPoly == NULL){
            fprintf(stderr, "Unknown material in obj file: %s\n", ligne);
            return;
        }
    }

    Face *faces;

    *points = malloc(sizeof(Point3D) * totalPoints);
    *normals = malloc(sizeof(VertexN) * totalNorms);

    if (hasTextures)
        *textures = malloc(sizeof(SDL_FPoint) * totalTextures);
    if (totalPolys == 0){
        *polygons = malloc(sizeof(Polygon));
        (*polygons)[0] = (Polygon) {malloc(sizeof(Face) * totalFaces), NULL, totalFaces, ""};
        faces = (*polygons)[0].faces;
        totalPolys = 1;
    }

    fseek(fichier, 0, SEEK_SET);

    float x, y, z;
    int v[3];
    int vt[3];
    int vn[3];
    int countP = 0;
    int countF = 0;
    int countN = 0;
    int countT = 0;
    currentPoly = NULL;

    int linecount = 1;

    Vertex maxPoint, minPoint;
    float rnorm;

    Vertex light = {1/SDL_sqrt(3), 1/SDL_sqrt(3), 1/SDL_sqrt(3)};

    while (fgets(ligne, sizeof(ligne), fichier) != NULL)
    {
        if (ligne[0] == 'v')
        {
            switch (ligne[1])
            {
            case ' ':
                sscanf(ligne, "v %f %f %f", &x, &y, &z);

                updateMax(&maxPoint, x, y, z);
                updateMin(&minPoint, x, y, z);

                (*points)[countP].x = x;
                (*points)[countP].y = y;
                (*points)[countP].z = z;
                (*points)[countP].dz = z;
                countP++;
                break;
            
            case 'n':
                sscanf(ligne, "vn %f %f %f", &x, &y, &z);
                (*normals)[countN] = (VertexN) {x, y, z};

                countN++;
                break;

            case 't':
                sscanf(ligne, "vt %f %f", &x, &y);
                (*textures)[countT] = (SDL_FPoint) {1-x, 1-y};
                countT++;
                break;            
            }

        } else if (ligne[0] == 'f')
        {
            if (!hasTextures && !hasNorms)
                sscanf(ligne, "f %d %d %d", v, v+1, v+2);
            else if (hasTextures && !hasNorms)
                sscanf(ligne, "f %d/%d %d/%d %d/%d", v, vt, v+1, vt+1, v+2, vt+2);
            else
                sscanf(ligne, "f %d/%d/%d %d/%d/%d %d/%d/%d", v, vt, vn, v+1, vt+1, vn+1, v+2, vt+2, vn+2);

            faces[countF].show = 1;

            faces[countF].p[0] = *points + v[0] - 1;
            faces[countF].p[1] = *points + v[1] - 1;
            faces[countF].p[2] = *points + v[2] - 1;
            if (hasNorms){
                faces[countF].n[0] = *normals + vn[0] - 1;
                faces[countF].n[1] = *normals + vn[1] - 1;
                faces[countF].n[2] = *normals + vn[2] - 1;
                for (int i=0; i<3; i++){
                    rnorm = Q_rsqrt(SQUARE(faces[countF].n[0]->x) + SQUARE(faces[countF].n[0]->y) + SQUARE(faces[countF].n[0]->z));
                    faces[countF].n[0]->x *= rnorm;
                    faces[countF].n[0]->y *= rnorm;
                    faces[countF].n[0]->z *= rnorm;
                }
            }
            if (hasTextures){
                faces[countF].t[0] = *textures + vt[0] - 1;
                faces[countF].t[1] = *textures + vt[1] - 1;
                faces[countF].t[2] = *textures + vt[2] - 1;
            }
            countF++;
        } else if ((strncmp(ligne, "usemtl", 6) == 0)){
            countF = 0;
            sscanf(ligne, "usemtl %s", polyName);
            currentPoly = searchPoly(*polygons, polyName);
            faces = currentPoly->faces;
        }
        linecount++;
    }

    for (int i=0; i<totalPolys; i++)
        calculateNorms((*polygons)[i].faces, (*polygons)[i].size);

    countN = 0;
    if (!hasNorms){
        uniformLight = 1;
        for (int i=0; i<totalPolys; i++){
            for (int j=0; j<(*polygons)[i].size; j++){
                (*normals)[countN] = (VertexN) {(*polygons)[i].faces[j].normal.x,
                                                (*polygons)[i].faces[j].normal.y,
                                                (*polygons)[i].faces[j].normal.z};
                (*polygons)[i].faces[j].n[0] = *normals + countN;
                (*polygons)[i].faces[j].n[1] = *normals + countN;
                (*polygons)[i].faces[j].n[2] = *normals + countN;
                countN++;
            }
        }
    }

    calculateLightD(*normals, &light);

    camera->center = (Vertex) {(maxPoint.x + minPoint.x)/2, (maxPoint.y + minPoint.y)/2, (maxPoint.z + minPoint.z)/2};
    setCameraAngle(camera, 0, 0);
    if (!orbit){
        camera->x = 0;
        camera->y = 0;
        camera->z = 0;
    }
    fclose(fichier);
}

void readMTL(char file[], Polygon **polygons){
    FILE *fichier = fopen(file, "r");
    char ligne[128];
    totalPolys = 0;

    while (fgets(ligne, sizeof(ligne), fichier) != NULL){
        if (strncmp(ligne, "newmtl", 6) == 0)
            totalPolys += 1;
    }
    if (totalPolys == 0) return;

    *polygons = malloc(sizeof(Polygon) * totalPolys);

    fseek(fichier, 0, SEEK_SET);

    int polyCount = -1;
    char texturePath[64];
    SDL_Surface *imageSurface;
    float r,g,b;
    char colorType;

    while (fgets(ligne, sizeof(ligne), fichier) != NULL){
        if (strncmp(ligne, "newmtl", 6) == 0){
            polyCount ++;
            sscanf(ligne, "newmtl %s", (*polygons)[polyCount].name);
        } else if (strncmp(ligne, "map_Kd", 6) == 0){
            sscanf(ligne, "map_Kd %s", texturePath);
            if ((imageSurface = IMG_Load(texturePath)) == NULL){
                fprintf(stderr, "Wrong path name in mtl file %s\n", texturePath);
            }
            (*polygons)[polyCount].surfaceT = SDL_ConvertSurface(imageSurface, SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888), 0);
        } else if (ligne[0] == 'K'){
            sscanf(ligne, "K%c %f %f %f", &colorType, &r, &g, &b);
            if (colorType=='d') (*polygons)[polyCount].phongModel.diffuse = (ColorFloat) {r, g, b};
            else if (colorType=='a') (*polygons)[polyCount].phongModel.ambiant = (ColorFloat) {r, g, b};
            else if (colorType=='s') (*polygons)[polyCount].phongModel.specular = (ColorFloat) {r, g, b};
        }
    }
    fclose(fichier);
}

Polygon* searchPoly(Polygon polygons[], char name[]){
    for (int i=0; i<totalPolys; i++){
        if (strcmp(name, polygons[i].name) == 0)
            return polygons + i;
    }
    return NULL;
}

void calculateDistance(Face faces[], int size)
{
    for (int i=0; i<size; i++){
        faces[i].distance = (faces[i].p[0]->dz + faces[i].p[1]->dz + faces[i].p[2]->dz) / 3;
        // faces[i].distance = faces[i].p[0]->dz;
        // for (int k=1; k<3; k++){
        //     if (faces[i].distance > faces[i].p[k]->dz)
        //         faces[i].distance = faces[i].p[k]->dz;
        // }
    }
}

void calculateNorms(Face faces[], int size)
{
    float norm;
    Vertex u, v;
    for (int i=0; i<size; i++){
        u = (Vertex) {
            faces[i].p[0]->x - faces[i].p[1]->x,
            faces[i].p[0]->y - faces[i].p[1]->y,
            faces[i].p[0]->z - faces[i].p[1]->z
        };
        v = (Vertex) {
            faces[i].p[2]->x - faces[i].p[1]->x,
            faces[i].p[2]->y - faces[i].p[1]->y,
            faces[i].p[2]->z - faces[i].p[1]->z
        };

        faces[i].normal = (Vertex) {
            u.y*v.z - u.z*v.y,
            u.z*v.x - u.x*v.z,
            u.x*v.y - u.y*v.x
        };

        norm = Q_rsqrt(faces[i].normal.x*faces[i].normal.x + faces[i].normal.y*faces[i].normal.y + faces[i].normal.z*faces[i].normal.z);

        faces[i].normal.x *= norm;
        faces[i].normal.y *= norm;
        faces[i].normal.z *= norm;
    }
}

void calculateLightD(VertexN norms[], Vertex *light)
{
    float scal;
    for (int i=0; i<totalNorms; i++){

        scal = norms[i].x*light->x + norms[i].y*light->y + norms[i].z*light->z;
        if (uniformLight){
            if (scal < 0) scal = 0;
            if (scal > 1) scal = 1;
        }
        norms[i].lightD = scal;
    }
}

void calculateLightS(VertexN norms[], Vertex *light, Camera *camera)
{
    float scal;
    Vertex view = {camera->sy*camera->cx,
                            - camera->sx,
                   camera->cx*camera->cy};

    Vertex h = {light->x - view.x, light->y - view.y, light->z - view.z};
    float normH = Q_rsqrt(h.x*h.x + h.y*h.y + h.z*h.z);
    if (normH != 0){
        h.x *= normH;
        h.y *= normH;
        h.z *= normH;
    }

    for (int i=0; i<totalNorms; i++){

        scal = norms[i].x*h.x + norms[i].y*h.y + norms[i].z*h.z;
        if (uniformLight){
            if (scal < 0) scal = 0;
        }

        norms[i].lightS = scal;
    }
}

void drawTriangle(Point3D *triangle[], VertexN *normals[], SDL_FPoint *texture[], Uint32 *screenPixels, SDL_Surface *Tsurface, PhongModel *phongModel, float zbuffer[], Camera *camera){
    for (int i=0; i<3; i++){
        if ((triangle[i]->dz - near) * (triangle[(i+1)%3]->dz - near) < 0){
            Point3D point;
            SDL_FPoint textureP;
            VertexN normal;
            float alpha = (triangle[i]->dz-near)/(triangle[i]->dz-triangle[(i+1)%3]->dz);
            interpolatePoint(triangle[i], triangle[(i+1)%3], &point, alpha);
            if (Tsurface != NULL) interpolate2DPoint(texture[i], texture[(i+1)%3], &textureP, alpha);
            interpolateNormal(normals[i], normals[(i+1)%3], &normal, alpha);
            projectPoint(&point, camera);
            point.dz = near;
            Point3D *newTriangle1[3] = {triangle[i], &point, triangle[(i+2)%3]};
            VertexN *newNormal1[3] = {normals[i], &normal, normals[(i+2)%3]};
            SDL_FPoint *newTexture1[3] = {texture[i], &textureP, texture[(i+2)%3]};
            Point3D *newTriangle2[3] = {&point, triangle[(i+1)%3], triangle[(i+2)%3]};
            VertexN *newNormal2[3] = {&normal, normals[(i+1)%3], normals[(i+2)%3]};
            SDL_FPoint *newTexture2[3] = {&textureP, texture[(i+1)%3], texture[(i+2)%3]};
            drawTriangle(newTriangle1, newNormal1, newTexture1, screenPixels, Tsurface, phongModel, zbuffer, camera);
            drawTriangle(newTriangle2, newNormal2, newTexture2, screenPixels, Tsurface, phongModel, zbuffer, camera);
            return;
        }
    }

    if (triangle[0]->dz < near || triangle[1]->dz < near || triangle[2]->dz < near) return;
    
    int minX = WIDTH, maxX = 0, minY = HEIGHT, maxY = 0;

    for (int i=0; i<3; i++){
        if (triangle[i]->point2D.x < minX) minX = triangle[i]->point2D.x;
        if (triangle[i]->point2D.x > maxX) maxX = triangle[i]->point2D.x;

        if (triangle[i]->point2D.y < minY) minY = triangle[i]->point2D.y;
        if (triangle[i]->point2D.y > maxY) maxY = triangle[i]->point2D.y;
    }

    if (maxX >= WIDTH) maxX = WIDTH - 1;
    if (maxY >= HEIGHT) maxY = HEIGHT - 1;
    if (minX < 0) minX = 0;
    if (minY < 0) minY = 0;

    if (maxX-minX <= 0 || maxY-minY <= 0) return;

    float d = triangle[0]->point2D.x*(triangle[1]->point2D.y-triangle[2]->point2D.y) + triangle[1]->point2D.x*(triangle[2]->point2D.y-triangle[0]->point2D.y) + triangle[2]->point2D.x*(triangle[0]->point2D.y-triangle[1]->point2D.y);

    float matrix[3][3] = {
        {(triangle[1]->point2D.y-triangle[2]->point2D.y)/d, (triangle[2]->point2D.x-triangle[1]->point2D.x)/d, (triangle[1]->point2D.x*triangle[2]->point2D.y-triangle[2]->point2D.x*triangle[1]->point2D.y)/d},
        {(triangle[2]->point2D.y-triangle[0]->point2D.y)/d, (triangle[0]->point2D.x-triangle[2]->point2D.x)/d, (triangle[2]->point2D.x*triangle[0]->point2D.y-triangle[0]->point2D.x*triangle[2]->point2D.y)/d},
        {(triangle[0]->point2D.y-triangle[1]->point2D.y)/d, (triangle[1]->point2D.x-triangle[0]->point2D.x)/d, (triangle[0]->point2D.x*triangle[1]->point2D.y-triangle[1]->point2D.x*triangle[0]->point2D.y)/d}
    };
    float lambda[3];

    int x, y, k, tx, ty, w, h;
    int pixel_coord;
    float Tx, Ty;
    float new_zbuff;
    float shiftX, shiftY;
    Uint32 *tpixels;
    float light;

    if (Tsurface != NULL){
        w = Tsurface->w;
        h = Tsurface->h;
        shiftX = (matrix[0][1]*texture[0]->x/triangle[0]->dz + matrix[1][1]*texture[1]->x/triangle[1]->dz + matrix[2][1]*texture[2]->x/triangle[2]->dz) * w;
        shiftY = (matrix[0][1]*texture[0]->y/triangle[0]->dz + matrix[1][1]*texture[1]->y/triangle[1]->dz + matrix[2][1]*texture[2]->y/triangle[2]->dz) * h;
        tpixels = Tsurface->pixels;
    }
    float minDz = triangle[0]->dz;
    if (triangle[1]->dz < minDz) minDz = triangle[1]->dz;
    if (triangle[2]->dz < minDz) minDz = triangle[2]->dz;
    Uint8 alpha;
    char inside = 0;
    int r, g, b;
    float spec;
    float diffuse;

    if (uniformLight){
        spec = SDL_pow(normals[0]->lightS, 50);
        r = (int) ((normals[0]->lightD * phongModel->diffuse.r * Kd + phongModel->ambiant.r * Ka + phongModel->specular.r * spec * Ks) * 255);
        g = (int) ((normals[0]->lightD * phongModel->diffuse.g * Kd + phongModel->ambiant.g * Ka + phongModel->specular.g * spec * Ks) * 255);
        b = (int) ((normals[0]->lightD * phongModel->diffuse.b * Kd + phongModel->ambiant.b * Ka + phongModel->specular.b * spec * Ks) * 255);
        if (r > 255) r = 255;
        if (g > 255) g = 255;
        if (b > 255) b = 255;
    }


    for (x = minX; x <= maxX; x++)
    {
        for (k=0; k<3; k++){
            lambda[k] = x*matrix[k][0] + (minY-1)*matrix[k][1] + matrix[k][2];
        }
        inside = 0;
        for (y = minY; y <= maxY; y++)
        {
            for (k=0; k<3; k++){
                lambda[k] += matrix[k][1];
            }
            if (zbuffer[y*WIDTH + x] < minDz) continue;

            if (lambda[0] >= 0 && lambda[0] <= 1 && lambda[1] >= 0 && lambda[1] <= 1 && lambda[2] >= 0 && lambda[2] <= 1){

                new_zbuff = 1/(lambda[0]/triangle[0]->dz + lambda[1]/triangle[1]->dz + lambda[2]/triangle[2]->dz);

                if (new_zbuff > zbuffer[y*WIDTH + x]){
                    continue;
                }
                zbuffer[y*WIDTH + x] = new_zbuff;

                if (Tsurface != NULL){
                    if (!inside){
                        Tx = (lambda[0]*texture[0]->x/triangle[0]->dz + lambda[1]*texture[1]->x/triangle[1]->dz + lambda[2]*texture[2]->x/triangle[2]->dz) * w;
                        Ty = (lambda[0]*texture[0]->y/triangle[0]->dz + lambda[1]*texture[1]->y/triangle[1]->dz + lambda[2]*texture[2]->y/triangle[2]->dz) * h;
                        inside = 1;
                    } else {
                        Tx += shiftX;
                        Ty += shiftY;
                    }
                
                    tx = (int) (Tx*new_zbuff) % w;
                    ty = (int) (Ty*new_zbuff) % h;

                    if (tx < 0) tx += w;
                    if (ty < 0) ty += h;

                    pixel_coord = ty*w + tx;

                    screenPixels[y*WIDTH + x] = tpixels[pixel_coord];
                } else {
                    if (!uniformLight){
                        spec = normals[0]->lightS*lambda[0] + normals[1]->lightS*lambda[1] + normals[2]->lightS*lambda[2];
                        if (spec < 0) spec = 0;
                        if (spec > 1) spec = 1;
                        else spec = SDL_pow(spec, 100);
                        diffuse = normals[0]->lightD*lambda[0] + normals[1]->lightD*lambda[1] + normals[2]->lightD*lambda[2];
                        if (diffuse < 0) diffuse = 0;
                        r = (int) ((diffuse * phongModel->diffuse.r * Kd + phongModel->ambiant.r * Ka + phongModel->specular.r * spec * Ks) * 255);
                        g = (int) ((diffuse * phongModel->diffuse.g * Kd + phongModel->ambiant.g * Ka + phongModel->specular.g * spec * Ks) * 255);
                        b = (int) ((diffuse * phongModel->diffuse.b * Kd + phongModel->ambiant.b * Ka + phongModel->specular.b * spec * Ks) * 255);
                        if (r > 255) r = 255;
                        if (g > 255) g = 255;
                        if (b > 255) b = 255;
                    }
                    screenPixels[y*WIDTH + x] = SDL_MapRGB(pixelFormat, r, g, b);
                }
            } else if (inside) break;
        }
    }
}

void setVisibleFaces(Face faces[], int size, Camera *camera){
    float x = camera->cx * camera->sy;
    float y = - camera->sx;
    float z = camera->cx * camera->cy;

    for (int i=0; i<size; i++){
        faces[i].show = faces[i].normal.x*x + faces[i].normal.y*y + faces[i].normal.z*z > 0;
    }
}

void swap(Face *a, Face *b)
{
    Face temp = *a;
    *a = *b;
    *b = temp;
}

void heapify(Face faces[], int n, int i)
{
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && faces[left].distance > faces[largest].distance)
        largest = left;

    if (right < n && faces[right].distance > faces[largest].distance)
        largest = right;

    if (largest != i)
    {
        swap(&faces[i], &faces[largest]);
        heapify(faces, n, largest);
    }
}

void sortFaces(Face faces[], int n)
{
    calculateDistance(faces, n);
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(faces, n, i);

    for (int i = n - 1; i > 0; i--)
    {
        swap(&faces[0], &faces[i]);
        heapify(faces, i, 0);
    }
}

float Q_rsqrt(float number)
{
    long i;
    float x2, y;
    const float threehalfs = 1.5F;

    x2 = number * 0.5F;
    y = number;
    i = *(long *)&y;           // evil floating point bit level hacking
    i = 0x5f3759df - (i >> 1); // what the fuck?
    y = *(float *)&i;
    y = y * (threehalfs - (x2 * y * y)); // 1st iteration
    // y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

    return y;
}

void calculateCameraPos(Camera * camera){
    camera->cx = SDL_cos(camera->Rx);
    camera->sx = SDL_sin(camera->Rx);
    camera->cy = SDL_cos(camera->Ry);
    camera->sy = SDL_sin(camera->Ry);

    if (orbit){
        camera->x = - camera->sy*camera->cx * camera->distance + camera->center.x;
        camera->y =   camera->sx * camera->distance + camera->center.y;
        camera->z = - camera->cx*camera->cy * camera->distance + camera->center.z;
    }
}

void setCameraAngle(Camera *camera, float Rx, float Ry){
    camera->Rx = Rx;
    camera->Ry = Ry;
    calculateCameraPos(camera);
}

void addCameraAngle(Camera *camera, float Rx, float Ry){
    camera->Rx += Rx;
    camera->Ry += Ry;
    calculateCameraPos(camera);
}

float min(float a, float b){
    if (a > b)  return b;
    return a;
}

float max(float a, float b){
    if (a < b)  return b;
    return a;
}

void updateMin(Vertex *point, float x, float y, float z){
    point->x = min(point->x, x);
    point->y = min(point->y, y);
    point->z = min(point->z, z);
}

void updateMax(Vertex *point, float x, float y, float z){
    point->x = max(point->x, x);
    point->y = max(point->y, y);
    point->z = max(point->z, z);
}