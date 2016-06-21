#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::ivec2;
using glm::vec2;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

struct Pixel
{
	int x;
	int y;
	float zinv;
	vec3 illumination;
	vec3 pos3d;
	int u;
	int v;
};

struct Vertex
{
	vec3 position;
	vec3 normal;
	vec3 reflectance;
};

const int SCREEN_WIDTH = 700;
const int SCREEN_HEIGHT = 700;
const int f = SCREEN_WIDTH/2;	//focal length
SDL_Surface* screen;
int t;
vec3 cameraPos( 0, 0, -2.001 );

// Light
vec3 lightPos(0,-0.5,-0.7);
vec3 lightPower = 16.1f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

vector<Triangle> triangles;

//Used in update
mat3 R;
float yaw = 0; // Yaw angle controlling camera rotation around y-axis

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

vec3 currentColor;
vec3 currentNormal;
vec3 currentReflectance;


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
vec3 Illumination( vec3 pos, vec3 reflectance, vec3 normal );
void VertexShader( const Vertex& v, Pixel& p );
void PixelShader( const Pixel& p );
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color );
void DrawPolygonEdges( const vector<vec3>& vertices );
void ComputePolygonRows( const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels );
void DrawPolygonRows( const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels );
void DrawPolygon( const vector<Vertex>& vertices );
void PolygonRowsTest();


int main( int argc, char* argv[] ) {

	LoadTestModel( triangles );
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
	
	//PolygonRowsTest();


	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;

}

void Update() {

	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	//cout << "Render time: " << dt << " ms." << endl;

		Uint8* keystate = SDL_GetKeyState( 0 );

	if( keystate[SDLK_UP] )
	{
	// Move camera forward
		cameraPos.y += 0.1;
	}
	if( keystate[SDLK_DOWN] )
	{
	// Move camera backward
		cameraPos.y -= 0.1;
	}
	if( keystate[SDLK_LEFT] )
	{
	// Move camera to the left
		cameraPos.x += 0.1;
	}
	if( keystate[SDLK_RIGHT] )
	{
	// Move camera to the right
		cameraPos.x -= 0.1;
	}

}

void Draw() {
	
	for( int y=0; y<SCREEN_HEIGHT; ++y )
		for( int x=0; x<SCREEN_WIDTH; ++x )
			depthBuffer[y][x] = -1;


	SDL_FillRect( screen, 0, 0 );
	
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	for( unsigned int i=0; i<triangles.size(); ++i ) {
		currentReflectance = triangles[i].color;
		currentNormal = triangles[i].normal;
		vector<Vertex> vertices(3);
		vertices[0].position = triangles[i].v0;
		vertices[0].normal = triangles[i].normal;
		vertices[0].reflectance = triangles[i].color;
		vertices[1].position = triangles[i].v1;
		vertices[1].normal = triangles[i].normal;
		vertices[1].reflectance = triangles[i].color;
		vertices[2].position = triangles[i].v2;
		vertices[2].normal = triangles[i].normal;
		vertices[2].reflectance = triangles[i].color;
		DrawPolygon( vertices );
	}

	if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

// Returns illumination
vec3 Illumination( vec3 pos, vec3 reflectance, vec3 normal ) {
	vec3 rhat = lightPos - pos;
	float r = pow(pow(rhat.x,2) + pow(rhat.y,2) + pow(rhat.z,2),0.5);
	vec3 nhat = normal;
	float dotProduct = dot(rhat,nhat);
	float A = 4 * M_PI * pow(r,2);
	vec3 B = lightPower/A;
	vec3 D = B * max(dotProduct,float(0));
	return reflectance * (D + indirectLightPowerPerArea);
}

void VertexShader( const Vertex& v, Pixel& p ) {

	vec3 pos = (v.position - cameraPos);//*R;
	if (pos.z != 0) p.zinv = 1/pos.z;
	else p.zinv = std::numeric_limits<float>::max();
	p.x = int(f * pos.x * p.zinv) + (SCREEN_WIDTH/2);
	p.y = int(f * pos.y * p.zinv) + (SCREEN_HEIGHT/2);
	p.pos3d = v.position;

	// Illumination
	//p.illumination = Illumination(v.position, v.reflectance, v.normal);
}

void PixelShader( const Pixel& p ) {
	int x = p.x;
	int y = p.y;
	if( p.zinv > depthBuffer[y][x]) {
		// Illumination
		vec3 illumination = Illumination(p.pos3d, currentReflectance, currentNormal);
		depthBuffer[y][x] = p.zinv;
		PutPixelSDL( screen, x, y, illumination );
	}
}


void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result ) {
	int N = result.size();
	vec2 step = vec2(b-a) / float(max(N-1,1));
	vec2 current( a );
	for( int i=0; i<N; ++i ) {
		result[i] = current;
		current += step;
	}
}

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result ) {
	
	int N = result.size();

	vec3 apos3d = a.pos3d*a.zinv;
	vec3 bpos3d = b.pos3d*b.zinv;

	float dx = (b.x-a.x) / float(max(N-1,1));
	float dy = (b.y-a.y) / float(max(N-1,1));
	float dzinv = (b.zinv-a.zinv) / float(max(N-1,1));
	vec3 dill = (b.illumination-a.illumination) / float(max(N-1,1));
	vec3 dpos = (bpos3d-apos3d) / float(max(N-1,1));
	float cx = a.x;
	float cy = a.y;
	float czinv = a.zinv;
	vec3 cill = a.illumination;
	vec3 cpos = apos3d;

	for( int i=0; i<N; ++i ) {
		result[i].x = (int)cx;
		result[i].y = (int)cy;
		result[i].zinv = czinv;
		result[i].illumination = cill;
		result[i].pos3d = cpos/czinv;
		cx += dx;
		cy += dy;
		czinv += dzinv;
		cill += dill;
		cpos += dpos;
	}
}

//Draws lines between vertices
void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color ) {

	ivec2 delta = glm::abs( a - b );
	int pixels = glm::max( delta.x, delta.y ) + 1;

	// Pixel positions along the line using interpolate
	vector<ivec2> line( pixels );
	Interpolate( a, b, line );

	for(unsigned int i=0; i < line.size(); i++){
  	ivec2 l = line[i];
		vec3 color(1,1,1);
		PutPixelSDL( screen, l.x, l.y, color );	
	}

}

/*
// Gives Figure 3, i.e outlines scene
void DrawPolygonEdges( const vector<vec3>& vertices ) {

	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<ivec2> projectedVertices( V );
	for( int i=0; i<V; ++i ) {
		VertexShader( vertices[i], projectedVertices[i] );
	}

	// Loop over all vertices and draw the edge from it to the next vertex:
	for( int i=0; i<V; ++i ) {
		int j = (i+1)%V; // The next vertex
		vec3 color( 1, 1, 1 );
		DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
	}

}*/

void ComputePolygonRows(const vector<Pixel>& vertexPixels,vector<Pixel>& leftPixels,vector<Pixel>& rightPixels ) {
	// 1. Find max and min y-value of the polygon
	//    and compute the number of rows it occupies.

	int ymax = -numeric_limits<int>::max();
	int ymin = +numeric_limits<int>::max();

	for (unsigned int i = 0; i < vertexPixels.size(); i++){

		if (vertexPixels[i].y < ymin )
			ymin = vertexPixels[i].y;
		
		if (vertexPixels[i].y > ymax )
			ymax = vertexPixels[i].y;

	}
	
	int rows = (ymax - ymin) + 1;
	//cout << "rows  =" << rows << "\n";
	//cout << "ymin =" << ymin << "\n";
	//cout << "ymax =" << ymax << "\n";

	// 2. Resize leftPixels and rightPixels
	//    so that they have an element for each row.
	leftPixels.resize(rows);
	rightPixels.resize(rows);

	int y = ymin;

	//Set y values 
	for ( int i=0; i<rows; ++i ) {
		leftPixels[i].y  = y;
		rightPixels[i].y = y;
		y++;
	}


	// 3. Initialize the x-coordinates in leftPixels
	//    to some really large value and the x-coordinates
	//    in rightPixels to some really small value.
	for( int i=0; i<rows; ++i ) {
		leftPixels[i].x  = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
	}

	// 4. Loop through all edges of the polygon and use
	//    linear interpolation to find the x-coordinate for
	//    each row it occupies. Update the corresponding
	//    values in rightPixels and leftPixels.

	int V = vertexPixels.size();
	for( int i=0; i<V; i++ ) {
		int j = (i+1)%V; // The next vertex
		Pixel a = vertexPixels[i];
		Pixel b = vertexPixels[j];
		
		int dx = abs(a.x-b.x);
		int dy = abs(a.y-b.y);			
		int pixels = glm::max( dx, dy ) + 1;

		// Pixel positions along the line using interpolate
		vector<Pixel> line( pixels );
		Interpolate( a, b, line );
		
		for(unsigned int i=0; i < line.size()-1; i++){
			Pixel l = line[i];
			//cout << "(" << l.x << "," << l.y << ")\n"; 

			if (l.y < ymin) {
				cout << "(" << l.x << "," << l.y << ")" << " l.y less than ymin\n";
				cout << "out by: " << (l.y-ymin);
			} 
			else {
				if (l.x < leftPixels[l.y-ymin].x) {
					leftPixels[l.y-ymin].x = l.x;
					leftPixels[l.y-ymin].zinv = l.zinv;
					leftPixels[l.y-ymin].illumination = l.illumination;
					leftPixels[l.y-ymin].pos3d = l.pos3d;
				}
				if (l.x > rightPixels[l.y-ymin].x) {
					rightPixels[l.y-ymin].x = l.x;
					rightPixels[l.y-ymin].zinv = l.zinv;
					rightPixels[l.y-ymin].illumination = l.illumination;
					rightPixels[l.y-ymin].pos3d = l.pos3d;
				}
			}
		}
	}
}

void DrawPolygonRows( vector<Pixel>& leftPixels, vector<Pixel>& rightPixels ) {

	int rows = leftPixels.size();

	for ( int i=0; i<rows; ++i ) {

		// Interpolate rows
		vector<Pixel> line( rightPixels[i].x - leftPixels[i].x + 1);
		Interpolate( leftPixels[i], rightPixels[i], line );

		// Draw each pixel on the row
		for (unsigned int j = 0; j < line.size(); j++) {
			PixelShader(line[j]);
		}
	}
}


void DrawPolygon( const vector<Vertex>& vertices ) {
	
	int V = vertices.size();
	vector<Pixel> vertexPixels( V );
	
	for( int i=0; i<V; ++i ) 
		VertexShader( vertices[i], vertexPixels[i] );
	

	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
	DrawPolygonRows( leftPixels, rightPixels );
}


/*
void PolygonRowsTest(){

	vector<ivec2> vertexPixels(3);
	vertexPixels[0] = ivec2(10, 5);
	vertexPixels[1] = ivec2( 5,10);
	vertexPixels[2] = ivec2(15,15);
	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
	ComputePolygonRows( vertexPixels, leftPixels, rightPixels );

	for(unsigned int row=0; row<leftPixels.size(); ++row ) {
		cout << "Start: ("
		<< leftPixels[row].x << ","
		<< leftPixels[row].y << "). "
		<< "End: ("
		<< rightPixels[row].x << ","
		<< rightPixels[row].y << "). " << endl;
	}

}
*/





