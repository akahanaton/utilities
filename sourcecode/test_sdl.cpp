#include <GL/gl.h>
#include <SDL/SDL.h>

void initGL()
{
	    SDL_Init(SDL_INIT_VIDEO);
		    SDL_SetVideoMode(600,300,16,SDL_OPENGL);

}
void destroyGL()
{
	    SDL_Quit();
}
void draw()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glBegin(GL_TRIANGLES);
    glColor3f(1,0,0);
    glVertex3f(0,0,0);
    glVertex3f(1,0,0);
    glVertex3f(0,1,0);
    glEnd();
	SDL_GL_SwapBuffers();

}
bool running=true;
void quit()
{
	    running=false;
}
void loop()
{
	SDL_Event event;
	while(running)
	{
		while(SDL_PollEvent(&event))
		 {
		    switch(event.type)
		    {
		       case SDL_QUIT:
		            quit();
		            break;
		    }
		 }
		draw();
		SDL_Delay(50);
	}

}
int main(int argc,char* argv[])
{
   initGL();
   loop();
   destroyGL();
   return 0;
}

