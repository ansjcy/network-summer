#ifndef _RGBPIXMAP
#define _RGBPIXMAP
#include <GLUT/glut.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>

using namespace std;
typedef struct
{
   //unsigned short type;
   unsigned int size;
   unsigned short reserve1;
   unsigned short reserve2;
   unsigned int offset;
}head;
typedef struct
{
   unsigned int infoSize;
   unsigned int width;
   unsigned int height;
   unsigned short planes;
   unsigned short bitCount;
   unsigned int  compressionType;
   unsigned int  imageSize;
   unsigned int  XPelsPerMeter;
   unsigned int  YPelsPerMeter;
   unsigned int  colorNum;
   unsigned int  colorImportantNum;
}information;
typedef struct
{
   unsigned char rgbB; //
   unsigned char rgbG; //
   unsigned char rgbR; //
   unsigned char rgbRES;
}palette;
class BMPHead{
public:
   head thisHead;
   information thisInfo;
   palette thisPalette[256];
   
   void getHead(FILE* fp);
   void getInfo(FILE* fp);
   void getPalette(FILE* fp);
   void showHead();
   void showInfo();
   void showPalette();
};

class mRGB
{
public:
   unsigned char r,g,b,a;
};

class RGBpixmap
{
private:
   mRGB* pixel; // array of pixels
public:
   int nRows, nCols; // dimensions of the pixmap
   int readBMPFile(string fname); // read BMP file into this pixmap
   void setTexture(GLuint textureName);
};
#endif



