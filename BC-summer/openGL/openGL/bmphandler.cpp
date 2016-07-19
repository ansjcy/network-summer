// RGBpixmap.cpp - routines to read a BMP file
#include "bmphandler.hpp"
void BMPHead::getHead(FILE* fp)
{
   fread(&thisHead,1,sizeof(head),fp);
}
void BMPHead::getInfo(FILE* fp)
{
   fread(&thisInfo,1,sizeof(information),fp);
}
void RGBpixmap :: setTexture(GLuint index)
{
   //bind the index on texture.. and retrieve it
   glBindTexture(GL_TEXTURE_2D, index);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,GL_NEAREST);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, nCols, nRows,0, GL_RGBA, GL_UNSIGNED_BYTE, pixel);
}
int RGBpixmap:: readBMPFile(string fname)
{
   BMPHead mybmphead;
   fstream stream;
   stream.open(fname.c_str(), ios::in|ios::binary);
   if(!stream){
       cout << " can't open file: " << fname << endl;
       return 0;
   }
   //is bmp or not
   short ch;
   FILE *fp = fopen(fname.c_str(), "rb");
   fread(&ch, 1, sizeof(short), fp);
   if(ch != 0x4d42){
       cout<<"the file is not a bmp file!"<<endl;
       return 0;
   }
   //get the information in this picture
   mybmphead.getHead(fp);
   mybmphead.getInfo(fp);
   
   int k, row, col, pad, totalInRow;
   unsigned long numCols = mybmphead.thisInfo.width;
   unsigned long numRows = mybmphead.thisInfo.height;
   
   totalInRow = ((3 * numCols + 3)/4) * 4;
   pad = totalInRow - 3 * numCols;
   nRows = numRows;
   nCols = numCols;
   pixel = new mRGB[nRows * nCols];
   if (!pixel)
       return 0;
   
   long count = 0;
   char dum;
   for(row = 0; row < nRows; row++)
   {
       for(col = 0; col < nCols; col++)
       {
           char r,g,b;
           stream.get(b);
           stream.get(g);
           stream.get(r);
           pixel[count].r = r;
           pixel[count].g = g;
           pixel[count].b = b;
           pixel[count++].a = 255;
       }
       for(k = 0; k < pad ; k++)
           stream >> dum;
   }
   stream.close();
   return 1;
}

