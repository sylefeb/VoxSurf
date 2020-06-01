#pragma once

typedef struct 
{
  uint   width;
  uint   height;
  uchar  depth;
  uchar *pixels;
} t_image_nfo;

t_image_nfo *ReadTGAFile(const char *filename);
