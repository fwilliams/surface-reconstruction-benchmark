/*
Benedict Brown
Princeton University

ppmio.c
Read a binary ppm file.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

float *read_ppm(char *fname, int *w, int *h) {
  FILE *f;

  unsigned char *data;
  float *img;
  char tmpstr[3];

  int i;

  f = fopen(fname, "r");

  // read the PPM header
  fscanf(f, " %3s ", tmpstr);

  // make sure this is really a binary PPM image
  if (strlen(tmpstr) != 2 || strcmp(tmpstr, "P6"))
    return NULL;

  // now eat any comments there may be
  while ((i = getc(f)) == '#')
    while(getc(f) != '\n');
  ungetc(i, f);

  // now read the image width and height
  fscanf(f, " %d %d %*d", w, h);
  getc(f);

  // now read the data
  data = (unsigned char *) malloc((*w) * (*h)  * 3);
  fread(data, 3, (*w) * (*h), f);
  fclose(f);

  img = malloc((*w) * (*h) * 3 * sizeof(float));
  for (i = 0; i < 3 * (*w) * (*h); i++)
    img[i] = (float) data[i] / 255.0f;
  free(data);

  return img;
}
