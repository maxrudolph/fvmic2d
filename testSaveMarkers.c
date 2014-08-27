#include "fdcode.h"
#include "io.h"

int main(){
  MarkerSet markerset;
  Marker markers[50];
  markerset.markers = &markers[0];

  /* compute offset to key parts of markers struct*/
  printf("offset to markers[0].VX = %d,%d\n",&(markers[0].VX), (markers[0]));

}
