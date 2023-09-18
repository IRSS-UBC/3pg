#ifndef HEADERFNAME_H
#define HEADERFNAME_H

//#ident "@(#)HeaderFname.hpp   1.1 98/04/24 12:12:35"

#ifdef WIN32
#define strncasecmp strnicmp
#endif

void makehdrfname(char *fname, char *hdrfname);
void makewldfname(char *fname, char *wldfname);

#endif
