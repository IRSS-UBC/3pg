#ifndef EXCEPTION_H
#define EXCEPTION_H

//#ident "@(#)Exception.hpp	1.2 98/07/15 13:47:14 JCG"


#include <string.h>


class Exception {
protected:
  char *msg;
  int fatal;

  void NoMemory(char *cp);
  
public:
  static char msgbuf[1000];

  Exception(char *cp = msgbuf, int isfatal = 0) : fatal(isfatal) {
    msg = new char[strlen(cp)+1];
    if (msg == 0)
      NoMemory(cp);
    else
      strcpy(msg, cp);
  }

  virtual ~Exception() {
    delete[] msg;
  }

  char *Message() {
    return msg;
  }

  int Fatal() {
    return fatal;
  }

};


class NoMemException : public Exception {
public:
  NoMemException(char *cp = 0);
};


#endif
