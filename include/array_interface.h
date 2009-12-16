#include "Python.h"

#define PYCOGENT_VERSION "1.4"

/* Array Interface flags */
#define CONTIGUOUS    0x001
#define FORTRAN       0x002
#define ALIGNED       0x100
#define NOTSWAPPED    0x200
#define WRITEABLE     0x400

typedef struct PyArrayInterface {
  int version;
  int nd;
  char typekind;
  int itemsize;
  int flags;
  Py_intptr_t *shape;
  Py_intptr_t *strides;
  void *data;
} PyArrayInterface;

